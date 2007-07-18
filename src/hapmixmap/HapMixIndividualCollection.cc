/* 
 *   HAPMIXMAP
 *   HapMixIndividualCollection.cc
 *   Container class for HapMixIndividual objects
 *   Copyright (c) 2006, 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "HapMixIndividualCollection.h"
#include "HapMixOptions.h"
#include "InputHapMixData.h"
//#include "HapMixIndividual.h"
#include "HapMixFreqs.h"
#include "bclib/Regression.h"

HapMixIndividualCollection
::HapMixIndividualCollection(const HapMixOptions* const options, 
			     InputHapMixData* const Data, Genome* Loci, const double* theta):
  IndividualCollection(Data->getNumberOfIndividuals(), options->getPopulations(), Loci->GetNumberOfCompositeLoci()),
  NumCaseControls(Data->getNumberOfCaseControlIndividuals()){

  SetNullValues();
  //NumCompLoci = Loci->GetNumberOfCompositeLoci();
  GlobalConcordanceCounts = 0;
  GlobalSumArrivalCounts = 0;
  ConcordanceCounts = new int[NumCompLoci*2*(options->getNumberOfBlockStates())];
  SumArrivalCounts = new int[NumCompLoci*options->getNumberOfBlockStates()];

  //Populations = options->getPopulations();
  //NumInd = Data->getNumberOfIndividuals();//number of individuals, including case-controls
  size = NumInd;
  //  NumCaseControls = Data->getNumberOfCaseControlIndividuals();
  
  //  Individual::SetStaticMembers(Loci, options);
  Individual::SetStaticMembers(Loci, options);

  _child = new Individual*[size];
  for (unsigned int i = 0; i < size; i++) {
    // _child[i] = new Individual(i+1, options, Data);//NB: first arg sets Individual's number
    //NB: first arg sets Individual's number
    HapMixChild.push_back(new HapMixIndividual(i+1, options, Data, theta, isMaskedIndividual(i, options->getMaskedIndividuals())));

    _child[i] = HapMixChild[i];
    if(!_child[i]->isHaploidIndividual())
      ++NumDiploidIndividuals;
  }
  
  if(options->OutputCGProbs())GPO.Initialise(options->GetNumMaskedIndividuals(), options->GetNumMaskedLoci());
}

HapMixIndividualCollection::~HapMixIndividualCollection(){
  delete[] ConcordanceCounts;
}
const HapMixIndividual* HapMixIndividualCollection::getHapMixIndividual(int num)const{
  if (num < (int)size){
    return HapMixChild[num];
  } else {
    throw string("ERROR in HMIC::getHapMixIndividual: index out of range");
  }
}

//indicates if an individual's genotypes have been masked
bool HapMixIndividualCollection::isMaskedIndividual(unsigned i, const vector<unsigned>& maskedIndividuals)const{
  const vector<unsigned>::const_iterator mi_end = maskedIndividuals.end();
  return (find(maskedIndividuals.begin(), mi_end, i+1) != mi_end);

}
void HapMixIndividualCollection::SampleHiddenStates(const HapMixOptions& options, unsigned iteration){

  fill(ConcordanceCounts, ConcordanceCounts + NumCompLoci*(options.getPopulations()+1), 0);
  fill(SumArrivalCounts, SumArrivalCounts + NumCompLoci*options.getPopulations(), 0);

  const vector<unsigned>& maskedIndividuals = options.getMaskedIndividuals();
  const vector<unsigned>::const_iterator mi_begin = maskedIndividuals.begin();
  const vector<unsigned>::const_iterator mi_end = maskedIndividuals.end();

  for(unsigned int i = 0; i < size; i++ ){

    /*do calculating of ordered probs first because it requires both forward and backward HMM updates
      and these can be done in parallel on a dual-core processor */
    if (
	(//If it's after the burnin 
	 (int)iteration > options.getBurnIn()
	 // and it's a case or control individual
	 && isCaseControl(i)
	 // and if the score tests are switched on
	 && options.getTestForAllelicAssociation()
	 // and the individual is diploid
	 && (! _child[i]->isHaploidIndividual())
	 )
	//or it's a masked individual
	// maskedIndividuals indices are 1-based, offset of 1 is needed
	|| (find(mi_begin, mi_end, i+1) != mi_end)
	)
    {
     HapMixChild[i]->calculateUnorderedGenotypeProbs(options);
    } 
    // ** Run HMM forward recursions and sample hidden states
    _child[i]->SampleHiddenStates(options);

    //calculate and store loglikelihoods, ready for energy accumulation
    _child[i]->getLogLikelihood(options, false, true);

    //accumulate sufficient statistics for update of arrival rates and mixture proportions
    HapMixChild[i]->AccumulateConcordanceCounts(ConcordanceCounts);
    HapMixChild[i]->SampleJumpIndicators(SumArrivalCounts);

    //set HMMisBad before moving to the next individual, but keep stored loglikelihood
    _child[i]->HMMIsBad(false); 
  }
}

const int* HapMixIndividualCollection::getConcordanceCounts()const{
  return ConcordanceCounts;
}

const int* HapMixIndividualCollection::getSumArrivalCounts()const{
  return SumArrivalCounts;
}

int HapMixIndividualCollection::getNumberOfIndividualsForScoreTests()const{
  if(NumCaseControls > 0)return NumCaseControls;
  else return size;
}

unsigned int HapMixIndividualCollection::getFirstScoreTestIndividualNumber()const{
  if(NumCaseControls > 0)return size - NumCaseControls;
  else return 0;
}

///determines if individual i is a case/control ie its genotype came from ccgenotypesfile
bool HapMixIndividualCollection::isCaseControl(unsigned i)const{
  return ( i > (size - NumCaseControls) );
}

void HapMixIndividualCollection::AccumulateConditionalGenotypeProbs(const HapMixOptions& options, const Genome& Loci){
  
  const std::vector<unsigned>& MaskedLoci = options.getMaskedLoci();
  const std::vector<unsigned>& MaskedIndividuals = options.getMaskedIndividuals();
  // vu_ci stands for vector<unsigned>::const_iterator
  typedef std::vector<unsigned>::const_iterator vu_ci;
  unsigned j = 0;
  //NB: indices count from 1 so must be offset by -1
  for(vu_ci locus_i = MaskedLoci.begin(); locus_i!= MaskedLoci.end(); ++j, ++locus_i) {
    if (*locus_i <= Loci.GetNumberOfCompositeLoci()){
      unsigned i = 0;
      for(vu_ci indiv_i = MaskedIndividuals.begin(); indiv_i!= MaskedIndividuals.end(); ++i, ++indiv_i) {
        if(*indiv_i <= size) {
          GPO.Update(i, j, HapMixChild[(*indiv_i) - 1]->getUnorderedProbs((*locus_i) - 1));
        }
      }
    }
  }
}
void HapMixIndividualCollection::OutputCGProbs(const char* filename, const Vector_s& MaskedLocusLabels){
  GPO.Output(filename, MaskedLocusLabels);
}

//TODO: use posterior means of mixture props if fixedmixtureprops=0
// double HapMixIndividualCollection
// ::getDevianceAtPosteriorMean(const Options& options, vector<bclib::Regression *> &R, Genome* Loci,
//                              LogWriter &Log, const double* const MixtureProps, const vector<double>& SumLogLambda){

//   //TODO: broadcast SumLogRho to workers
//   //TODO: set mixture props to posterior means

//   //SumRho = ergodic sum of global sumintensities
//   const int iterations = options.getTotalSamples()-options.getBurnIn();
//   int NumDiploid = 0;
  
//   NumDiploid = getNumDiploidIndividuals();

//   //update chromosomes using globalrho, for globalrho model
//   vector<double> LambdaBar(Loci->GetNumberOfCompositeLoci());

//   for(unsigned i = 0; i < Loci->GetNumberOfCompositeLoci(); ++i)
//     LambdaBar[i] = (exp(SumLogLambda[i] / (double)iterations));

//   //set locus correlation

//   Loci->SetLocusCorrelation(LambdaBar);

//   for( unsigned int j = 0; j < Loci->GetNumberOfChromosomes(); j++ )
//     //set global state arrival probs in hapmixmodel
//     //TODO: can skip this if xonly analysis with no females
//     //NB: assumes always diploid in hapmixmodel
//     Loci->getChromosome(j)->HMM->SetStateArrivalProbs(MixtureProps, options.isRandomMatingModel(),
// 						      (NumDiploid)>0);
  
//   //set haplotype pair probs to posterior means 
//   if(NumDiploid)
//     for( unsigned int j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ )
//       (*Loci)(j)->SetHapPairProbsToPosteriorMeans(iterations);
  
//   //set genotype probs using happair probs calculated at posterior means of allele freqs 
//   //setGenotypeProbs(Loci, A);

//   //accumulate deviance at posterior means for each individual
//   double Lhat = 0.0; // Lhat = loglikelihood at estimates
//   for(unsigned int i = 0; i < size; i++ ){
//     Lhat += _child[i]->getLogLikelihood(options, false, false);
//   }
  
//   Log << Quiet << "DevianceAtPosteriorMean(parameters)" << -2.0*Lhat << "\n";
//   for(unsigned c = 0; c < R.size(); ++c){
//     double RegressionLogL = R[c]->getLogLikelihoodAtPosteriorMeans(iterations, getOutcome(c));
//     Lhat += RegressionLogL;
//     Log << "DevianceAtPosteriorMean(Regression " << c+1 << ")"
// 	<< -2.0*RegressionLogL << "\n";
//   }
//   return(-2.0*Lhat);
// }
