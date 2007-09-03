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
  NumTestIndividuals(Data->getNumberOfTestIndividuals()){

  SetNullValues();
  //NumCompLoci = Loci->GetNumberOfCompositeLoci();
  GlobalConcordanceCounts = 0;
  GlobalSumArrivalCounts = 0;
  ConcordanceCounts = new int[NumCompLoci*2*(options->getNumberOfBlockStates())];
  SumArrivalCounts = new int[NumCompLoci*options->getNumberOfBlockStates()];

  size = NumInd;
  
  //  Individual::SetStaticMembers(Loci, options);
  Individual::SetStaticMembers(Loci, options);

  _child = new Individual*[size];
  for (unsigned int i = 0; i < size; i++) {
    // _child[i] = new Individual(i+1, options, Data);//NB: first arg sets Individual's number
    //NB: first arg sets Individual's number
    HapMixChild.push_back(new HapMixIndividual(i+1, options, Data, theta));

    _child[i] = HapMixChild[i];
    if(!_child[i]->isHaploidIndividual())
      ++NumDiploidIndividuals;
  }
  
  if(options->OutputCGProbs())GPO.Initialise(Data->getNumberOfTestIndividuals(), Loci->GetNumberOfCompositeLoci()-Data->getNumTypedLoci());
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

void HapMixIndividualCollection::SampleHiddenStates(const HapMixOptions& options, unsigned iteration){

  fill(ConcordanceCounts, ConcordanceCounts + NumCompLoci*(options.getPopulations()+1), 0);
  fill(SumArrivalCounts, SumArrivalCounts + NumCompLoci*options.getPopulations(), 0);

  for(unsigned int i = 0; i < size; i++ ){

    /*do calculating of ordered probs first because it requires both forward and backward HMM updates
      and these can be done in parallel on a dual-core processor */
    if (
	(//If it's after the burnin 
	 (int)iteration > options.getBurnIn()
	 // and it's a test individual
	 && isTestIndividual(i)
	 // and if the score tests are switched on
	 && options.getTestForAllelicAssociation()
	 // and the individual is diploid
	 && (! _child[i]->isHaploidIndividual())
	 )
	//or it's a masked individual
	// maskedIndividuals indices are 1-based, offset of 1 is needed
	|| (options.OutputCGProbs() && isTestIndividual(i))
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
  if(NumTestIndividuals > 0)return NumTestIndividuals;
  else return size;
}

unsigned int HapMixIndividualCollection::getFirstScoreTestIndividualNumber()const{
  if(NumTestIndividuals > 0)return size - NumTestIndividuals;
  else return 0;
}

///determines if individual i is a case/control ie its genotype came from ccgenotypesfile
bool HapMixIndividualCollection::isTestIndividual(unsigned i)const{
  return ( i >= (size - NumTestIndividuals) );
}

void HapMixIndividualCollection::AccumulateConditionalGenotypeProbs(const HapMixOptions& options, 
								    const InputHapMixData &Data, unsigned NumCompositeLoci){
  
  unsigned j = 0;
  for(unsigned locus_i = 0; locus_i < NumCompositeLoci; ++locus_i) {
    if (!Data.isTypedLocus(locus_i)){
      unsigned i = 0;
      for(unsigned indiv_i = size - NumTestIndividuals; indiv_i < size; ++i, ++indiv_i) {
	GPO.Update(i, j, HapMixChild[indiv_i]->getUnorderedProbs(locus_i));
      }
      ++j;
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
