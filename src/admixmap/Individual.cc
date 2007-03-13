/** 
 *   Individual.cc 
 *   Class to represent an individual and update individual-level parameters
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "Individual.h"
#include "regression/Regression.h" 
//#include "utils/misc.h"
//#include "utils/dist.h"
//#include "utils/linalg.h"
//#include <algorithm>
//#include <limits>
#include <sstream>
#include <exception>

#define PR(x) cout << #x << " = " << x << endl;

// ******* Static Member Declarations
bool Individual::Xdata;
unsigned Individual::X_posn;
unsigned int Individual::numChromosomes;
Genome *Individual::Loci;
int Individual::Populations;

//******** Constructors **********
Individual::Individual() {
  // initialise pointers here
  PossibleHapPairs = NULL;
  GenotypesMissing = NULL;
  missingGenotypes = NULL;
  LocusAncestry = NULL;
  Theta = NULL;
}

// Individual::Individual(int number, const Options* const options, const InputData* const Data)
// : myNumber(number)
// {
//   missingGenotypes = NULL;//allocated later, if needed
//   Initialise(number, options, Data);
// }

void Individual::Initialise(int number, const Options* const options, const InputData* const Data){
  myNumber = number;
  if( options->isRandomMatingModel() && !isHaploid) NumGametes = 2;
  else NumGametes = 1;
  
  // Read sex value if present.
  SexIsFemale = false;
  if (options->getgenotypesSexColumn() == 1) {//if sex column in genotypesfile
    SexIsFemale = Data->isFemale(myNumber);
  }
  
  double L = Loci->GetLengthOfGenome();
  double LX = 0.0;
  if(Xdata) LX = Loci->GetLengthOfXchrm();
  // effective length of genome is L + 0.5*LX if there is an X chrm: i.e. if g=1 or sex is female
  if(SexIsFemale) {
    EffectiveL[0] = L + 0.5*LX;
  } else {
    EffectiveL[0] = L;
  }
  EffectiveL[1] = L + 0.5*LX;
  
  int numCompositeLoci = Loci->GetNumberOfCompositeLoci();
  
  // vector of possible haplotype pairs - 2 integers per locus if diploid, 1 if haploid
  PossibleHapPairs = new vector<hapPair>[numCompositeLoci];
  
  LocusAncestry = new int*[ numChromosomes ]; // array of matrices in which each col stores 2 integers   
  for (unsigned chrNo = 0; chrNo < numChromosomes; ++chrNo) {
    LocusAncestry[chrNo] = NULL;
  }
  //initialise genotype probs array and array of indicators for genotypes missing at locus
  
  size_t AncestrySize = 0;  // set size of locus ancestry array
  //gametes holds the number of gametes for each chromosome, either 1 or 2
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    if(isHaploid || (!SexIsFemale && Loci->isXChromosome(j))){//haploid on this chromosome
      AncestrySize = Loci->GetSizeOfChromosome(j) ;
      gametes.push_back(1);
    }
    else{
      AncestrySize = 2 * Loci->GetSizeOfChromosome(j) ;
      gametes.push_back(2);
    }
    LocusAncestry[j] = new int[ AncestrySize];
    for(unsigned i = 0; i < AncestrySize; ++i) LocusAncestry[j][i] = 0;
  }
  
  logLikelihood.value = 0.0;
  logLikelihood.ready = false;
  logLikelihood.HMMisOK = false;
  
  // Allocate space for unordered genotype probs
  // They have form of vector of vectors of vectors of doubles.
  if(options->getHapMixModelIndicator() && options->getTestForAllelicAssociation()){
    vector<double> v1 = vector<double>(1);
    vector<vector<double> > v3 = vector<vector<double> >(3, v1);
    UnorderedProbs = vector<vector<vector<double> > >(numCompositeLoci, v3);
  }
}

//********** Destructor **********
Individual::~Individual() {
  if (PossibleHapPairs) delete[] PossibleHapPairs;

  if(GenotypesMissing){
    for(unsigned j = 0; j < numChromosomes; ++j) {
      delete[] GenotypesMissing[j];
    }
    delete[] GenotypesMissing;
  }
  if(LocusAncestry){
    for(unsigned j = 0; j < numChromosomes; ++j) {
      delete[] LocusAncestry[j];
    }
    delete[] LocusAncestry;
  }

  delete[] missingGenotypes;
}

void Individual::SetUniformAdmixtureProps() {
  size_t K = Populations;
  for( unsigned g = 0; g < NumGametes; ++g ) { 
    for(size_t k = 0; k < K; ++k)
      Theta[g*K+k] = 1.0 / K;
  }
}

void Individual::setOutcome(double* Y){
  Outcome = Y;
}
void Individual::setCovariates(double* X){
  Covariates = X;
}

void Individual::setGenotypesToMissing(){
  for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c)
    for(unsigned j = 0; j < Loci->GetSizeOfChromosome(c); ++j)
      GenotypesMissing[c][j] = true;
}

void Individual::DeleteGenotypes(){
  unsigned noCompositeLoci = Loci->GetNumberOfCompositeLoci();
  for(unsigned j = 0; j < noCompositeLoci; ++j){
#ifdef PARALLEL
      genotypes[j][0].clear();
#else
    int noLoci = Loci->getNumberOfLoci(j);
    for(int k = 0; k < noLoci; ++k)
      genotypes[j][k].clear();
#endif
    genotypes[j].clear();
  }
  genotypes.clear();
}

/// sets static members, including allocation and deletion of static objects for score tests
void Individual::SetStaticMembers(Genome* const pLoci, const Options* const options){
  Loci = pLoci;
  numChromosomes = Loci->GetNumberOfChromosomes();
  Populations = options->getPopulations();
  Xdata = Loci->isX_data();
  X_posn = 9999; //position of the X chromosome in the sequence of chromosomes in the input data
  if(Xdata) X_posn = Loci->GetChrNumOfLocus(Loci->getFirstXLocus());//too clunky, should simplify
}

/// Sets genome (Loci)
/// Needed for tests which need to substitute only this one static
/// variable.
// void Individual::setGenome(Genome *pLoci)
// {
//   Loci = pLoci;
// }

// void Individual::setPopulations(int p)
// {
//   Populations = p;
// }

void Individual::HMMIsBad(bool loglikisbad) {
  logLikelihood.HMMisOK = false;
  if(loglikisbad)logLikelihood.ready = false;
}

//******************** Accessors ***********************************************************
const double* Individual::getAdmixtureProps()const {
  return Theta;
}

const std::vector<hapPair > &Individual::getPossibleHapPairs(unsigned int locus)const{
  return PossibleHapPairs[locus];
}

const int* Individual::getSampledHapPair(int locus)const{
  return sampledHapPairs[locus].haps;
}

/** 
    Get the locus ancestry by the absolute locus index,ingoring the chromosome.
 */
void Individual::GetLocusAncestry(int locus, int Ancestry[2])const{
  unsigned c, l;
  Loci->GetChrmAndLocus(locus, &c, &l);
  GetLocusAncestry(c, l, Ancestry);
}

void Individual::GetLocusAncestry(int chrm, int locus, int Ancestry[2])const{
  Ancestry[0]  = LocusAncestry[chrm][locus];
  if(isHaploid || ((unsigned)chrm == X_posn && !SexIsFemale))
    Ancestry[1] = Ancestry[0];
  else
    Ancestry[1] = LocusAncestry[chrm][Loci->GetSizesOfChromosomes()[chrm]  + locus];
}

///returns value of LocusAncestry at a locus for a particular gamete
int Individual::GetLocusAncestry(int chrm, int gamete, int locus)const{
  int g = (gametes[chrm] == 2) ? gamete : 0; //so that gamete = 1 works when gametes[chrm] = 1;
  return LocusAncestry[chrm][g * Loci->GetSizesOfChromosomes()[chrm]  + locus] ;
}

///Indicates whether genotype is missing at all simple loci within a composite locus
bool Individual::GenotypeIsMissing(unsigned int locus)const {
  unsigned c, l;
  Loci->GetChrmAndLocus(locus, &c, &l);
  return GenotypesMissing[c][l];
}
///Indicates whether genotype is missing at a simple locus
//used by HW score test
bool Individual::simpleGenotypeIsMissing(unsigned locus)const{
  if(!missingGenotypes)throw string("missingGenotypes not allocated");
  return missingGenotypes[locus];
}

bool Individual::isHaploidatLocus(unsigned j)const{
  return (bool)(isHaploid || (!SexIsFemale && Loci->isXLocus(j)));
}
bool Individual::isHaploidIndividual()const{
  return isHaploid;
}
//****************** Log-Likelihoods **********************
// public function: 
// calls private function to get log-likelihood at current parameter values, and stores it either as loglikelihood.value or as loglikelihood.tempvalue
// store should be false when calculating energy for an annealed run, or when evaluating proposal for global sum-intensities
double Individual::getLogLikelihood( const Options* const options, const bool forceUpdate, const bool store) {

  if (!logLikelihood.ready || forceUpdate) {
    logLikelihood.tempvalue = getLogLikelihood(options, Theta, _rho, true);
    if(store) {  
      logLikelihood.value = logLikelihood.tempvalue; 
      logLikelihood.ready = false; //true;
      logLikelihood.HMMisOK = false; //true; //because forward probs now correspond to current parameter values 
    }                               //and call to UpdateHMM has set this to false
    return logLikelihood.tempvalue;   
  } else return logLikelihood.value; // nothing was changed
}

// private function: gets log-likelihood at parameter values specified as arguments, but does not update loglikelihoodstruct
double Individual::getLogLikelihood(const Options* const options, const double* const theta, 
				    const vector<double > rho,  bool updateHMM) {
  double LogLikelihood = 0.0;
  for( unsigned int j = 0; j < numChromosomes; j++ ) {
    if(updateHMM){// force update of forward probs 
      UpdateHMMInputs(j, options, theta, rho);
    }
    LogLikelihood += Loci->getChromosome(j)->getLogLikelihood( !isHaploid && (!Loci->isXChromosome(j) || SexIsFemale) );
  }
  return LogLikelihood; // if HMM update not required, can just use stored log-likelihood  
}

void Individual::storeLogLikelihood(const bool setHMMAsOK) { // to call if a Metropolis proposal is accepted
    logLikelihood.value = logLikelihood.tempvalue; 
    logLikelihood.ready = true;
    if(setHMMAsOK) logLikelihood.HMMisOK = true; 
}                               

double Individual::getLogLikelihoodAtPosteriorMeans(const Options* const options) {
  // should set allele freqs also to posterior means, and recalculate prob genotypes at these freqs before calling getloglikelihood 
  double LogLikelihood = 0.0;
  for( unsigned int j = 0; j < numChromosomes; j++ ) {
    UpdateHMMInputs(j, options, Theta, _rho); 
    LogLikelihood += Loci->getChromosome(j)->getLogLikelihood( !isHaploid && (!Loci->isXChromosome(j) || SexIsFemale) );
  }
  return LogLikelihood;
}

//************** Updating (Public) **********************************************************
void Individual::SampleLocusAncestry(const Options* const options){
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    Chromosome *C = Loci->getChromosome(j);
    // update of forward probs here is unnecessary if SampleTheta was called and proposal was accepted  
    //Update Forward/Backward probs in HMM
    if( !logLikelihood.HMMisOK ) {
      UpdateHMMInputs(j, options, Theta, _rho);
    }
    // sampling locus ancestry can use current values of forward probability vectors alpha in HMM 
    C->SampleLocusAncestry(LocusAncestry[j], (!isHaploid && (!Loci->isXChromosome(j) || SexIsFemale)));
  } //end chromosome loop
}

/**
 * Return unordered probs as a vector of vectors of doubles.
 * Function AncestryAssocTest::Update() wants them this way.
 */
vector<vector<double> >& Individual::getUnorderedProbs(
  const unsigned int j)
{
  return UnorderedProbs[j];
}

/** 
    function to calculate genotype probs as an average
    over conditional probs of hidden states.
    
    call this function from IndividualCollection.cc just after SampleLocusAncestry
    call new function in Chromosome to get state probs
    unchanged state probs from HMM
 */

void Individual::calculateUnorderedGenotypeProbs(void){
  unsigned int numberCompositeLoci = Loci->GetNumberOfCompositeLoci();
  for(unsigned int j = 0; j < numberCompositeLoci; ++j) {
    calculateUnorderedGenotypeProbs(j);
  }
  return;
}

/**
   Same as Individual::calculateUnorderedProbs(void),
   but for j^th locus only
 */
void Individual::calculateUnorderedGenotypeProbs(unsigned j){
  if (isHaploidIndividual()) {
    string s = "Individual::calculateUnorderedGenotypeProbs(int j) not implemented for haploid individuals";
    throw(s);
  }
  vector<double> orderedGenotypeProbs(4);
  vector<double> orderedStateProbs = vector<double>(Populations * Populations);
  vector<hapPair> hp;
  int anc[2];
  
  if (not (Loci->GetNumberOfStates(j) == 2)) {
    throw string("Trying to calculate UnorderedProbs but Loci->GetNumberOfStates(j) != 2");
    return;
  }
  
  const int numberOfHiddenStates = Populations;
  
  // code below should be executed as a loop over all K^2 states
  // of anc
  // GetLocusAncestry(j, anc);
  // chromosome, locus = Loci->GetChrmAndLocus(j);
  orderedStateProbs = getStateProbs(not this->isHaploidIndividual(), Loci->getChromosomeNumber(j), Loci->getRelativeLocusNumber(j));

  hp = PossibleHapPairs[j];
  // set UnorderedProbs[j][*][0] to 0;
    vector<vector<double> >::iterator gi;
    for (gi = UnorderedProbs[j].begin(); gi != UnorderedProbs[j].end(); ++gi) {
      (*gi)[0] = 0;
    }
 
    int ospIdx;
    orderedGenotypeProbs.assign(4, 0);

    for (anc[0] = 0; anc[0] < numberOfHiddenStates; ++anc[0]) {
      for (anc[1] = 0; anc[1] < numberOfHiddenStates; ++anc[1]) {
        ospIdx = anc[0] * numberOfHiddenStates + anc[1];
        
        /* Possible optimization: if the probability of the state
       * (orderedStateProbs[ospIdx]) is close to zero, it might have
       * a very little effect on the results, so this state could
       * be skipped. Unfortunately, the threshold of 1e-7 is
       * still too high.
       * // if (orderedStateProbs[ospIdx] > 1e-7) {...}
       */
       
      (*Loci)(j)->getConditionalHapPairProbs(orderedGenotypeProbs, hp, anc);

      /* multiply result by conditional probs of anc and accumulate
       * result in array genotype probs (size 3 x number of loci) */
      for (int ogpi = 0; ogpi < 4; ++ogpi) {
        orderedGenotypeProbs[ogpi] *= orderedStateProbs[ospIdx];
      }
      
      // TODO: Rename UnorderedProbs
   	  //P(no copies of allele2)
      UnorderedProbs[j][0][0] += orderedGenotypeProbs[0];
      //P(1 copy of allele2)
   	  UnorderedProbs[j][1][0] += (orderedGenotypeProbs[1] + orderedGenotypeProbs[2]);
      //P(2 copies of allele2)
   	  UnorderedProbs[j][2][0] += orderedGenotypeProbs[3];
    }
  }
}

vector<double> Individual::getStateProbs(const bool isDiploid,const int chromosome,const int locus)const{
  return Loci->getChromosome(chromosome)->getHiddenStateProbs(isDiploid, locus);
}

/**
   Accumulate counts of arrivals of each state. ConcordanceCounts is a L * 2K  array, where the first K elements
   in each row are counts of discordant loci and the remaining K are counts of concordant loci
*/
void Individual::AccumulateConcordanceCounts(int* ConcordanceCounts)const{
  unsigned locus = 0;
  const unsigned  K = Populations;
  // const unsigned KSq = Populations * Populations;
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    const Chromosome* C = Loci->getChromosome(j);
    ++locus;//skip first locus on each chromosome
    for(unsigned locus = 1; locus < C->GetSize(); ++locus){
      
      //first gamete
      if( LocusAncestry[j][locus-1] != LocusAncestry[j][locus]) //discordant loci
        ++ConcordanceCounts[ locus*2*K + LocusAncestry[j][locus] ];
      else//concordant loci
        ++ConcordanceCounts[ locus*2*K + K + LocusAncestry[j][locus] ];
      
      //second gamete
      if(!j==X_posn || SexIsFemale){
        if( LocusAncestry[j][C->GetSize() + locus-1] != LocusAncestry[j][C->GetSize() + locus])
          ++ConcordanceCounts[locus*2*K + LocusAncestry[j][C->GetSize()+locus]];
        else
          ++ConcordanceCounts[locus*2*K + K + LocusAncestry[j][C->GetSize()+locus]];
 	
      }
      ++locus;
    }
  } //end chromosome loop
}

#ifdef PARALLEL
void Individual::SampleHapPair(unsigned j, unsigned jj, unsigned locus, AlleleFreqs *A, 
			       bool skipMissingGenotypes, bool annealthermo, bool UpdateCounts,
			       const double* const AlleleProbs){
  if( !skipMissingGenotypes || !GenotypesMissing[j][jj]){
    int ancestry[2];//to store ancestry states
    GetLocusAncestry(j,jj,ancestry);
    //might be a shortcut for haploid data since there is only one compatible hap pair, no need to sample
    if(PossibleHapPairs[locus].size()> 1){// no need to sample if only one possible hap pair

     //code taken from CompositeLocus. workers have no CompositeLocus objects so this must be done here
      const unsigned NumberOfStates = Loci->GetNumberOfStates(locus);
      double* Probs = new double[PossibleHapPairs[locus].size()];//1way array of hap.pair probs
      double* p = Probs;
      happairiter end = PossibleHapPairs[locus].end();
      for(happairiter hiter = PossibleHapPairs[locus].begin() ; hiter != end ; ++hiter, ++p) {
	*p = AlleleProbs[ancestry[0] * NumberOfStates + hiter->haps[0] ] 
	  * AlleleProbs[ancestry[1] * NumberOfStates + hiter->haps[1] ];
      }
      
      int h = Rand::SampleFromDiscrete(Probs, PossibleHapPairs[locus].size());
      delete[] Probs;
      
      sampledHapPairs[locus].haps[0] = PossibleHapPairs[locus][h].haps[0];
      sampledHapPairs[locus].haps[1] = PossibleHapPairs[locus][h].haps[1];

    }
    if(UpdateCounts && !GenotypesMissing[j][jj])
      //now update allelecounts in AlleleFreqs using sampled hap pair
      A->UpdateAlleleCounts(locus, sampledHapPairs[locus].haps, ancestry, (gametes[j]==2), annealthermo);
  }
}
#else
void Individual::SampleHapPair(unsigned j, unsigned jj, unsigned locus, AlleleFreqs *A, bool skipMissingGenotypes, bool annealthermo, bool UpdateCounts){
  if( !skipMissingGenotypes || !GenotypesMissing[j][jj]) {
    int anc[2];//to store ancestry states
    GetLocusAncestry(j,jj,anc);
    if(PossibleHapPairs[locus].size() > 1 && !annealthermo) { 
      // no need to sample if only one possible hap pair or if annealing for thermo integration
      (*Loci)(locus)->SampleHapPair(&(sampledHapPairs[locus]), PossibleHapPairs[locus], anc);
    }
    // now update allelecounts in AlleleFreqs using sampled hap pair
    // UpdateAlleleCounts does nothing if annealthermo and > 2 alleles 
    if(UpdateCounts && !GenotypesMissing[j][jj])
      A->UpdateAlleleCounts(locus, sampledHapPairs[locus].haps, anc, (gametes[j]==2), annealthermo);
  }
}
#endif

void Individual::UpdateAlleleCounts(unsigned j, unsigned jj, unsigned locus, AlleleFreqs *A, bool annealthermo)const{
  if(!GenotypesMissing[j][jj]){
    int anc[2];//to store ancestry states
    GetLocusAncestry(j,jj,anc);
    A->UpdateAlleleCounts(locus, sampledHapPairs[locus].haps, anc, (gametes[j]==2), annealthermo);
  }
}

void Individual::SampleMissingOutcomes(DataMatrix *Outcome, const vector<Regression*>& R){
  int NumOutcomes = Outcome->nCols();
  // sample missing values of outcome variable
  for( int k = 0; k < NumOutcomes; k++ ){
    if( Outcome->isMissing( myNumber-1, k ) ){
      if( R[k]->getRegressionType() == Linear)
	Outcome->set( myNumber-1, k, Rand::gennor( R[k]->getExpectedOutcome(myNumber-1), 1 / sqrt( R[k]->getlambda() ) ));
      else{
	if( Rand::myrand() * R[k]->getExpectedOutcome(myNumber-1) < 1 )
	  Outcome->set( myNumber-1, k, 1);
	else
	  Outcome->set( myNumber-1, k, 0);
      }
    }
  }
}

// const int Individual::getNumberOfHiddenStates()
// {
//   return Populations;
// }
