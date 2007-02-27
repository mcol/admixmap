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

#define PR(x) cout << #x << " = " << x << endl;

// ******* Static Member Declarations
bool Individual::Xdata;
unsigned Individual::X_posn;
unsigned int Individual::numChromosomes;
Genome *Individual::Loci;
int Individual::Populations;

///for printing a happair
std::ostream& operator<<(std::ostream& os, const hapPair &h){
  os << h.haps[0] << " " << h.haps[1];
  return os;
}

//******** Constructors **********
Individual::Individual() {//should initialise pointers here
}
Individual::Individual(int number, const Options* const options, const InputData* const Data) {
  Initialise(number, options, Data);
}

void Individual::Initialise(int number, const Options* const options, const InputData* const Data) {
  myNumber = number;
  NumGametes = 1;
  GenotypesMissing = new bool*[numChromosomes];
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    GenotypesMissing[j] = new bool[ Loci->GetSizeOfChromosome(j) ];
  }  
  missingGenotypes = 0;//allocated later, if needed
  //retrieve genotypes
  Data->GetGenotype(myNumber, options->getgenotypesSexColumn(), *Loci, &genotypes, GenotypesMissing);
  isHaploid = (bool)(genotypes[0][0].size()==1);//note: assumes at least one autosome before X-chr
  if( options->isRandomMatingModel() && !isHaploid) NumGametes = 2;
  
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

  Theta = new double[ Populations * NumGametes ];
  
  SetUniformAdmixtureProps();
  
  // vector of possible haplotype pairs - 2 integers per locus if diploid, 1 if haploid 
  PossibleHapPairs = new vector<hapPair>[numCompositeLoci];

  LocusAncestry = new int*[ numChromosomes ]; // array of matrices in which each col stores 2 integers   
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
}

//********** Destructor **********
Individual::~Individual() {
  delete[] PossibleHapPairs;

  if(GenotypesMissing){
    for(unsigned j = 0; j < numChromosomes; ++j)delete[] GenotypesMissing[j];
    delete[] GenotypesMissing;
  }
  if(LocusAncestry){
    for(unsigned j = 0; j < numChromosomes; ++j)delete[] LocusAncestry[j];
    delete[] LocusAncestry;
  }

  delete[] Theta;
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
///sets possible hap pairs for a single SNP
void Individual::SetPossibleHaplotypePairs(const vector<vector<unsigned short> > Genotype, vector<hapPair> &PossibleHapPairs){
  if(Genotype.size()!=1)throw string("Invalid call to Individual::SetPossibleHapPairs()");
  hapPair hpair;
  PossibleHapPairs.clear();
  if(Genotype[0][0] == 0 || Genotype[0][1]==0){//missing genotype
    hpair.haps[0] = 0; hpair.haps[1] = 0;
    PossibleHapPairs.push_back(hpair);//(1,1)
    hpair.haps[1] = 1;
    PossibleHapPairs.push_back(hpair);//(1,2)
    hpair.haps[0] = 1; hpair.haps[1] = 0;
    PossibleHapPairs.push_back(hpair);//(2,1)
    hpair.haps[1] = 1;
    PossibleHapPairs.push_back(hpair);//(2,2)
  }
  //case of homozygote - only one possible happair
  else if(Genotype[0][0] == Genotype[0][1]){
    hpair.haps[0] = hpair.haps[1] = Genotype[0][0]-1;
    PossibleHapPairs.push_back(hpair);
  }
  else{//heterozygote - two possibilities
    hpair.haps[0] = 0; hpair.haps[1] = 1;
    PossibleHapPairs.push_back(hpair);//(1,2)
    hpair.haps[0] = 1; hpair.haps[1] = 0;
    PossibleHapPairs.push_back(hpair);//(2,1)
  }
}

void Individual::setGenotypesToMissing(){
  for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c)
    for(unsigned j = 0; j < Loci->GetSizeOfChromosome(c); ++j)
      GenotypesMissing[c][j] = true;
}

void Individual::DeleteGenotypes(){
  for(unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); ++j){
#ifdef PARALLEL
      genotypes[j][0].clear();
#else
    for(int k = 0; k < Loci->getNumberOfLoci(j); ++k)
      genotypes[j][k].clear();
#endif
    genotypes[j].clear();
  }
  genotypes.clear();
}

void Individual::SetMissingGenotypes(){
  //allocates and sets an array of bools indicating whether genotypes at each locus are missing
  //used in HW score test; NB call before genotypes are deleted
  if(genotypes.size()==0)throw string("determining missing genotypes after genotypes have been deleted");
  missingGenotypes = new bool[Loci->GetTotalNumberOfLoci()];
  unsigned index = 0;
  for(unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); ++j)
    for(int k = 0; k < Loci->getNumberOfLoci(j); ++k){
      missingGenotypes[index++] = (genotypes[j][k][0] == 0);
    }
}

//********** sets static members, including allocation and deletion of static objects for score tests
void Individual::SetStaticMembers(Genome* const pLoci, const Options* const options){
  Loci = pLoci;
  numChromosomes = Loci->GetNumberOfChromosomes();
  Populations = options->getPopulations();
  Xdata = Loci->isX_data();
  X_posn = 9999; //position of the X chromosome in the sequence of chromosomes in the input data
  if(Xdata) X_posn = Loci->GetChrNumOfLocus(Loci->getFirstXLocus());//too clunky, should simplify
}

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

void Individual::GetLocusAncestry(int locus, int Ancestry[2])const{
  unsigned c, l;
  Loci->GetChrmAndLocus(locus, &c, &l);
  GetLocusAncestry(c, l, Ancestry);
}
void Individual::GetLocusAncestry(int chrm, int locus, int Ancestry[2])const {
  Ancestry[0]  = LocusAncestry[chrm][locus];
  if(isHaploid || ((unsigned)chrm == X_posn && !SexIsFemale))Ancestry[1] = Ancestry[0];
  else Ancestry[1] = LocusAncestry[chrm][Loci->GetSizesOfChromosomes()[chrm]  + locus];
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
    Chromosome* C = Loci->getChromosome(j);
    // update of forward probs here is unnecessary if SampleTheta was called and proposal was accepted  
      //Update Forward/Backward probs in HMM
      if( !logLikelihood.HMMisOK ) {
	UpdateHMMInputs(j, options, Theta, _rho);
      }
      // sampling locus ancestry can use current values of forward probability vectors alpha in HMM 
      C->SampleLocusAncestry(LocusAncestry[j], (!isHaploid && (!Loci->isXChromosome(j) || SexIsFemale)));
  } //end chromosome loop
}

//ADDHERE - new method to calculate genotype probs as an average over conditional probs of hidden states
// call this method from IndividualCollection.cc just after SampleLocusAncestry
//
// call new method in Chromosome to get state probs
// unchanged state probs from HMM


// code below should be executed as a loop over all K^2 states of anc


// 	    (*Lociptr)(j)->getConditionalHapPairProbs(OrderedProbs, ind->getPossibleHapPairs(j), anc);
// #else
//             //TODO: write alternative for parallel version
// #endif
// 	    UnorderedProbs[0] = vector<double>(1, OrderedProbs[0]);//P(no copies of allele2)
// 	    UnorderedProbs[1] = vector<double>(1, OrderedProbs[1] + OrderedProbs[2]);//P(1 copy of allele2)
// 	    UnorderedProbs[2] = vector<double>(1, OrderedProbs[3]);//P(2 copies of allele2)
// now multiply result by conditional probs of anc and accumulate result in array genotype probs (size 3 x number of loci)
// also add a public method to get the genotype probs for this individual

void Individual::AccumulateAncestry(int* SumAncestry){
  unsigned locus = 0;
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    Chromosome* C = Loci->getChromosome(j);
    ++locus;//skip first locus on each chromosome
      for(unsigned l = 1; l < C->GetSize(); ++l){
	if( LocusAncestry[j][l-1] != LocusAncestry[j][l])//first gamete
	  ++SumAncestry[locus*2];
	else
	  ++SumAncestry[locus*2 + 1];
	if(!isHaploid && ((j!=X_posn) || SexIsFemale)){//second gamete
	if( LocusAncestry[j][C->GetSize() + l-1] != LocusAncestry[j][C->GetSize() + l])
	  ++SumAncestry[locus*2];
	else
	  ++SumAncestry[locus*2 + 1];
	}
	++locus;
      }
  } //end chromosome loop
}
#ifdef PARALLEL
void Individual::SampleHapPair(unsigned j, unsigned jj, unsigned locus, AlleleFreqs *A, bool skipMissingGenotypes, bool annealthermo, bool UpdateCounts,
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
	*p = AlleleProbs[ancestry[0] * NumberOfStates + hiter->haps[0] ] * AlleleProbs[ancestry[1] * NumberOfStates + hiter->haps[1] ];
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

