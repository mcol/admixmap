/** 
 *   ADMIXMAP
 *   Individual.cc 
 *   Class to represent an individual and update individual-level parameters
 *   Copyright (c) 2002-2006 LSHTM
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "Individual.h"
#include "StringConvertor.h"
#include <algorithm>
#include <limits>
#include <sstream>

#define PR(x) cout << #x << " = " << x << endl;

// ******* Static Member Declarations
double *Individual::LikRatio1;
double *Individual::LikRatio2;
double *Individual::AffectedsScore = 0;
double *Individual::AffectedsVarScore = 0;
double *Individual::AffectedsInfo = 0;
double **Individual::AncestryScore = 0;
double **Individual::AncestryInfo = 0;
double **Individual::AncestryVarScore = 0;
double **Individual::AncestryInfoCorrection = 0;
double *Individual::B;
double *Individual::PrevB;
double *Individual::Xcov;
bool Individual::Xdata;
unsigned Individual::X_posn;
unsigned Individual::NumIndGametes;

//#define TruncationPt 99 // upper truncation point for sum intensities parameter rho

unsigned int Individual::numChromosomes;
Genome *Individual::Loci;
int Individual::Populations;

void operator<<(std::ostream& os, const hapPair &h){
  os << h.haps[0] << " " << h.haps[1];
}

//******** Constructors **********
Individual::Individual() {//should initialise pointers here
}

Individual::Individual(int number, const AdmixOptions* const options, const InputData* const Data,  
		       bool undertest=false) {
  myNumber = number;
  IAmUnderTest = undertest;
  
  double init=0.0;
  if( !options->isGlobalRho() && options->getIndAdmixHierIndicator()) {//model with individual- or gamete-specific sumintensities
    //set prior mean as initial value for rho 
    if(options->getRhobetaShape() > 1) { 
      init = options->getRhoalpha() * options->getRhobetaRate() / (options->getRhobetaShape() - 1 );
    } else {
      init = options->getRhoalpha() * options->getRhobetaRate() / options->getRhobetaShape() ;//conditional prior mean
    } 
  } else { // no hierarchical model or globalrho prior
    init = options->getRhoalpha() / options->getRhobeta();
  }
  _rho.assign(NumIndGametes, 0.0); // set to 0 for unadmixed gametes 
  rhohat.assign(NumIndGametes, 0.0); // set to 0 for unadmixed gametes 
  
  for(unsigned g = 0; g < NumIndGametes; ++g) {
    if(options->isAdmixed(g)) {
      _rho[g] = init;
      rhohat[g] = init;
    }
  }
  // for global rho, Latent should assign initial value 
  sumlogrho.assign(_rho.size(), 0.0);
  
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

  Theta = 0;
  thetahat = 0;
  SumSoftmaxTheta = 0;
  ThetaProposal = 0;
  dirparams = 0;
  // SumLocusAncestry is sum of locus ancestry states over loci at which jump indicator xi is 1 
  SumLocusAncestry = 0; 
  SumLocusAncestry_X = 0;
  
  Theta = new double[ Populations * NumIndGametes ];
  thetahat = new double[ Populations * NumIndGametes ];
  SumSoftmaxTheta = new double[ Populations * NumIndGametes ];
  fill(SumSoftmaxTheta, SumSoftmaxTheta + Populations*NumIndGametes, 0.0);
  
  if(!options->getHapMixModelIndicator()) {
    dirparams = new double[Populations]; //to hold dirichlet parameters for conjugate updates of theta
    ThetaProposal = new double[ Populations * NumIndGametes ];
    SumLocusAncestry = new int[Populations * 2];
    SumLocusAncestry_X = new int[Populations * 2];
    if(!options->isGlobalRho())SumNumArrivals.resize(2*numCompositeLoci);  
  } else { 
    SetDefaultAdmixtureProps();
  }
  //LogPrior = 0.0;
  
  // vector of possible haplotype pairs - 2 integers per locus if diploid, 1 if haploid 
  PossibleHapPairs = new vector<hapPair>[numCompositeLoci];
  
  GenotypesMissing = new bool*[numChromosomes];  
  missingGenotypes = 0;//allocated later, if needed
  LocusAncestry = new int*[ numChromosomes ]; // array of matrices in which each col stores 2 integers   
  
  size_t AncestrySize = 0;  // set size of locus ancestry array
  //gametes holds the number of gametes for each chromosome, either 1 or 2
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    if( !Loci->isXChromosome(j) ){ // if not X chromosome, set number of elements to 2 * num loci
      AncestrySize = 2 * Loci->GetSizeOfChromosome(j) ;
      gametes.push_back(2);
    } else if( !SexIsFemale ) {//male or missing
      AncestrySize = Loci->GetSizeOfChromosome(j) ;
      gametes.push_back(1);
      //X_posn = j;
    } else{//female
      AncestrySize = 2 * Loci->GetSizeOfChromosome(j) ;
      gametes.push_back(2);
      //X_posn = j;
    } 
    GenotypesMissing[j] = new bool[ Loci->GetSizeOfChromosome(j) ];
    LocusAncestry[j] = new int[ AncestrySize];
    for(unsigned i = 0; i < AncestrySize; ++i) LocusAncestry[j][i] = 0;
  }
  //retrieve genotypes
  Data->GetGenotype(myNumber, options->getgenotypesSexColumn(), *Loci, &genotypes, GenotypesMissing);
  // loop over composite loci to set possible haplotype pairs compatible with genotype 
  for(unsigned j = 0; j < (unsigned)numCompositeLoci; ++j) {
#ifdef PARALLEL
    SetPossibleHaplotypePairs(genotypes[j], PossibleHapPairs[j]); 
    //NOTE: X data not yet supported in parallel version
#else
    //     if( (Loci->GetChrNumOfLocus(j)==X_posn) && (sex==male)){
    //       (*Loci)(j)->setPossibleXHaplotypes(genotypes[j], PossibleHapPairs[j]);
    //     }
    //     else
    (*Loci)(j)->setPossibleHaplotypePairs(genotypes[j], PossibleHapPairs[j]);
#endif
    
    //initialise sampledHapPairs with the first of the possible happairs. Then, if there is only one, sampling of hap pair can be skipped.
    sampledHapPairs.push_back(PossibleHapPairs[j][0]);
  }
  //Now the PossibleHapPairs have ben determined and missing genotype indicators have been set, 
  //the genotypes are deleted as they are no longer needed 
  if( options->getHWTestIndicator())SetMissingGenotypes();
  DeleteGenotypes();
  
  //initialise genotype probs array and array of indicators for genotypes missing at locus
  GenotypeProbs = new double*[numChromosomes];
  // unsigned locus = 0;
  for(unsigned j = 0; j < numChromosomes; ++j) {
    if( (j==X_posn) && !SexIsFemale)
      GenotypeProbs[j] = new double[Loci->GetSizeOfChromosome(j)*Populations];
    else
      GenotypeProbs[j] = new double[Loci->GetSizeOfChromosome(j)*Populations*Populations];
    //     for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ){
    //       SetGenotypeProbs(j, jj, locus, false);
    //       locus++;
    //     }
  }
  
  // ** set up StepSizeTuner object for random walk updates of admixture **
  NumberOfUpdates = 0;
  w = 1;
  step0 = 0.3; // initial sd of random walk proposal distribution 
  step = step0;
  ThetaTuner.SetParameters( step0, 0.00, 10.0, 0.44);  
  
  logLikelihood.value = 0.0;
  logLikelihood.ready = false;
  logLikelihood.HMMisOK = false;
}

//********** Destructor **********
Individual::~Individual() {
  delete[] PossibleHapPairs;
  if(GenotypeProbs){
    for(unsigned j = 0; j < numChromosomes; ++j) delete[] GenotypeProbs[j];
    delete[] GenotypeProbs; 
  }
  if(GenotypesMissing){
    for(unsigned j = 0; j < numChromosomes; ++j)delete[] GenotypesMissing[j];
    delete[] GenotypesMissing;
  }
  if(LocusAncestry){
    for(unsigned j = 0; j < numChromosomes; ++j)delete[] LocusAncestry[j];
    delete[] LocusAncestry;
  }

  delete[] Theta;
  delete[] SumSoftmaxTheta;
  delete[] missingGenotypes;
  delete[] dirparams;
  delete[] ThetaProposal;
  delete[] SumLocusAncestry;
  delete[] SumLocusAncestry_X;
}

/// draw initial values for admixture proportions theta from Dirichlet prior 
void Individual::drawInitialAdmixtureProps(const std::vector<std::vector<double> > &alpha) {
  size_t K = Populations;
  for( unsigned g = 0; g < NumIndGametes; ++g ) { 
    double sum = 0.0;
    for(size_t k = 0; k < K; ++k) { 
      sum += alpha[g][k];
    }
    for(size_t k = 0; k < K; ++k) { 
      thetahat[g*K+k] = alpha[g][k] / sum; // set thetahat to prior mean
      dirparams[k] = alpha[g][k]; 
    }
    // draw theta from Dirichlet with parameters dirparams
    Rand::gendirichlet(K, dirparams, Theta+g*K ); 
  }  
}

void Individual::SetDefaultAdmixtureProps() {
  size_t K = Populations;
  for( unsigned g = 0; g < NumIndGametes; ++g ) { 
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

#ifdef PARALLEL
//this version can also be used in non-parallel version
void Individual::SetGenotypeProbs(int j, int jj, unsigned locus, const double* const AlleleProbs){
  if( !GenotypesMissing[j][jj] ){
    double *q = GenotypeProbs[j]+jj*Populations*Populations;
    
    happairiter end = PossibleHapPairs[locus].end();
    const unsigned NumberOfStates = Loci->GetNumberOfStates(locus);
    for(int k0 = 0; k0 < Populations; ++k0)
      for(int k1 = 0; k1 < Populations; ++k1) {
	*q = 0.0;
	happairiter h = PossibleHapPairs[locus].begin();
	for( ; h != end ; ++h) {
	  *q += AlleleProbs[k0*NumberOfStates+h->haps[0]] * AlleleProbs[k1*NumberOfStates+h->haps[1]];
	}
	q++;
      }
    
  } else {
    for( int k = 0; k < Populations*Populations; ++k ) GenotypeProbs[j][jj*Populations*Populations + k] = 1.0;
  }
}
#endif

void Individual::SetGenotypeProbs(int j, int jj, unsigned locus, bool chibindicator=false){
  //chibindicator is passed to CompositeLocus object.  If set to true, CompositeLocus will use HapPairProbsMAP
  //instead of HapPairProbs when allelefreqs are not fixed.
  if( !GenotypesMissing[j][jj] ) {
    if( j!=(int)X_posn || SexIsFemale) { //diploid genotype
      (*Loci)(locus)->GetGenotypeProbs(GenotypeProbs[j]+jj*Populations*Populations, PossibleHapPairs[locus], 
				       chibindicator);
    } else {//haploid genotype
      (*Loci)(locus)->GetHaploidGenotypeProbs(GenotypeProbs[j]+jj*Populations, PossibleHapPairs[locus], 
					      chibindicator);
    }
  } else {
    if( j!=(int)X_posn || SexIsFemale)  //diploid genotype
      for( int k = 0; k < Populations*Populations; ++k ) GenotypeProbs[j][jj*Populations*Populations + k] = 1.0;
    else //haploid genotype
      for( int k = 0; k < Populations; ++k ) GenotypeProbs[j][jj*Populations + k] = 1.0;
  }
}

void Individual::AnnealGenotypeProbs(int j, const double coolness) {
  // called after energy has been evaluated, before updating model parameters
  int locus = Loci->getChromosome(j)->GetLocus(0);
  for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ){ // loop over composite loci
    if( !GenotypesMissing[j][jj] ) { 
      if( j!=(int)X_posn || SexIsFemale)  //diploid genotype
	for(int k = 0; k < Populations*Populations; ++k) // loop over ancestry states
	  GenotypeProbs[j][jj*Populations*Populations+k] = pow(GenotypeProbs[j][jj*Populations*Populations+k], coolness); 
      else //haploid genotype
	for(int k = 0; k < Populations; ++k) // loop over ancestry states
	  GenotypeProbs[j][jj*Populations+k] = pow(GenotypeProbs[j][jj*Populations+k], coolness); 
    }
    locus++;
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
void Individual::SetStaticMembers(Genome* const pLoci, const AdmixOptions* const options){
  Loci = pLoci;
  numChromosomes = Loci->GetNumberOfChromosomes();
  Populations = options->getPopulations();
  Xdata = Loci->isX_data();
  X_posn = 9999; //position of the X chromosome in the sequence of chromosomes in the input data
  if(Xdata) X_posn = Loci->GetChrNumOfLocus(Loci->getFirstXLocus());//too clunky, should simplify
  NumIndGametes = 1;
  if( options->isRandomMatingModel() ) NumIndGametes = 2;
  int K = Populations;
  int L = Loci->GetNumberOfCompositeLoci();

  LikRatio1 = 0;
  LikRatio2 = 0;
  AncestryScore = 0;
  AncestryInfo = 0;
  AncestryVarScore = 0;
  AncestryInfoCorrection = 0;
  AffectedsScore = 0;
  AffectedsVarScore = 0;
  AffectedsInfo = 0;
  B = 0;
  PrevB = 0;
  Xcov = 0;

  if( options->getTestForAffectedsOnly() ){
    if(Populations == 2) K = 1;
    AffectedsScore = new double[L * K];
    AffectedsVarScore = new double[L * K];
    AffectedsInfo = new double[L * K];
    LikRatio1 = new double[L*K];
    LikRatio2 = new double[L*K];
    fill(LikRatio1, LikRatio1+L*K, 0.0);
    fill(LikRatio2, LikRatio2+L*K, 0.0);
  }
  if( options->getTestForLinkageWithAncestry() ){
    AncestryScore = alloc2D_d(L, 2*Populations);
    AncestryInfo = alloc2D_d(L, 4*Populations*Populations);
    AncestryVarScore = alloc2D_d(L, Populations);
    AncestryInfoCorrection = alloc2D_d(L, Populations);
    B = new double[Populations * Populations];
    PrevB = new double[Populations * Populations];
    Xcov = new double[Populations];
  }
}

void Individual::DeleteStaticMembers(){
  delete[] AffectedsScore;
  delete[] AffectedsInfo;
  delete[] AffectedsVarScore;
  delete[] LikRatio1;
  delete[] LikRatio2;
  delete[] B;
  delete[] PrevB;
  delete[] Xcov;
  free_matrix(AncestryScore, Loci->GetNumberOfCompositeLoci());
  free_matrix(AncestryInfo, Loci->GetNumberOfCompositeLoci());  
  free_matrix(AncestryVarScore, Loci->GetNumberOfCompositeLoci());
  free_matrix(AncestryInfoCorrection, Loci->GetNumberOfCompositeLoci());
}

//********** set admixture proportions *********

void Individual::setAdmixtureProps(const double* const a, size_t size) {
  //TODO: size arg not necessary should equal NumIndGametes*K
  for(unsigned i = 0; i < size; ++i)  {
    Theta[i] = a[i];
  }
}

void Individual::HMMIsBad(bool loglikisbad) {
  logLikelihood.HMMisOK = false;
  if(loglikisbad)logLikelihood.ready = false;
}

//******************** Accessors ***********************************************************

const std::vector<hapPair > &Individual::getPossibleHapPairs(unsigned int locus)const{
  return PossibleHapPairs[locus];
}

const int* Individual::getSampledHapPair(int locus)const{
  return sampledHapPairs[locus].haps;
}

const double* Individual::getAdmixtureProps()const {
  return Theta;
}

///returns sum of sumintensities over gametes
double Individual::getSumrho()const { 
  double sumrho = 0;
  for( unsigned int g = 0; g < _rho.size(); g++ )
    sumrho += _rho[g];
  return sumrho;
}

const vector<double> Individual::getRho()const {
   return _rho;
}

void Individual::GetLocusAncestry(int locus, int Ancestry[2])const{
  unsigned c, l;
  Loci->GetChrmAndLocus(locus, &c, &l);
  GetLocusAncestry(c, l, Ancestry);
}
void Individual::GetLocusAncestry(int chrm, int locus, int Ancestry[2])const {
  Ancestry[0]  = LocusAncestry[chrm][locus];
  if((unsigned)chrm == X_posn && !SexIsFemale)Ancestry[1] = Ancestry[0];
  else Ancestry[1] = LocusAncestry[chrm][Loci->GetSizesOfChromosomes()[chrm]  + locus];
}

///returns value of LocusAncestry at a locus for a particular gamete
int Individual::GetLocusAncestry(int chrm, int gamete, int locus)const{
  int g = (gametes[chrm] == 2) ? gamete : 0; //so that gamete = 1 works when gametes[chrm] = 1;
  return LocusAncestry[chrm][g * Loci->GetSizesOfChromosomes()[chrm]  + locus] ;
}

const int *Individual::getSumLocusAncestry()const{
  return SumLocusAncestry;
}
const int *Individual::getSumLocusAncestryX()const{
  return SumLocusAncestry_X;
}
  //returns number of arrivals across genome
const vector<unsigned>Individual::getSumNumArrivals()const{
  vector<unsigned> SumN(2,0);
  unsigned locus = 0;
  for(unsigned i = 0; i < numChromosomes; ++i){
      for(unsigned j = 0; j < Loci->GetSizeOfChromosome(i); ++j, ++locus)    
	if(i != X_posn){//skip X chromosome
	  //accumulate sums over loci
	  SumN[0] += SumNumArrivals[locus*2];//first gamete
	  SumN[1] += SumNumArrivals[locus*2+1];//second gamete
	}
  }
  return SumN;
}

const vector<unsigned>Individual::getSumNumArrivals_X()const{
  vector<unsigned> SumN(2,0);
  if(Xdata){//accumulate sums over Loci on X chromosome
    unsigned offset = Loci->getFirstXLocus();
    for(unsigned i = 0; i < Loci->GetSizeOfChromosome(X_posn); ++i){
      SumN[0] += SumNumArrivals[2*(i + offset)];
      SumN[1] += SumNumArrivals[2*(i + offset)+1];
    }
  }
  return SumN;
}
void Individual::getSumNumArrivals(std::vector<unsigned> *sum)const{
  //should check size of argument
  for(unsigned i = 0; i < sum->size(); ++i){//sum over gametes
    (*sum)[i] += SumNumArrivals[2*i] + SumNumArrivals[2*i + 1];
  }
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
//****************** Log-Likelihoods **********************
// public function: 
// calls private function to get log-likelihood at current parameter values, and stores it either as loglikelihood.value or as loglikelihood.tempvalue
// store should be false when calculating energy for an annealed run, or when evaluating proposal for global sum-intensities
double Individual::getLogLikelihood( const AdmixOptions* const options, const bool forceUpdate, const bool store) {

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
double Individual::getLogLikelihood(const AdmixOptions* const options, const double* const theta, 
				    const vector<double > rho,  bool updateHMM) {
  double LogLikelihood = 0.0;
  if(Populations == 1) LogLikelihood = getLogLikelihoodOnePop();
  else { 
    for( unsigned int j = 0; j < numChromosomes; j++ ){
	if(updateHMM){// force update of forward probs 
	    UpdateHMMInputs(j, options, theta, rho);
	}
	LogLikelihood += Loci->getChromosome(j)->getLogLikelihood( !Loci->isXChromosome(j) || SexIsFemale );
    }
  }
  
//   cout << "From getLogLikelihood: Populations " << Populations << endl;
//   for(int g = 0; g < NumIndGametes; ++g) {
//       cout << "rho " << rho[g] << " theta";
//       for(int k=0; k < Populations; ++k) {
// 	cout << theta[g*Populations+k] << " ";
//       }
//   }
//   cout << endl;
//   cout << "LogL " << LogLikelihood << endl;
  
  return LogLikelihood; // argument updateHMM is unnecessary - why call this function unless you want an HMM update
  // if HMM update not required, can just use stored log-likelihood  
}

void Individual::storeLogLikelihood(const bool setHMMAsOK) { // to call if a Metropolis proposal is accepted
    logLikelihood.value = logLikelihood.tempvalue; 
    logLikelihood.ready = true;
    if(setHMMAsOK) logLikelihood.HMMisOK = true; 
}                               

double Individual::getLogLikelihoodAtPosteriorMeans(const AdmixOptions* const options) {
  // should set allele freqs also to posterior means, and recalculate prob genotypes at these freqs before calling getloglikelihood 
  //obtain ergodic averages of (inv_softmax)admixture props and (log)sumintensities and transform back
  //to original scales
  unsigned size = Populations * NumIndGametes;
  for(unsigned i = 0; i < size; ++i)SumSoftmaxTheta[i] /= (options->getTotalSamples() - options->getBurnIn());
  for(unsigned i = 0; i < _rho.size(); ++i)sumlogrho[i] = exp(sumlogrho[i]/(options->getTotalSamples() - options->getBurnIn()));

  //apply softmax transformation to obtain thetabar
  double* ThetaBar = new double[NumIndGametes*Populations];
  bool* b = new bool[Populations];
  for( unsigned int g = 0; g < NumIndGametes; g++ ){
    for(int k = 0; k < Populations; ++k)if(Theta[g*Populations + k] > 0.0){
      b[k] = true; //to skip elements set to zero
    } else b[k] = false;
    softmax(Populations, ThetaBar+g*Populations, SumSoftmaxTheta+g*Populations, b);
  }

  double LogLikelihood = 0.0;
  if(Populations == 1) LogLikelihood = getLogLikelihoodOnePop();
  else {
    for( unsigned int j = 0; j < numChromosomes; j++ ) {
      UpdateHMMInputs(j, options, ThetaBar, sumlogrho); // sumlogrho is posterior mean of rho
      // should replace 2nd sumlogrho by half its value
      LogLikelihood += Loci->getChromosome(j)->getLogLikelihood( !Loci->isXChromosome(j) || SexIsFemale );
    }
  }
  delete[] b;
  delete[] ThetaBar;
  return LogLikelihood;
}
double Individual::getLogLikelihoodOnePop(){ //convenient for a single population as no arguments required
  double LogLikelihood = 0.0;
  double *Prob;
  Prob = new double[1];//one pop so 1x1 array
  for( unsigned int j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ) {// loop over composite loci
    if(!GenotypeIsMissing(j)){
      (*Loci)(j)->GetGenotypeProbs(Prob,getPossibleHapPairs(j), false);
      LogLikelihood += log( Prob[0] );
    }
  }
  delete[] Prob;
  return LogLikelihood;
}

//************** Updating (Public) ***************************************************************************************
void Individual::ResetSufficientStats(){
  if(Populations>1) {
    // ** reset SumLocusAncestry to zero
    for(int j = 0; j < Populations *2; ++j)SumLocusAncestry[j] = 0;
    // if(Loci->isX_data() ){
    int J = Populations;
    if(SexIsFemale) J *=2;
    for(int j = 0; j < J ;++j) SumLocusAncestry_X[j] = 0;
    //  }

    //SumNumArrivals is the number of arrivals between each pair of adjacent loci
    fill(SumNumArrivals.begin(), SumNumArrivals.end(), 0);
  }
}

void Individual::SampleLocusAncestry(const AdmixOptions* const options){
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    Chromosome* C = Loci->getChromosome(j);
    // update of forward probs here is unnecessary if SampleTheta was called and proposal was accepted  
      //Update Forward/Backward probs in HMM
      if( !logLikelihood.HMMisOK ) {
	UpdateHMMInputs(j, options, Theta, _rho /* , _rho_X*/);
      }
      // sampling locus ancestry can use current values of forward probability vectors alpha in HMM 
      C->SampleLocusAncestry(LocusAncestry[j], (!Loci->isXChromosome(j) || SexIsFemale));
  } //end chromosome loop
}

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
	if((j!=X_posn) || SexIsFemale){//second gamete
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
void Individual::SampleHapPair(unsigned j, unsigned jj, unsigned locus, AlleleFreqs *A, bool hapmixmodel, bool anneal, 
			       const double* const AlleleProbs){
  if( hapmixmodel || !GenotypesMissing[j][jj]){
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

    }//now update allelecounts in AlleleFreqs using sampled hap pair
    A->UpdateAlleleCounts(locus, sampledHapPairs[locus].haps, ancestry, (gametes[j]==2), anneal);
  }
}
#else
void Individual::SampleHapPair(unsigned j, unsigned jj, unsigned locus, AlleleFreqs *A, bool hapmixmodel, bool anneal){
  if( hapmixmodel || !GenotypesMissing[j][jj]){
    int anc[2];//to store ancestry states
    GetLocusAncestry(j,jj,anc);
    //might be a shortcut for haploid data since there is only one compatible hap pair, no need to sample
    if(PossibleHapPairs[locus].size() > 1){// no need to sample if only one possible hap pair
      (*Loci)(locus)->SampleHapPair(&(sampledHapPairs[locus]), PossibleHapPairs[locus], anc);
    }//now update allelecounts in AlleleFreqs using sampled hap pair
    A->UpdateAlleleCounts(locus, sampledHapPairs[locus].haps, anc, (gametes[j]==2), anneal);
  }
}
#endif

void Individual::SampleJumpIndicators(bool sampleArrivals){
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    Chromosome* C = Loci->getChromosome(j);
    // don't need to sample jump indicators if globalrho and no conjugate update of admixture this iteration
    //sample number of arrivals, update SumNumArrivals and SumLocusAncestry
    if( !Loci->isXChromosome(j) )
      C->SampleJumpIndicators(LocusAncestry[j], gametes[j], SumLocusAncestry, SumNumArrivals, 
			      sampleArrivals);
    else 
      C->SampleJumpIndicators(LocusAncestry[j], gametes[j], SumLocusAncestry_X, SumNumArrivals, 
			      sampleArrivals);
  } //end chromosome loop
}

///uses an EM algorithm to search for posterior modes of individual parameters theta and rho
// uses current values of allele freqs 
void Individual::FindPosteriorModes(const AdmixOptions* const options, const vector<vector<double> > &alpha,  
				    double rhoalpha, double rhobeta, ofstream &modefile) {
  unsigned numEMiters = 10;
  unsigned NumEstepiters = 10; 
  double LogUnnormalizedPosterior = - numeric_limits<double>::max( );
  bool isadmixed = options->isAdmixed(0);
  if(NumIndGametes ==2) isadmixed = isadmixed | options->isAdmixed(1);//indicates if individual is admixed
  
  //use current parameter values as initial values
  double *SumLocusAncestryHat = new double[2*Populations];
  for(unsigned EMiter = 0; EMiter < numEMiters; ++EMiter) {
    double SumNumArrivalsHat[2] = {0,0}; //, SumNumArrivals_XHat[2] = {0,0};
    fill(SumLocusAncestryHat, SumLocusAncestryHat + 2*Populations, 0.0);
    
    //E-step: fix theta and rho, sample Locus Ancestry and Number of Arrivals
    NumEstepiters *= 2;
    for(unsigned Estepiters = 0; Estepiters < (unsigned)NumEstepiters ; ++Estepiters) {
      if(Populations >1){
	ResetSufficientStats();
	SampleLocusAncestry(options);
	SampleJumpIndicators((!options->isGlobalRho() || options->getHapMixModelIndicator()));
      }
      vector<unsigned> SumN = getSumNumArrivals();
      vector<unsigned> SumN_X = getSumNumArrivals_X(); // element 0 should be zero in a male
      // accumulate sums
      SumNumArrivalsHat[0] += SumN[0] + SumN_X[0];
      SumNumArrivalsHat[1] += SumN[1] + SumN_X[1];
      transform(SumLocusAncestry, SumLocusAncestry+2*Populations, 
		SumLocusAncestryHat, SumLocusAncestryHat, std::plus<double>());
    }
    // set SumLocusAncestry and SumNumArrivals to their averages over current E step
    for(int i = 0; i < 2*Populations; ++i) {
      SumLocusAncestryHat[i] /= (double)NumEstepiters; 
    }
    cout << "E step\t";
    //     for(unsigned int g = 0; g < NumIndGametes; ++g) {
    //       for(int k = 0; k < Populations; ++k) {
    // 	cout << SumLocusAncestryHat[g*Populations + k] << "\t";
    //       }
    //     }
    // cout << "\t" << SumNumArrivals[0] << "\t" << SumNumArrivals[1] << "\t";
    SumNumArrivalsHat[0] /= (double)NumEstepiters;   
    SumNumArrivalsHat[1] /= (double)NumEstepiters; 
    
    //M-step: set mode for log rho and for theta in softmax basis - no minus ones in exponents
    // use vector rhohat and array ThetaProposal to store updates before accept/reject 
    unsigned gg = 0; // indexes which admixture prior (alpha) to use: use alpha[1] if second gamete and no indadmixhiermodel
    for(unsigned g = 0; g < NumIndGametes; ++g) {
      if( options->isAdmixed(g) ) {
	if(!options->isGlobalRho()) {
	  rhohat[g] = ( rhoalpha + SumNumArrivalsHat[g] ) / ( rhobeta + EffectiveL[g] );
	}
	if(g==1 && !options->getIndAdmixHierIndicator()) {
	  gg = 1;
	}
	double sum = accumulate(alpha[gg].begin(), alpha[gg].end(),0.0, std::plus<double>())
	  + accumulate(SumLocusAncestryHat+g*Populations, 
		       SumLocusAncestryHat+(g+1)*Populations, 0.0, std::plus<double>());
	for(int k = 0; k < Populations; ++k) {
	  if(alpha[gg][k]>0.0) ThetaProposal[g*Populations+k] = (alpha[gg][k]+SumLocusAncestryHat[g*Populations+k]) / sum;
	}
      } else {//unadmixed gamete - set ThetaProposal to (fixed) values in thetahat, which will be something like 1,0,0 
	copy(thetahat+g*Populations, thetahat+(g+1)*Populations, ThetaProposal+g*Populations);
      } 
    } // end loop over gametes
    double logpriorhat =  LogPriorTheta_Softmax(ThetaProposal, options, alpha) + 
      LogPriorRho_LogBasis(rhohat, options, rhoalpha, rhobeta);
    double loglikhat = getLogLikelihood(options, ThetaProposal, rhohat, false);
    double LogUnnormalizedPosteriorHat  = logpriorhat + loglikhat;
    if(LogUnnormalizedPosteriorHat > LogUnnormalizedPosterior) { //accept update only if density increases 
      if(isadmixed) {
	setAdmixtureProps(ThetaProposal, NumIndGametes * Populations);
	copy(rhohat.begin(), rhohat.end(), _rho.begin());
      }
      LogUnnormalizedPosterior = LogUnnormalizedPosteriorHat;
    } 
    cout << "LogPrior " << logpriorhat << "\tLogLikelihood " << loglikhat << "\tLogUnNormalizedPosterior " << 
      LogUnnormalizedPosterior << endl << flush;
  } //end EM outer loop
  delete[] SumLocusAncestryHat;
   
  //print values to file
  modefile<<setiosflags(ios::fixed)<<setprecision(3);
  modefile << myNumber << "\t";
  if(!options->isGlobalRho()) {
    cout << "rhohat\t"; 
    for(unsigned i = 0; i < NumIndGametes; ++i) {
      modefile<<_rho[i]<<"\t ";
      cout << _rho[i] << "\t";
    }
    cout << endl;
  }
  for(unsigned i = 0; i < NumIndGametes; ++i) { //loop over populations within gametes
    for(int k = 0; k < Populations; ++k) modefile<<Theta[i*Populations +k]<<"\t ";
  }
  
  if(myNumber==1 && options->getChibIndicator()){ // copy modes into hat arrays to use in Chib algorithm
    for(unsigned k = 0; k < Populations*NumIndGametes; ++k){
      thetahat[k] = Theta[k];
    }
    //check for zeros where not specified by prior
    for(unsigned i = 0; i < NumIndGametes; ++i) {
      unsigned gg = i; // index of alpha to use: for second gamete and no indadmixhiermodel use 1
      if(options->getIndAdmixHierIndicator()) gg = 0;
      double sum = 0.0;
      for(int k = 0; k < Populations; ++k) { // if prior not zero, and mode < 0.001, use 0.001
	if(alpha[gg][k]>0.0 && thetahat[i*Populations+k]<0.001) {
	    thetahat[i*Populations+k] = 0.001;
	}
      }
      for(int k = 0; k < Populations; ++k) {
	sum += thetahat[i*Populations+k];
      }
      cout << "thetahat" << i << " "; 
      for(int k = 0; k < Populations; ++k) { // re-normalize
	thetahat[i*Populations+k] /= sum;
	cout << thetahat[i*Populations+k] << " ";
      }
     } //end gamete loop
    cout << endl;
    copy(_rho.begin(), _rho.end(), rhohat.begin());
  }
}


void Individual::SampleTheta( int iteration, double *SumLogTheta, const DataMatrix* const Outcome, 
			      const DataType* const OutcomeType, const vector<double> lambda, int NumCovariates,
			      DataMatrix *Covariates, const vector<const double*> beta, const double* const poptheta,
			      const AdmixOptions* const options, const vector<vector<double> > &alpha, 
			      double DInvLink, double dispersion, bool RW, bool anneal=false)
// samples individual admixture proportions
// called with RW true for a random-walk proposal, false for a conjugate proposal
{
  double logpratio = 0.0;
  try{
    if(RW) {
      NumberOfUpdates++;
      logpratio += ProposeThetaWithRandomWalk(options, alpha); 
    } else ProposeTheta(options, alpha, SumLocusAncestry, SumLocusAncestry_X);       
  }
  catch(string s){
    string err = "Error encountered while generating proposal individual admixture proportions:\n" + s;
    throw err;
  }
  int K = Populations;

  //calculate Metropolis acceptance probability ratio for proposal theta    
  if(!options->getTestForAdmixtureAssociation()){
    RegressionType RegType;
    int NumOutcomes = Outcome->nCols();
    for( int k = 0; k < NumOutcomes; k++ ){
      if(OutcomeType[k] == Binary)RegType = Logistic; else RegType = Linear;
      logpratio +=  LogAcceptanceRatioForRegressionModel( RegType, options->isRandomMatingModel(), K, NumCovariates, 
							  Covariates, beta[k], //ExpectedY[ k ][myNumber-1], 
							  Outcome->get( myNumber-1, k ), poptheta, lambda[k]);
    }
  }
  
  //Accept or reject proposed value - if conjugate update and no regression model, proposal will be accepted because logpratio = 0
  Accept_Reject_Theta(logpratio, /*Loci->isX_data(), */ K, options->isRandomMatingModel(), RW );
  
  // update the value of admixture proportions used in the regression model  
  if( options->getNumberOfOutcomes() > 0 )
    UpdateAdmixtureForRegression(K, NumCovariates, poptheta, options->isRandomMatingModel(), Covariates);
  
  if(!anneal && iteration > options->getBurnIn()){ // accumulate sums in softmax basis for calculation of posterior means 
    
    for( unsigned int g = 0; g < NumIndGametes; g++ ){
      bool* b = new bool[Populations];
      double* a = new double[Populations];
      for(int k = 0; k < Populations; ++k)if(Theta[g*Populations + k] > 0.0){
	b[k] = true; //to skip elements set to zero
      } else b[k] = false;
      inv_softmax(Populations, Theta+g*Populations, a, b);
      transform(a, a+Populations, SumSoftmaxTheta+g*Populations, SumSoftmaxTheta+g*Populations, std::plus<double>());
	delete[] b;
	delete[] a;
    }
  }
  if(!IAmUnderTest){
    for( int k = 0; k < K; k++ ){
      SumLogTheta[ k ] += log( Theta[ k ] );
      if(options->isRandomMatingModel() )
	SumLogTheta[ k ] += log( Theta[ K + k ] );
    }
    
    //increment B using new Admixture Props
    //Xcov is a vector of admixture props as covariates as in UpdateScoreForAncestry
    if(iteration >= options->getBurnIn() && options->getTestForLinkageWithAncestry()){
      double* admixtureCovars = new double[Populations-1];
      for(int t = 0; t < Populations-1; ++t)admixtureCovars[t] = Covariates->get(myNumber-1, Covariates->nCols()-Populations+1+t);
      UpdateB(DInvLink, dispersion, admixtureCovars);
      delete[] admixtureCovars;
    }
  }
}

// ****** End Public Interface *******

double Individual::ProposeThetaWithRandomWalk(const AdmixOptions* const options, const vector<vector<double> > &alpha) {
  double LogLikelihoodRatio = 0.0;
  double LogPriorRatio = 0.0;
  
  //generate proposals
  for( unsigned int g = 0; g < NumIndGametes; g++ ){
    if(options->isAdmixed(g)){
      // inverse softmax transformation from proportions to numbers on real line that sum to 0
      bool* b = new bool[Populations];
      double* a = new double[Populations]; // should be at class scope
      for(int k = 0; k < Populations; ++k) {
	if(Theta[g*Populations + k] > 0.0) {
	  b[k] = true; //to skip elements set to zero
	} else b[k] = false;
      }
      inv_softmax(Populations, Theta+g*Populations, a, b);
      //random walk step - on all elements of array a
      for(int k = 0; k < Populations; ++k) {
	if( b[k] ) a[k] = Rand::gennor(a[k], step);  
      }
      //reverse transformation from numbers on real line to proportions 
      softmax(Populations, ThetaProposal+g*Populations, a, b);
      //compute contribution of this gamete to log prior ratio
      for(int k = 0; k < Populations; ++k) {
	if( b[k] ) { 
	  //	  LogPriorRatio += (alpha[g][k] - 1.0)*(log(ThetaProposal[g*Populations+k]) - 
	  //					 log(Theta[g*Populations+k])); 
	  // prior densities must be evaluated in softmax basis
	  LogPriorRatio += alpha[g][k]*(log(ThetaProposal[g*Populations+k]) - 
						 log(Theta[g*Populations+k])); 
	}
      }
      //     //compute contribution of this gamete to log prior ratio 
      //       LogPriorRatio += getDirichletLogDensity_Softmax(alpha[g], ThetaProposal+g*Populations) - 
      // 	getDirichletLogDensity_Softmax(alpha[g], Theta+g*Populations);
      delete[] a;
      delete[] b; 
    }
    else
      copy(Theta+g*Populations, Theta+(g+1)*Populations, ThetaProposal+g*Populations);
  }// end loop over gametes

  //get log likelihood at current parameter values - do not force update, store result of update
  LogLikelihoodRatio -= getLogLikelihood(options, false, true); 
  
  //get log likelihood at proposal theta and current rho - force update 
  // store result in loglikelihood.tempvalue, and accumulate loglikelihood ratio   
  logLikelihood.tempvalue = getLogLikelihood(options, ThetaProposal, _rho, true);
  LogLikelihoodRatio += logLikelihood.tempvalue;
  return LogLikelihoodRatio + LogPriorRatio;// log ratio of full conditionals
}

// Proposes new values for individual admixture proportions 
// as conjugate Dirichlet posterior conditional on prior parameter vector alpha and 
// multinomial likelihood given by sampled values of ancestry at loci where jump indicator xi is 1 (SumLocusAncestry)
// proposes new values for both gametes if random mating model 
void Individual::ProposeTheta(const AdmixOptions* const options, const vector<vector<double> > &alpha,
			      int* sumLocusAncestry, int* sumLocusAncestry_X){
  size_t K = Populations;
  if( options->isRandomMatingModel() ){ //random mating model
    for( unsigned int g = 0; g < 2; g++ ) {
      if(options->isAdmixed(g)) {
	unsigned gg = g; // index of which prior (alpha) to use - use alpha[1] if second gamete and no indadmixhiermodel
	if(options->getIndAdmixHierIndicator()) gg = 0;
	for(size_t k = 0; k < K; ++k) { // parameters of conjugate Dirichlet update
	  dirparams[k] = alpha[gg][k] + (double)(sumLocusAncestry[k + K*g] + sumLocusAncestry_X[g*K + k]);
	  // if male, paternal elements (0 to K-1) of sumLocusAncestry_X will be zero
	}
	//generate proposal theta from Dirichlet with parameters dirparams
	Rand::gendirichlet(K, dirparams, ThetaProposal+g*K );
      } else { // gamete unadmixed
	copy(Theta+g*Populations, Theta+(g+1)*Populations, ThetaProposal+g*Populations);
      }
    } // end loop over gametes
    
  } else { //assortative mating model
    for(size_t k = 0; k < K; ++k) {
      dirparams[k] = alpha[0][k] + double(sumLocusAncestry[k] + sumLocusAncestry_X[k] + 
					  sumLocusAncestry[k + K] + sumLocusAncestry_X[K + k]);
    }
    Rand::gendirichlet(K, dirparams, ThetaProposal );
  }
}

double Individual::LogAcceptanceRatioForRegressionModel( RegressionType RegType, bool RandomMatingModel, 
							 int Populations, int NumCovariates, 
							 const DataMatrix* const Covariates, const double* beta, 
							 const double Outcome, const double* const poptheta, const double lambda) {
  // returns log of ratio of likelihoods of new and old values of population admixture
  // in regression models.  individual admixture theta is standardized about the mean poptheta calculated during burn-in. 
  double logprobratio = 0.0, XBeta = 0.0, currentXBeta = 0.0;
  vector<double> avgtheta(Populations);avgtheta[0] = 0.0;
  vector<double> currentavgtheta(Populations);currentavgtheta[0] = 0.0;
  if( RandomMatingModel )
    for(int k = 1;k < Populations; ++k){
      avgtheta[k] = (ThetaProposal[k] + ThetaProposal[k + Populations ])/ 2.0 - poptheta[k];
      currentavgtheta[k] = (Theta[k] + Theta[k + Populations ])/ 2.0 - poptheta[k];
    }
  else
    for(int k = 1;k < Populations; ++k){
      avgtheta[k] = ThetaProposal[k]  - poptheta[k];
      currentavgtheta[k] = Theta[k]  - poptheta[k];
    }
  
  for( int jj = 0; jj < NumCovariates - Populations + 1; jj++ ){
    XBeta += Covariates->get( myNumber-1, jj ) * beta[jj];
    currentXBeta += Covariates->get( myNumber-1, jj ) * beta[jj];
  }
  for( int k = 1; k < Populations; k++ ){
    XBeta += avgtheta[ k ] * beta[NumCovariates - Populations + k ];
    currentXBeta += currentavgtheta[ k ] * beta[NumCovariates - Populations + k]; 
  }
  if(RegType == Linear) {
    logprobratio = 0.5 * lambda * (XBeta - currentXBeta) * (Outcome + Outcome - currentXBeta - XBeta);
    // ( currentXBeta - Outcome )^2 - ( XBeta - Outcome )^2 factorized 
  }
  else if(RegType == Logistic) {
    if( Outcome == 1 ) {
      logprobratio =  log( ( 1.0 + exp( -currentXBeta ) ) / ( 1.0 + exp( -XBeta ) ) );
    } else { 
      logprobratio =  log( ( 1.0 + exp( currentXBeta ) ) / ( 1.0 + exp( XBeta ) ) );
    }
  }
  return( logprobratio );
}

// update individual admixture values (mean of both gametes) used in the regression model
void Individual::UpdateAdmixtureForRegression( int Populations, int NumCovariates,
                                               const double* const poptheta, bool RandomMatingModel, 
					       DataMatrix *Covariates) {
  vector<double> avgtheta(Populations);
  if(RandomMatingModel )//average over gametes
    for(int k = 0; k < Populations; ++k) avgtheta[k] = (Theta[k] + Theta[k + Populations]) / 2.0;    
  else
    for(int k = 0; k < Populations; ++k) avgtheta[k] = Theta[k];
  for( int k = 1; k < Populations ; k++ )
    Covariates->set( myNumber-1, NumCovariates - Populations + k, avgtheta[ k ] - poptheta[ k ] );
}

void Individual::Accept_Reject_Theta( double logpratio, /*bool xdata, */ int Populations, bool RandomMatingModel, bool RW ) {
  // Metropolis update for admixture proportions theta, taking log of acceptance probability ratio as argument
  // uses log ratio because this is easier than ratio to calculate for linear regression model
  // if conjugate update and no regression model, logpratio remains set to 0, so all proposals are accepted 
  bool test = true;
  bool accept = false;
  double AccProb = 1.0; 
  // loop over populations: if any element of proposed Dirichlet parameter vector is too small, reject update without test step
  for( int k = 0; k < Populations; k++ ) {
    if( Theta[ k ] > 0.0 && ThetaProposal[ k ] < 0.0001 ) {
      test = false;
    }
    else if( RandomMatingModel && Theta[k + Populations] > 0.0 && ThetaProposal[ k + Populations ] < 0.0001 ) {
      test = false;
    }
  }

  if(test) { // generic Metropolis step
    if( logpratio < 0 ) {
      AccProb = exp(logpratio); 
      if( Rand::myrand() < AccProb ) accept=true;
    } else {
      accept = true;
    }
  }
  
  if(accept) { // set proposed values as new values    
    setAdmixtureProps(ThetaProposal, NumIndGametes * Populations);
    if(RW) { //if random-walk update, store the temp log-likelihood and set loglikelihood.HMMisOK to true
      storeLogLikelihood(true); 
    } else { // conjugate update of parameters invalidates both HMM forward probs and stored loglikelihood
      logLikelihood.HMMisOK = false;
      logLikelihood.ready = false;
    }
  } // if RW proposal is rejected, loglikelihood.HMMisOK is already set to false, and stored log-likelihood is still valid 

  if(RW) { //update step size in tuner object every w updates
    if( !( NumberOfUpdates % w ) ) {
      step = ThetaTuner.UpdateStepSize( AccProb );
    }
  }
}

void Individual::resetStepSizeApproximator(int k) {
  ThetaTuner.resetApproximator(k);
}

void Individual::UpdateHMMInputs(unsigned int j, const AdmixOptions* const options, 
				 const double* const theta, const vector<double> rho) {
  //Updates inputs to HMM for chromosome j
  //also sets Diploid flag in Chromosome (last arg of SetStateArrivalProbs)
  Chromosome* C = Loci->getChromosome(j);
  C->SetGenotypeProbs(GenotypeProbs[j], GenotypesMissing[j]);

  if(!options->getHapMixModelIndicator() && !options->isGlobalRho()){
    //set locus correlation, f, if individual- or gamete-specific rho
    C->SetLocusCorrelation(rho, !options->isRandomMatingModel(), options->isRandomMatingModel());
  }
  C->SetStateArrivalProbs(theta, options->isRandomMatingModel(), (j!=X_posn || SexIsFemale));
  logLikelihood.HMMisOK = false;//because forward probs in HMM have been changed
}

void Individual::SampleRho(const AdmixOptions* const options, double rhoalpha, double rhobeta, 
			   bool updateSumLogRho) {
  vector<unsigned> sumNumArrivals = getSumNumArrivals();
  vector<unsigned> sumNumArrivals_X = getSumNumArrivals_X();
  // rho_X is set to 0.5*rho - equivalent to setting effective length of X chromosome as half length in morgans
  // conjugate gamma update for rho includes arrivals on X chromosome, and 0.5 * length of X chromosome
  if(options->isRandomMatingModel() ) {
    // SumNumArrivals_X has length 2, and SumNumArrivals_X[0] remains fixed at 0 if male 
    for( unsigned int g = 0; g < 2; g++ )if(options->isAdmixed(g)) {
      //do {
      _rho[g] = Rand::gengam( rhoalpha + (double)(sumNumArrivals[g] + sumNumArrivals_X[g]), rhobeta + EffectiveL[g] );
      //} while( _rho[g] > TruncationPt || _rho[g] < 1.0 );
    } 
  } else {//assortative mating
    _rho[0] = Rand::gengam( rhoalpha + (double)(sumNumArrivals[0] + sumNumArrivals[1] + sumNumArrivals_X[0] + sumNumArrivals_X[1]), 
			rhobeta + EffectiveL[0] + EffectiveL[1] );
  }
  //now that rho has changed, current stored value of loglikelihood is no longer valid and 
  //HMMs will need to be updated before getting loglikelihood
  logLikelihood.HMMisOK = false;
  logLikelihood.ready = false;

  if(updateSumLogRho) { // accumulate log rho for calculation of posterior mean
    for(unsigned i = 0; i < _rho.size(); ++i) sumlogrho[i] += log(_rho[i]);
  }
}

void Individual::SampleMissingOutcomes(DataMatrix *Outcome, const DataType* const OutcomeType, 
				       const double* const* ExpectedY, const vector<double> lambda){
  int NumOutcomes = Outcome->nCols();
  // sample missing values of outcome variable
  for( int k = 0; k < NumOutcomes; k++ ){
    if( Outcome->isMissing( myNumber-1, k ) ){
      if( OutcomeType[k] == Continuous )
	Outcome->set( myNumber-1, k, Rand::gennor( ExpectedY[k][myNumber-1], 1 / sqrt( lambda[k] ) ));
      else{
	if( Rand::myrand() * ExpectedY[k][myNumber-1] < 1 )
	  Outcome->set( myNumber-1, k, 1);
	else
	  Outcome->set( myNumber-1, k, 0);
      }
    }
  }
}

//********************** Score Tests ***************************
void Individual::ResetScores(const AdmixOptions* const options){
  int KK = Populations;
  if(Populations == 2)KK = 1;

  if( options->getTestForAffectedsOnly() ){
    for(unsigned j = 0; j < Loci->GetNumberOfCompositeLoci()*KK; ++j){
      AffectedsScore[j] = 0.0;
      AffectedsVarScore[j] = 0.0;
      AffectedsInfo[j] = 0.0;
    }
  }
  if( options->getTestForLinkageWithAncestry() ){
    for(unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); ++j){
      for(int k = 0; k < 2*Populations; ++k)
	AncestryScore[j][k] = 0.0;
      for(int k = 0; k < 4*Populations*Populations; ++k)
	AncestryInfo[j][k] = 0.0;
      for(int k = 0; k < Populations; ++k){
	AncestryInfoCorrection[j][k] = 0.0;
	AncestryVarScore[j][k] = 0.0;
      }
    }
    for(int k = 0; k < Populations*Populations; ++k){
      PrevB[k] = B[k];           //PrevB stores the sum for the previous iteration
      B[k] = 0.0;                //while B accumulates the sum for the current iteration 
    }
    for(int k = 0; k < Populations; ++k){
      Xcov[k] = 0.0;
    }
  }
}

void Individual::UpdateScores(const AdmixOptions* const options, DataMatrix *Outcome, const DataType* const OutcomeType, 
			      DataMatrix *Covariates, double DInvLink, double dispersion,const double* const * ExpectedY){//merge with updatescoretests
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    Chromosome* C = Loci->getChromosome(j);
    // update of forward probs here is unnecessary if SampleTheta was called and proposal was accepted  
      //Update Forward/Backward probs in HMM
      if( !logLikelihood.HMMisOK ) {
	UpdateHMMInputs(j, options, Theta, _rho);
      }
      //update of score tests for linkage with ancestry requires update of backward probs
      double* admixtureCovars = 0;
      if(options->getTestForLinkageWithAncestry()){
	admixtureCovars = new double[Populations-1];
	for(int t = 0; t < Populations-1; ++t)admixtureCovars[t] = Covariates->get(myNumber-1, Covariates->nCols()-Populations+1+t);
      }
      UpdateScoreTests(options, admixtureCovars, Outcome, OutcomeType, C, DInvLink, dispersion, ExpectedY);
      delete[] admixtureCovars;
  } //end chromosome loop
}

void Individual::UpdateScoreTests(const AdmixOptions* const options, const double* admixtureCovars, DataMatrix *Outcome, 
				  const DataType* const OutcomeType,
				  Chromosome* chrm, double DInvLink, double dispersion, const double* const* ExpectedY){
  bool IamAffected = false;
  try {
    if( options->getTestForAffectedsOnly()){
      //determine which regression is logistic, in case of 2 outcomes
      unsigned col = 0;
      if(options->getNumberOfOutcomes() >1 && OutcomeType[0]!=Binary)col = 1;
      //check if this individual is affected
      if(options->getNumberOfOutcomes() == 0 || Outcome->get(myNumber-1, col) == 1) IamAffected = true;
    }
    
    //we don't bother computing scores for the first population when there are two
    int KK = Populations,k0 = 0;
    if(Populations == 2) {KK = 1;k0 = 1;}
    
    int locus;
    for( unsigned int jj = 0; jj < chrm->GetSize(); jj++ ){
      locus = chrm->GetLocus(jj); 
      //retrieve AncestryProbs from HMM
      std::vector<std::vector<double> > AProbs = 
	chrm->getAncestryProbs(!chrm->isXChromosome() || SexIsFemale,  jj );
      
      //Update affecteds only scores      
      if(IamAffected){
	UpdateScoreForLinkageAffectedsOnly(locus, KK, k0, options->isRandomMatingModel(), AProbs );
      }
      
      //update ancestry score tests
      if( options->getTestForLinkageWithAncestry() ){
	UpdateScoreForAncestry(locus, admixtureCovars, dispersion, Outcome->get(myNumber-1, 0) - ExpectedY[0][myNumber-1], 
			       DInvLink, AProbs);
      }
      ++locus;
    }//end within-chromosome loop
  } catch (string msg) {
    cout << msg;
  }
}

void Individual::UpdateScoreForLinkageAffectedsOnly(unsigned int locus, int Pops, int k0, bool RandomMatingModel, 
						    const vector<vector<double> > AProbs){
  // values of ancestry risk ratio at which likelihood ratio is evaluated
  double r1 = 0.5;
  double r2 = 2.0;//hard-coding these for now, can make them vary later
  if( SexIsFemale  || (Loci->GetChrNumOfLocus(locus) != X_posn) ) { // diploid case
    double theta[2];//paternal and maternal admixture proportions
    double Pi[3];//probs of 0,1,2 copies of Pop k given admixture
    for( int k = 0; k < Pops; k++ ){
      theta[0] = Theta[ k+k0 ];
      if( RandomMatingModel )
	theta[1] = Theta[ Populations + k+k0 ];
      else
	theta[1] = theta[0];
      
      //accumulate score, score variance, and info
      AffectedsScore[locus *Pops + k]+= 0.5*( AProbs[1][k+k0] + 2.0*AProbs[2][k+k0] - theta[0] - theta[1] );
      AffectedsVarScore[locus * Pops + k]+= 0.25 *( AProbs[1][k+k0]*(1.0 - AProbs[1][k+k0]) + 4.0*AProbs[2][k+k0]*AProbs[0][k+k0]); 
      AffectedsInfo[locus * Pops +k]+= 0.25* ( theta[0]*( 1.0 - theta[0] ) + theta[1]*( 1.0 - theta[1] ) );
      
      //probs of 0,1,2 copies of Pop k given admixture
      Pi[2] = theta[0] * theta[1];
      Pi[1] = theta[0] * (1.0 - theta[1]) + theta[1] * (1.0 - theta[0]);
      Pi[0] = (1.0 - theta[0]) * (1.0 - theta[1]);
      
      //compute contribution to likelihood ratio
      LikRatio1[locus *Pops + k] += (AProbs[0][k+k0] + sqrt(r1)*AProbs[1][k+k0] + r1 * AProbs[2][k+k0]) / 
	(Pi[0] + sqrt(r1)*Pi[1] + r1*Pi[2]);
      LikRatio2[locus *Pops + k] += (AProbs[0][k+k0] + sqrt(r2)*AProbs[1][k+k0] + r2 * AProbs[2][k+k0]) / 
	(Pi[0] + sqrt(r2)*Pi[1] + r2*Pi[2]);
    }
  } else { // haploid - effect of one extra copy from pop k0 is equivalent to two extra copies in diploid case 
    double theta;//paternal and maternal admixture proportions
    double Pi[2];//probs of 0,1 copies of Pop k given admixture
    for( int k = 0; k < Pops; k++ ){
      theta = Theta[ k+k0 ];
      
      //accumulate score, score variance, and info
      AffectedsScore[locus *Pops + k] += AProbs[1][k+k0] - theta;
      AffectedsVarScore[locus * Pops + k] += AProbs[0][k+k0] * AProbs[1][k+k0]; 
      AffectedsInfo[locus * Pops +k]+= theta * (1.0 - theta);
      
      //probs of 0,1 copies of Pop k given admixture
      Pi[1] = theta;
      Pi[0] = 1.0 - theta;
      
      //compute contribution to likelihood ratio - check this formula
      LikRatio1[locus *Pops + k] += (AProbs[0][k+k0] + r1*AProbs[1][k+k0]) / (Pi[0] + r1*Pi[1]);
      LikRatio2[locus *Pops + k] += (AProbs[0][k+k0] + r2*AProbs[1][k+k0]) / (Pi[0] + r2*Pi[1]);
    }
  }
}

void Individual::UpdateScoreForAncestry(int locus, const double* admixtureCovars, double phi, double YMinusEY, double DInvLink, 
					const vector<vector<double> > AProbs) {
  //Updates score stats for test for association with locus ancestry
  //now use Rao-Blackwellized estimator by replacing realized ancestries with their expectations
  //Notes: 1/phi is dispersion parameter
  //       = lambda[0] for linear regression, = 1 for logistic
  //       YMinusEY = Y - E(Y) = Y - g^{-1}(\eta_i)
  //       VarX = Var(X)
  //       DInvLink = {d  g^{-1}(\eta)} / d\eta = derivative of inverse-link function
  // admixtureCovars are the centred admixture proportions (except first one) used in regression model
  //Xcov is a vector of covariates
  //Note that only the intercept and admixture proportions are used.
  // X is (A, cov)'  
  
  double *X = new double [2 * Populations], *Xcopy = new double[2*Populations], *XX = new double[4*Populations*Populations];
  //Xcopy is an exact copy of X; We need two copies as one will be destroyed
  double xBx[1];
  double* BX = new double[Populations];
  double* VarA = new double[Populations];
 
  X[ 2*Populations - 1] = 1;//intercept
  Xcov[Populations-1] = 1;
  //set covariates, admixture props for pops 2 to K 
  for( int k = 0; k < Populations - 1; k++ ){
    X[ Populations + k] = admixtureCovars[k];//Theta[ k+1 ];
    BX[k] = Xcov[k] = admixtureCovars[k];//Theta[ k+1 ];
  }

  for( int k = 0; k < Populations ; k++ ){
    Xcopy[k] = X[k] = AProbs[1][k] + 2.0 * AProbs[2][k];//Conditional expectation of ancestry
    VarA[k] = AProbs[1][k]*(1.0 - AProbs[1][k]) + 4.0*AProbs[2][k]*AProbs[0][k];//conditional variances
  }
  //KLUDGE: need to reset Xcopy each time since destroyed in computation of score
  Xcopy[2*Populations-1] = 1;
  for( int k = 0; k < Populations-1; k++ )Xcopy[k + Populations] = admixtureCovars[k];//Theta[ k+1 ];

  try{
    // ** compute expectation of score **
    scale_matrix(Xcopy, YMinusEY*phi, 2*Populations, 1);      //Xcopy *= YMinusEY *phi
    add_matrix(AncestryScore[locus], Xcopy, 2*Populations, 1);//AncestryScore[locus] += Xcopy
    // ** compute uncorrected info **
    matrix_product(X, X, XX, 2*Populations, 1, 2*Populations);        //XX = X'X
    scale_matrix(XX, DInvLink*phi, 2*Populations, 2*Populations);     //XX = DInvLink * phi * X'X
    add_matrix(AncestryInfo[locus], XX, 2*Populations, 2*Populations);//AncestryInfo[locus] += XX
    // ** compute variance of score and correction term for info **    
    HH_solve(Populations, PrevB, Xcov, BX);          //BX = inv(PrevB) * Xcov
    matrix_product(Xcov, BX, xBx, 1, Populations, 1);//xBx = Xcov' * BX
  }
  catch(string s){
    delete[] X; delete[] Xcopy;
    delete[] XX;
    delete[] BX;
    delete[] VarA;
    std::string error_string = "Error in Individual::UpdateScoreForAncestry:\n";
    error_string.append(s);
    throw(error_string);//throw error message to top level
  }
  for( int k = 0; k < Populations ; k++ ){
    AncestryInfoCorrection[locus][k] += VarA[k] * (DInvLink *phi - phi * phi * DInvLink * DInvLink * xBx[0]); 
    AncestryVarScore[locus][k] += VarA[k] * phi * phi * YMinusEY * YMinusEY;
  }
  delete[] X; delete[] Xcopy;
  delete[] XX;
  delete[] BX;
  delete[] VarA;
}

void Individual::UpdateB(double DInvLink, double dispersion, const double* admixtureCovars){
  //increment B using new Admixture Props
  //Xcov is a vector of admixture props as covariates as in UpdateScoreForAncestry
    Xcov[Populations-1] = 1;//last entry is intercept
    for( int k = 0; k < Populations - 1; k++ ){
      Xcov[k] = admixtureCovars[k]; //centred admixture covariates
    }
    double *temp = new double[Populations*Populations];
    matrix_product(Xcov, Xcov, temp, Populations, 1, Populations);
    scale_matrix(temp, DInvLink*dispersion, Populations, Populations);
    add_matrix(B, temp, Populations, Populations);
    delete[] temp;
}

void Individual::SumScoresForLinkageAffectedsOnly(int j, double *SumAffectedsScore, 
				      double *SumAffectedsVarScore, double *SumAffectedsScore2, double *SumAffectedsInfo){
  int KK = Populations;
  if(KK == 2) KK = 1;

  for( int k = 0; k < KK; k++ ){
    SumAffectedsScore[j*KK +k] += AffectedsScore[j*KK + k];
    SumAffectedsVarScore[j*KK +k] += AffectedsVarScore[j * KK +k];
    SumAffectedsInfo[j*KK +k] += AffectedsInfo[j * KK +k];
    SumAffectedsScore2[j*KK +k] +=  AffectedsScore[j*KK +k] * AffectedsScore[j*KK +k];
  }
}

void Individual::SumScoresForAncestry(int j, double *SumAncestryScore, double *SumAncestryInfo, double *SumAncestryScore2,
				      double *SumAncestryVarScore){

  double *score = new double[Populations], *info = new double[Populations*Populations];
  try{
    CentredGaussianConditional(Populations, AncestryScore[j], AncestryInfo[j], score, info, 2*Populations );
  }
  catch(string s){
    string error = "Error centring ancestry association scores in Individual:\n";
    error.append(s);
    throw(error);
  }
  
  //accumulate over iterations
  //for two populations, accumulate the scores for second population only
  int KK = Populations, k1 = 0;
  if(Populations == 2){KK = 1; k1 = 1;}
     
  for( int k = 0; k < KK ; k++ ){
    SumAncestryScore[j*KK +k] += score[k+k1];
    SumAncestryInfo[j*KK +k] += info[(k+k1)*Populations +k+k1] + AncestryInfoCorrection[j][k+k1];
    SumAncestryScore2[j*KK +k] += score[k+k1] * score[k+k1];
    SumAncestryVarScore[j*KK +k] += AncestryVarScore[j][k+k1];
  }
  delete[] score;
  delete[] info;
}

// this function does three things:
// 1. sets allelefreqsMAP to current values of allelefreqs
// 2. calculates log-likelihood and log prior at thetahat, rhohat, allelefreqsMAP
void Individual::setChibNumerator(const AdmixOptions* const options, const vector<vector<double> > &alpha, 
				  double rhoalpha, double rhobeta, chib *MargLikelihood, AlleleFreqs* A) {
  
  // 1. set allelefreqsMAP in AlleleFreqs object
  if(A->IsRandom() ) {
    A->setAlleleFreqsMAP(); 
    /** does three things: 
	1. allocates array for AlleleFreqsMAP
	2. sets elements of array to current value
	3. loops over composite loci to set AlleleProbsMAP to values in AlleleFreqsMAP
    **/

    //   do we need to set HapPairProbs in composite loci 

    //set HapPairProbsMAP to current values of HapPairProbs
    for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
      (*Loci)(j)->setHapPairProbsMAP(); // sets HapPairProbsMAP to HapPairProbs
    }
    
    // now set genotype probs using HapPairProbsMAP and AlleleProbsMAP 
    for(unsigned j = 0; j < Loci->GetNumberOfChromosomes(); ++j){
      unsigned locus = Loci->getChromosome(j)->GetLocus(0);
      for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ) {
  	SetGenotypeProbs(j, jj, locus, true);  
      }
      ++locus;
    }
  }    
  
  // 2. calculate log-likelihood at MAP parameter values
  MargLikelihood->setLogLikelihood(getLogLikelihood( options, thetahat, rhohat, true));
  
  // 3. calculate log prior at MAP parameter values
  double LogPrior = LogPriorTheta_Softmax(thetahat, options, alpha)
    + LogPriorRho_LogBasis(rhohat, options, rhoalpha, rhobeta);
  if( A->IsRandom() ) {
    double LogPriorFreqs = 0.0;
    for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
      for( int k = 0; k < Populations; k++ ){
	vector<double> args = A->GetPriorAlleleFreqs(j,k);
	LogPriorFreqs += getDirichletLogDensity( A->GetPriorAlleleFreqs(j, k), A->getAlleleFreqsMAP(j,k) );
      }
    }
    LogPrior += LogPriorFreqs;
  }
  MargLikelihood->addLogPrior(LogPrior);
} 

void Individual::updateChib(const AdmixOptions* const options, const vector<vector<double> > &alpha, 
			    double rhoalpha, double rhobeta, 
			    // double *thetahat, vector<double> &rhohat, 
			    chib *MargLikelihood, AlleleFreqs* A){
  // *** After BurnIn *** - accumulate samples of posterior ordinates for theta, rho, allelefreqs separately
  double LogPosterior = 0.0;
  double LP = 0.0;
  if( Populations > 1 ){
    LP = LogPosteriorTheta_Softmax(options, thetahat, alpha);
    logPosterior[0].push_back(LP);
    LogPosterior += LP;
    LP = LogPosteriorRho_LogBasis(options, rhohat, rhoalpha, rhobeta);
    logPosterior[1].push_back(LP);
    LogPosterior += LP;
  }
  if( A->IsRandom() ){
    double LogPosteriorFreqs = 0.0;
    for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
      for( int k = 0; k < Populations; k++ ){
	vector<double> args = A->GetPriorAlleleFreqs(j,k);
	vector<int> counts = A->GetAlleleCounts(j,k);
	transform(counts.begin(), counts.end(), args.begin(), args.begin(), plus<double>());//PriorAlleleFreqs + AlleleCounts
	LogPosteriorFreqs += getDirichletLogDensity( args, A->getAlleleFreqsMAP(j, k));//LogPosterior for Allele Freqs
      }
    }
    LogPosterior += LogPosteriorFreqs;
    logPosterior[2].push_back( LogPosteriorFreqs  );
  }
  MargLikelihood->addLogPosteriorObs( LogPosterior );
}

// TODO: fix these two functions to work with assortative mating
double Individual::LogPriorTheta_Softmax(const double* const theta, const AdmixOptions* const options, 
					 const vector<vector<double> > &alpha) const {
  double LogPrior=0.0;
  for(int g = 0; g < 2; ++g) { //loop over gametes
    if( options->isAdmixed(g) ){
      LogPrior += getDirichletLogDensity_Softmax( alpha[g], theta + g*Populations);
    }
  }
  return LogPrior;
}

double Individual::LogPosteriorTheta_Softmax(const AdmixOptions* const options, const double* const theta, 
					     const vector<vector<double> > &alpha) const{
  // calculates log full conditional at theta, conditional on realized locus ancestry states and jump indicators
  double LogPosterior = 0.0;
  vector<double> alphaparams(Populations); // , alphaparams1(Populations);  // to be set to alpha + SumLocusAncestry
  for(int g = 0; g < 2; ++g) { //loop over gametes
    if( options->isAdmixed(g) ) {
      transform(alpha[g].begin(), alpha[g].end(), SumLocusAncestry + g*Populations, 
		alphaparams.begin(), std::plus<double>());
      LogPosterior += getDirichletLogDensity_Softmax(alphaparams, theta + g*Populations);
    }
  }
  //   if(  options->isAdmixed(1) ){//admixed second gamete
  //     transform(alpha[1].begin(), alpha[1].end(), SumLocusAncestry+Populations, alphaparams1.begin(), std::plus<double>());
  //     LogPosterior += getDirichletLogDensity_Softmax(alphaparams1, theta+Populations);
  //   }
  return LogPosterior;
}

double Individual::LogPriorRho_LogBasis(const vector<double> rho, const AdmixOptions* const options, double rhoalpha, 
					double rhobeta) const {
  // computes log prior density in log rho basis at supplied parameter values
  double LogPrior=0.0;
  if( options->isRandomMatingModel() ) {
    for(int g = 0; g < 2; ++g) { //loop over gametes
      if( options->isAdmixed(g) ){
	LogPrior += getGammaLogDensity_LogBasis( rhoalpha, rhobeta, rho[g] );
      }
    }
  } else { // assortative mating: rho assumed same on both gametes
    LogPrior = getGammaLogDensity_LogBasis( rhoalpha, rhobeta, rho[0] );
  }
  return LogPrior;
}

double Individual::LogPosteriorRho_LogBasis(const AdmixOptions* const options, const vector<double> rho, 
					    double rhoalpha, double rhobeta)const{
  // calculates log full conditional density at sum-intensities rho, conditional on realized number of arrivals
  double LogPosterior = 0.0;
  vector<unsigned> SumN = getSumNumArrivals();
  vector<unsigned> SumN_X = getSumNumArrivals_X();
  if(options->isRandomMatingModel() ) { // SumNumArrivals_X has length 2, and SumNumArrivals_X[0] remains fixed at 0 if male 
    for( unsigned int g = 0; g < 2; g++ ) {
      if(options->isAdmixed(g)) {
	LogPosterior += getGammaLogDensity_LogBasis( rhoalpha + (double)(SumN[g] + SumN_X[g]), rhobeta + EffectiveL[g], rho[g] );
      } 
    }
  } else {//assortative mating, rho assumed same on both gametes
    if(options->isAdmixed(0)) {
      LogPosterior+= getGammaLogDensity_LogBasis( rhoalpha + (double)(SumN[0] + SumN[1] + SumN_X[0] + SumN_X[1]), 
					 rhobeta + EffectiveL[0] + EffectiveL[1], rho[0] );
    }
  }
  return LogPosterior;
}

// these three functions can be used to monitor stability of estimates of components of posterior ordinate
double Individual::getLogPosteriorTheta()const{
  if(Populations > 1){
    std::vector<double>::const_iterator max = max_element(logPosterior[0].begin(), logPosterior[0].end());
    return AverageOfLogs(logPosterior[0], *max);
  }
  else return 0.0;
}
double Individual::getLogPosteriorRho()const{
  if(Populations > 1){
    std::vector<double>::const_iterator max = max_element(logPosterior[1].begin(), logPosterior[1].end());
    return AverageOfLogs(logPosterior[1], *max);
  }
  else return 0.0;
}
double Individual::getLogPosteriorAlleleFreqs()const{
  std::vector<double>::const_iterator max = max_element(logPosterior[2].begin(), logPosterior[2].end());
  return AverageOfLogs(logPosterior[2], *max);
}

//******************* Likelihood Ratios (Affecteds-only score test) ************************

static string double2R( double x )
{
  if( isnan(x) )
    return "NaN";
  else{
    stringstream ret;
    ret << x;
    return( ret.str() );
  }
}

void Individual::OutputLikRatios(const char* const filename, int iterations, const std::string* const PopLabels){
  //open outut file
  std::ofstream outputstream(filename);

  //start writing R object
  outputstream<< "structure(.Data=c(" << endl;

  //output ergodic averages of LikRatios
  int KK = Populations, k1 = 0;
  if(KK == 2 ){
    KK = 1;k1 = 1;
  }
  double L1, L2;

  for(unsigned int j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
    for( int k = 0; k < KK; k++ ){//end at 1 for 2pops
      outputstream << "\"" << (*Loci)(j)->GetLabel(0) << "\",";
      outputstream << "\""<<PopLabels[k+k1] << "\","; //need offset to get second poplabel for 2pops
      
      L1 = LikRatio1[ j*KK + k] / ( iterations );
      L2 = LikRatio2[ j*KK + k] / ( iterations );
      
      outputstream << double2R(L1)<< ","
		   << double2R(L2)<< ","<<endl;
    }
  }
  
  vector<int> dim(2,0);
  dim[0] = 4;
  dim[1] = Loci->GetNumberOfCompositeLoci() * KK;
  
  vector<string> labels(4,"");
  labels[0] = "Locus";
  labels[1] = "Population";
  labels[2] = "L1";
  labels[3] = "L2";
  
  outputstream << ")," << endl;
  outputstream << ".Dim = c(";
  for(unsigned int i=0;i<dim.size();i++){
    outputstream << dim[i];
    if(i != dim.size() - 1){
      outputstream << ",";
    }
  }
  outputstream << ")," << endl;
  outputstream << ".Dimnames=list(c(";
  for(unsigned int i=0;i<labels.size();i++){
    outputstream << "\"" << labels[i] << "\"";
    if(i != labels.size() - 1){
      outputstream << ",";
    }
  }
  outputstream << "), character(0)))" << endl;
  outputstream.close();
}
