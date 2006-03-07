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

unsigned int Individual::numChromosomes;
Genome *Individual::Loci;
int Individual::Populations;

//******** Constructors **********
Individual::Individual() {//should initialise pointers here
}

Individual::Individual(int number, const AdmixOptions* const options, const InputData* const Data,  
		       bool undertest=false) {
  myNumber = number;
  IAmUnderTest = undertest;
  NumIndGametes = 1;
  if( options->isRandomMatingModel() ) NumIndGametes = 2;

  double init=0.0;
  TruncationPt = options->getTruncPt(); 
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
  _rho.assign(NumIndGametes, init); // for global rho, Latent should assign initial value 
  _rho_X.assign(NumIndGametes, 0.5*init); // set to half the autosomal value 
  
  sumlogrho.assign(_rho.size(), 0.0);
  
  // Read sex value if present.
  sex = male;
  if (options->getgenotypesSexColumn() == 1) {
    sex = Data->GetSexValue(myNumber);
  }
  // this bool will be used to replace the sex struct
  SexIsFemale = false;
  if(sex == female) SexIsFemale = true;
  Xdata = Loci->isX_data();
  
  int numCompositeLoci = Loci->GetNumberOfCompositeLoci();
  Theta = 0;
  ThetaX = 0;
  ThetaProposal = 0;
  ThetaXProposal = 0;
  SumSoftmaxTheta = 0;

  // SumLocusAncestry is sum of locus ancestry states over loci at which jump indicator xi is 1  
  SumLocusAncestry = new int[options->getPopulations()*2];

  SumNumArrivals.resize(2*numCompositeLoci);  

  dirparams = new double[Populations]; //to hold dirichlet parameters for conjugate updates of theta
  
  ThetaProposal = new double[ Populations * NumIndGametes ];
  Theta = new double[ Populations * NumIndGametes ];
  SumSoftmaxTheta = new double[ Populations * NumIndGametes ];
  fill(SumSoftmaxTheta, SumSoftmaxTheta + Populations*NumIndGametes, 0.0);
  
  // X chromosome objects
  SumLocusAncestry_X = 0;    
  // if(Loci->isX_data() ){
  if( !SexIsFemale ){
    ThetaXProposal = new double[ Populations];
    ThetaX = new double[Populations];
    // ThetaXHat = new double[Populations];
    SumLocusAncestry_X = new int[Populations]; 
  } else {
    ThetaXProposal = new double[ Populations * 2 ];
    ThetaX = new double[ Populations * 2 ];
    //ThetaXHat = new double[ Populations * 2 ];
    SumLocusAncestry_X = new int[Populations * 2 ];
  }

  if(options->getHapMixModelIndicator()){
    SetDefaultAdmixtureProps();
  }
    
  // vector of possible haplotype pairs - 2 integers per locus if diploid, 1 if haploid 
  PossibleHapPairs = new vector<hapPair>[numCompositeLoci];
  
  LocusAncestry = new int*[ numChromosomes ]; // array of matrices in which each col stores 2 integers   
  X_posn = 9999; //position of the X chromosome in the sequence of chromosomes in the input data
  size_t AncestrySize = 0;  // set size of locus ancestry array
  //gametes holds the number of gametes for each chromosome, either 1 or 2
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    if( !Loci->isXChromosome(j) ){ // if not X chromosome, set number of elements to 2 * num loci
      AncestrySize = 2 * Loci->GetSizeOfChromosome(j) ;
      gametes.push_back(2);
    } else if( !SexIsFemale ) {//male or missing
      AncestrySize = Loci->GetSizeOfChromosome(j) ;
      gametes.push_back(1);
      X_posn = j;
    } else{//female
      AncestrySize = 2 * Loci->GetSizeOfChromosome(j) ;
      gametes.push_back(2);
      X_posn = j;
    } LocusAncestry[j] = new int[ AncestrySize];
    for(unsigned i = 0; i < AncestrySize; ++i) LocusAncestry[j][i] = 0;
  }
  //retrieve genotypes
  Data->GetGenotype(myNumber, options->getgenotypesSexColumn(), *Loci, &genotypes);
  // loop over composite loci to set possible haplotype pairs compatible with genotype 
  for(int j = 0; j < numCompositeLoci; ++j) {
    (*Loci)(j)->setPossibleHaplotypePairs(genotypes[j].alleles, PossibleHapPairs[j]);
    hapPair h;
    sampledHapPairs.push_back(h);
  }
  //initialise genotype probs array and array of indicators for genotypes missing at locus
  GenotypeProbs = new double*[numChromosomes];
  GenotypesMissing = new bool*[numChromosomes];
  for(unsigned j = 0; j < numChromosomes; ++j) {
    GenotypeProbs[j] = new double[Loci->GetSizeOfChromosome(j)*Populations*Populations];
    GenotypesMissing[j] = new bool[ Loci->GetSizeOfChromosome(j) ];
    SetGenotypeProbs(j, false);
    int nlocus = Loci->getChromosome(j)->GetLocus(0);  
    for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ) {
      if( !(GenotypeIsMissing(nlocus)) ) GenotypesMissing[j][jj] = false;
      else GenotypesMissing[j][jj] = true;
      nlocus++;
    }
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
  delete[] GenotypeProbs; //TODO: delete properly
  delete[] GenotypesMissing;
  delete[] LocusAncestry;
  delete[] SumLocusAncestry;
  delete[] SumLocusAncestry_X;
  delete[] dirparams;
  delete[] Theta;
  delete[] ThetaX;
  delete[] ThetaProposal;
  delete[] ThetaXProposal;
  delete[] SumSoftmaxTheta;
}

void Individual::drawInitialAdmixtureProps(const std::vector<std::vector<double> > &alpha) {
  // draw initial values for admixture proportions theta from Dirichlet prior 
  size_t K = Populations;
  for( unsigned g = 0; g < NumIndGametes; ++g ) { 
    double sum = 0.0;
    for(size_t k = 0; k < K; ++k) { //loop over array of Dirichlet params
      dirparams[k] = alpha[g][k]; // default value if elements of alpha sum to <= 1
      sum += alpha[g][k];
    }
   for(size_t k = 0; k < K; ++k) { //loop over array of Dirichlet params
      Theta[g*K+k] = alpha[g][k] / sum; // default value if elements of alpha sum to <= 1
      // if(alpha[g][k]>0.0) ++sum;
   }
    //generate proposal theta from Dirichlet with parameters dirparams
    if(sum>1)gendirichlet(K, dirparams, Theta+g*K );
  } // end loop over gametes: 
  if(Xdata) 
    setAdmixturePropsX(Theta, K*NumIndGametes);
}
void Individual::SetDefaultAdmixtureProps(){
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

void Individual::SetGenotypeProbs(int j, bool chibindicator=false){
  //chibindicator is passed to CompositeLocus object.  If set to true, CompositeLocus will use HapPairProbsMAP
  //instead of HapPairProbs when allelefreqs are not fixed.
  int locus = Loci->getChromosome(j)->GetLocus(0);
  for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ){
    if( !(GenotypeIsMissing(locus)) ){
      (*Loci)(locus)->GetGenotypeProbs(GenotypeProbs[j]+jj*Populations*Populations, PossibleHapPairs[locus], 
				       chibindicator);
    } else {
      for( int k = 0; k < Populations*Populations; ++k ) GenotypeProbs[j][jj*Populations*Populations + k] = 1.0;
    }
    locus++;
  }
}

void Individual::AnnealGenotypeProbs(int j, const double coolness) {
  // called after energy has been evaluated, before updating model parameters
  int locus = Loci->getChromosome(j)->GetLocus(0);
  for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ){ // loop over composite loci
    if( !(GenotypeIsMissing(locus)) ) { 
      for(int k = 0; k < Populations*Populations; ++k) // loop over ancestry states
	GenotypeProbs[j][jj*Populations*Populations+k] = pow(GenotypeProbs[j][jj*Populations*Populations+k], coolness); 
    }
    locus++;
  }
}

void Individual::setGenotypesToMissing(){
  for(unsigned i = 0; i < genotypes.size(); ++i)genotypes[i].missing=true;
}

//********** Allocation and deletion of static objects for score tests
void Individual::SetStaticMembers(Genome* const pLoci, const AdmixOptions* const options){
  Loci = pLoci;
  numChromosomes = Loci->GetNumberOfChromosomes();
  Populations = options->getPopulations();
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
  if(Xdata) 
    setAdmixturePropsX(Theta, size);
}

void Individual::setAdmixturePropsX(const double* const a, size_t size) {
  if(SexIsFemale || size == (unsigned)Populations) { // if female or size_t=K, assign theta_X[i] = theta[i] 
    for(unsigned i = 0; i < size; ++i)  {
      ThetaX[i] = a[i];
    }
  } else { // if male and size_t = 2K, assign thetaX = maternal gamete admixture props
    for(unsigned i = 0; i < (unsigned)Populations; ++i) {
      ThetaX[i] = a[Populations + i];
    }
  }
}

void Individual::HMMIsBad(bool loglikisbad) {
  logLikelihood.HMMisOK = false;
  if(loglikisbad)logLikelihood.ready = false;
}

//******************** Accessors ***********************************************************
const vector<vector<unsigned short> > Individual::getGenotype(unsigned int locus)const{
  return genotypes[locus].alleles;
}

const std::vector<hapPair > &Individual::getPossibleHapPairs(unsigned int locus)const{
  return PossibleHapPairs[locus];
}

// returns an indicator for non-missing genotype at simple locus, and updates allele count array  
bool Individual::GetAlleleCountsAtLocus(int complocus, int locus, int allele, int* count)const
{
  int notmissing = true;
  if(genotypes[complocus].alleles[locus][0]==0 || genotypes[complocus].alleles[locus][1]==0)
    notmissing = false;
  else{
    (*count) = (genotypes[complocus].alleles[locus][0] == allele) + (genotypes[complocus].alleles[locus][1] == allele);
  }
  return notmissing;
}

const int* Individual::getSampledHapPair(int locus)const{
  return sampledHapPairs[locus].haps;
}

const double* Individual::getAdmixtureProps()const {
  return Theta;
}

Sex Individual::getSex()const {
   return sex;
}

double Individual::getSumrho()const { //returns sum of sumintensities over gametes
  double sumrho = 0;
  for( unsigned int g = 0; g < _rho.size(); g++ )
    sumrho += _rho[g];
  return sumrho;
}

const vector<double> Individual::getRho()const {
   return _rho;
}

void Individual::GetLocusAncestry(int locus, int Ancestry[2])const{
  const vector<int> locustable = Loci->GetChrmAndLocus(locus);
  GetLocusAncestry(locustable[0], locustable[1], Ancestry);
}
void Individual::GetLocusAncestry(int chrm, int locus, int Ancestry[2])const {
  Ancestry[0]  = LocusAncestry[chrm][locus];
  if((unsigned)chrm == X_posn)Ancestry[1] = Ancestry[0];
  else Ancestry[1] = LocusAncestry[chrm][Loci->GetSizesOfChromosomes()[chrm]  + locus];
}

//returns value of LocusAncestry at a locus for a particular gamete
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
const vector<unsigned>Individual::getSumNumArrivals()const{
  //returns number of arrivals on each across genome
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

//Indicates whether genotype is missing at all simple loci within a composite locus
bool Individual::GenotypeIsMissing(unsigned int locus)const {
  return genotypes[locus].missing;
}

//****************** Log-Likelihoods **********************
// public method: 
// calls private method to get log-likelihood at current parameter values, and stores it either as loglikelihood.value or as loglikelihood.tempvalue
// store should be false when calculating energy for an annealed run, or when evaluating proposal for global sum-intensities
double Individual::getLogLikelihood( const AdmixOptions* const options, const bool forceUpdate, const bool store) {

  if (!logLikelihood.ready || forceUpdate) {
    logLikelihood.tempvalue = getLogLikelihood(options, Theta, ThetaX, _rho, _rho_X, true);
    if(store) {  
      logLikelihood.value = logLikelihood.tempvalue; 
      logLikelihood.ready = false; //true;
      logLikelihood.HMMisOK = false; //true; //because forward probs now correspond to current parameter values 
    }                               //and call to UpdateHMMForwardProbs has set this to false
    return logLikelihood.tempvalue;   
  } else return logLikelihood.value; // nothing was changed
}

// private method: gets log-likelihood at parameter values specified as arguments, but does not update loglikelihoodstruct
// arguments rho and rho_X are ignored if globalrho
double Individual::getLogLikelihood(const AdmixOptions* const options, const double* const theta, 
				    const double* const thetaX, const vector<double > rho, const vector<double> rho_X, 
				    bool updateHMM) {
  double LogLikelihood = 0.0;
  if(Populations == 1) LogLikelihood = getLogLikelihoodOnePop();
  else { 
    //   if(myNumber==1) cout << "\nprivate method: theta[0] " << theta[0] << " " << theta[1] << " thetaX[0] " << thetaX[0] << " " << thetaX[1] 
    //                        << " rho[0] " << rho[0] << " rhoX[0] " << rho_X[0] << endl; 
    for( unsigned int j = 0; j < numChromosomes; j++ ){
      if(updateHMM){// force update of forward probs 
	UpdateHMMForwardProbs(j, options, theta, thetaX, rho, rho_X);
      }
      LogLikelihood += Loci->getChromosome(j)->getLogLikelihood( !Loci->isXChromosome(j) || SexIsFemale );
      // if(myNumber==1) cout << "\n " << j << " " << LogLikelihood << endl;
    }
  }
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
  //TODO: X chromosome objects
  //obtain ergodic averages of (inv_softmax)admixture props and (log)sumintensities and transform back
  //to original scales
  unsigned size = Populations * NumIndGametes;
  for(unsigned i = 0; i < size; ++i)SumSoftmaxTheta[i] /= (options->getTotalSamples() - options->getBurnIn());
  for(unsigned i = 0; i < _rho.size(); ++i)sumlogrho[i] = exp(sumlogrho[i]/(options->getTotalSamples() - options->getBurnIn()));

  //apply softmax transformation to obtain thetabar
  double ThetaBar[NumIndGametes*Populations];
  for( unsigned int g = 0; g < NumIndGametes; g++ ){
    bool b[Populations];
    for(int k = 0; k < Populations; ++k)if(Theta[g*Populations + k] > 0.0){
      b[k] = true; //to skip elements set to zero
    } else b[k] = false;
    softmax(Populations, ThetaBar+g*Populations, SumSoftmaxTheta+g*Populations, b);
  }

  double LogLikelihood = 0.0;
  if(Populations == 1) LogLikelihood = getLogLikelihoodOnePop();
  else {
    for( unsigned int j = 0; j < numChromosomes; j++ ) {
      UpdateHMMForwardProbs(j, options, ThetaBar, ThetaBar, sumlogrho, sumlogrho); // ?sumlogrho is posterior mean of rho
      // should replace 2nd sumlogrho by half its value
      LogLikelihood += Loci->getChromosome(j)->getLogLikelihood( !Loci->isXChromosome(j) || SexIsFemale );
    }
  }
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
	UpdateHMMForwardProbs(j, options, Theta, ThetaX, _rho, _rho_X);
      }
      // sampling locus ancestry can use current values of forward probability vectors alpha in HMM 
      C->SampleLocusAncestry(LocusAncestry[j], Theta);
  } //end chromosome loop
}

void Individual::SampleHapPair(AlleleFreqs *A){
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    Chromosome* C = Loci->getChromosome(j);
    //loop over loci on current chromosome and update allele counts
    for( unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ){
      int locus =  C->GetLocus(jj);
      if( !(GenotypeIsMissing(locus)) ){
	int anc[2];//to store ancestry states
	GetLocusAncestry(j,jj,anc);
	//might be a shortcut for haploid data since there is only one compatible hap pair, no need to sample
	(*Loci)(locus)->SampleHapPair(sampledHapPairs[locus].haps, PossibleHapPairs[locus], anc);
	A->UpdateAlleleCounts(locus, sampledHapPairs[locus].haps, anc, C->isDiploid());
      }
    }   
  } //end chromosome loop
}

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

void Individual::FindPosteriorModes(const AdmixOptions* const options, const vector<vector<double> > &alpha,  
				    double rhoalpha, double rhobeta, //const vector<double> sigma, 
				    ofstream &modefile,
				    double *thetahat, double *thetahatX, vector<double> &rhohat, vector<double> &rhohatX){

  //uses an EM algorithm to search for posterior modes of individual parameters
  unsigned numEMiters = 10;
  unsigned NumEstepiters = 10; 
  double LogUnnormalizedPosterior = - numeric_limits<double>::max( );  

  //use current parameter values as initial values
  double *SumLocusAncestryHat = new double[2*Populations];
  int *SumLocusAncestry_XHat = 0;
  for(unsigned EMiter = 0; EMiter < numEMiters; ++EMiter){
    double SumNumArrivalsHat[2] = {0,0}, SumNumArrivals_XHat[2] = {0,0};
    fill(SumLocusAncestryHat, SumLocusAncestryHat + 2*Populations, 0.0);
    //TODO: SumLocusAncestry_XHat

    //E-step: fix theta and rho, sample Locus Ancestry and Number of Arrivals
    NumEstepiters *= 2;
    for(unsigned Estepiters = 0; Estepiters < (unsigned)NumEstepiters ; ++Estepiters) {
 
      if(Populations >1){
	ResetSufficientStats();
	SampleLocusAncestry(options);
	SampleJumpIndicators((!options->isGlobalRho() || options->getHapMixModelIndicator()));
      }

      vector<unsigned> SumN = getSumNumArrivals();
      vector<unsigned> SumN_X = getSumNumArrivals_X();

      // accumulate sums
      SumNumArrivalsHat[0] += SumN[0];	SumNumArrivals_XHat[0] += SumN_X[0];
      SumNumArrivalsHat[1] += SumN[1];	SumNumArrivals_XHat[1] += SumN_X[1];
      transform(SumLocusAncestry, SumLocusAncestry+2*Populations, 
		SumLocusAncestryHat, SumLocusAncestryHat, std::plus<double>());
      //TODO: line for SumLocusAncestry_X;
    }
    // set SumLocusAncestry and SumNumArrivals to their averages
    for(int i = 0; i < 2*Populations; ++i) {
      SumLocusAncestryHat[i] /= (double)NumEstepiters; 
      //TODO: X 
    }
    
    cout << "E step\t";
    for(unsigned int g = 0; g < NumIndGametes; ++g) {
      for(int k = 0; k < Populations; ++k) {
	// cout << SumLocusAncestryHat[g*Populations + k] << "\t";
      }
    }
    //cout << "\t" << SumNumArrivals[0] << "\t" << SumNumArrivals[1] << "\t";
    
    SumNumArrivalsHat[0] /= (double)NumEstepiters;   
    SumNumArrivalsHat[1] /= (double)NumEstepiters; 
    SumNumArrivals_XHat[0] /= NumEstepiters; 
    SumNumArrivals_XHat[1] /= NumEstepiters;
    
    //M-step: set mode for log rho and for theta in softmax basis - no minus ones in exponents
    // use vector rhohat and array ThetaProposal to store updates before accept/reject 
    for(unsigned g = 0; g < NumIndGametes; ++g) {
      if( options->isAdmixed(g) ) {
	if(!options->isGlobalRho()){

	  rhohat[g] = (rhoalpha + SumNumArrivalsHat[g]) / (rhobeta + Loci->GetLengthOfGenome());
	  // TODO - fix this for rhoX = 0.5*rho
	  if( Loci->isX_data() && !options->isXOnlyAnalysis() )
	    _rho_X[g] = (rhoalpha + SumNumArrivals_XHat[g]) / (rhobeta + Loci->GetLengthOfXchrm());
	}
	unsigned gg = g; // index of prior (alpha) to use: use alpha[1] if second gamete and no indadmixhiermodel
	if(options->getIndAdmixHierIndicator()) gg = 0;
	double sum = accumulate(alpha[gg].begin(), alpha[gg].end(),0.0, std::plus<double>())
	  + accumulate(SumLocusAncestryHat+g*Populations, 
		       SumLocusAncestryHat+(g+1)*Populations, 0.0, std::plus<double>());
	double sum_X = 0.0;
	if( Loci->isX_data() && !options->isXOnlyAnalysis() )
	  sum_X = accumulate(alpha[gg].begin(), alpha[gg].end(),0.0, std::plus<double>())
	    + accumulate(SumLocusAncestry_XHat+g*Populations, 
			 SumLocusAncestry_XHat+(g+1)*Populations, 0.0, std::plus<double>());
	for(int k = 0; k < Populations; ++k) {
 
	  if(alpha[gg][k]>0.0) ThetaProposal[g*Populations+k] = (alpha[gg][k]+SumLocusAncestryHat[g*Populations+k]) / sum;

	  if( Loci->isX_data() && !options->isXOnlyAnalysis() )
	    ThetaX[g*Populations+k] = (alpha[gg][k]+SumLocusAncestry_XHat[g*Populations+k]) / sum_X;
	}
      } //end isadmixed
    } // end loop over gametes
    double logpriorhat =  LogPriorTheta_Softmax(ThetaProposal, ThetaX, options, alpha) + 
      LogPriorRho_LogBasis(rhohat, _rho_X, options, rhoalpha, rhobeta);
    double loglikhat = getLogLikelihood(options, ThetaProposal, ThetaX, rhohat, _rho_X, false);
    double LogUnnormalizedPosteriorHat  = logpriorhat + loglikhat;
    if(LogUnnormalizedPosteriorHat > LogUnnormalizedPosterior) { //accept update only if density increases 
      setAdmixtureProps(ThetaProposal, NumIndGametes * Populations);
      copy(rhohat.begin(), rhohat.end(), _rho.begin());
      LogUnnormalizedPosterior = LogUnnormalizedPosteriorHat;
    }
    cout << "LogPrior " << logpriorhat << "\tLogLikelihood " << loglikhat << "\tLogUnNormalizedPosterior  " << 
      LogUnnormalizedPosterior << endl << flush;
  } //end EM outer loop


  { //print values to file
    modefile<<setiosflags(ios::fixed)<<setprecision(3);
    modefile << myNumber << "\t";
    if(!options->isGlobalRho()) {
      for(unsigned i = 0; i < NumIndGametes; ++i) modefile<<_rho[i]<<"\t ";
    }
    for(unsigned i = 0; i < NumIndGametes; ++i) { //loop over populations within gametes
      for(int k = 0; k < Populations; ++k) modefile<<Theta[i*Populations +k]<<"\t ";
    }
  }
  delete[] SumLocusAncestryHat;
  
  if(myNumber==1 && options->getChibIndicator()){ // copy modes into hat arrays to use in chib algorithm
    for(unsigned k = 0; k < Populations*NumIndGametes; ++k){
      thetahat[k] = Theta[k];
    }
    //check for zeros where not specified by prior
    for(unsigned i = 0; i < NumIndGametes; ++i){
     unsigned gg = i; // index of alpha to use, 1 only for second gamete and no indadmixhiermodel
      if(options->getIndAdmixHierIndicator()) gg = 0;
      unsigned NumAdmixed = 0;
      double sum = 0.0;
      for(int k = 0; k < Populations; ++k){
	if(alpha[gg][k]>0.0)++NumAdmixed;
	sum += Theta[i*Populations+k];
      }
      int k0;
      for(int k = 0; k < Populations; ++k){
	thetahat[i*Populations+k] /= sum;//correct to sum to 1
	if(alpha[gg][k]>0.0 && thetahat[i*Populations+k]<0.001){
	  k0 = k;
	  for(int kk = 0; kk < Populations; ++kk)if(alpha[gg][kk]>0.0 && kk!=k0)
	    thetahat[i*Populations+kk]-= 0.001 / (double)(NumAdmixed);
	  thetahat[i*Populations+k] = 0.001;
	}
      }
    }//end gamete loop
    copy(_rho.begin(), _rho.end(), rhohat.begin());
    if( Loci->isX_data() ){
      for(unsigned k = 0; k < Populations*NumIndGametes; ++k)thetahatX[k] = ThetaX[k];
      copy(_rho_X.begin(), _rho_X.end(), rhohatX.begin());
    }
  }
}


// ****** End Public Interface *******

void Individual::SampleTheta( int iteration, double *SumLogTheta, 
			      const DataMatrix* const Outcome, 
			      const DataType* const OutcomeType, //const double* const* ExpectedY, 
			      const vector<double> lambda, int NumCovariates,
			      DataMatrix *Covariates, const vector<const double*> beta, const double* const poptheta,
			      const AdmixOptions* const options, const vector<vector<double> > &alpha, //const vector<double> sigma,
			      double DInvLink, double dispersion, bool RW, bool anneal=false)
// samples individual admixture proportions
// called with RW true for a random-walk proposal, false for a conjugate proposal
{
  double logpratio = 0.0;
  if(RW) {
    NumberOfUpdates++;
    logpratio += ProposeThetaWithRandomWalk(options, alpha); 
  } else ProposeTheta(options, /*sigma,*/ alpha, SumLocusAncestry, SumLocusAncestry_X);       

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
  //   //case of X Chromosome and not X only data
  //   if( Loci->isX_data() && !options->isXOnlyAnalysis() )
  //     logpratio += LogAcceptanceRatioForTheta_XChrm( sigma, K);
  
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
	//transform(ThetaX, ThetaX+size_admix, ThetaXHat, ThetaXHat, std::plus<double>());
    }
  }
  if(!IAmUnderTest){
    for( int k = 0; k < K; k++ ){
      SumLogTheta[ k ] += log( Theta[ k ] );
      if(options->isRandomMatingModel() && !options->isXOnlyAnalysis() )
	SumLogTheta[ k ] += log( Theta[ K + k ] );
    }
    
    //increment B using new Admixture Props
    //Xcov is a vector of admixture props as covariates as in UpdateScoreForAncestry
    if(iteration >= options->getBurnIn() && options->getTestForLinkageWithAncestry()){
      double admixtureCovars[Populations-1];
      for(int t = 0; t < Populations-1; ++t)admixtureCovars[t] = Covariates->get(myNumber-1, Covariates->nCols()-Populations+1+t);
      UpdateB(DInvLink, dispersion, admixtureCovars);
    }
  }
}

double Individual::ProposeThetaWithRandomWalk(const AdmixOptions* const options, const vector<vector<double> > &alpha) {
  //TODO: X-chromosome case
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
	if( b[k] ) a[k] = gennor(a[k], step);  
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
  logLikelihood.tempvalue = getLogLikelihood(options, ThetaProposal, ThetaXProposal,_rho, _rho_X, true);
  LogLikelihoodRatio += logLikelihood.tempvalue;
  return LogLikelihoodRatio + LogPriorRatio;// log ratio of full conditionals
}

// Proposes new values for individual admixture proportions 
// as conjugate Dirichlet posterior conditional on prior parameter vector alpha and 
// multinomial likelihood given by sampled values of ancestry at loci where jump indicator xi is 1 (SumLocusAncestry)
// proposes new values for both gametes if random mating model 
void Individual::ProposeTheta(const AdmixOptions* const options, /*const vector<double> sigma, */ const vector<vector<double> > &alpha,
			      int* sumLocusAncestry, int* sumLocusAncestry_X){
  size_t K = Populations;
  if( options->isRandomMatingModel() ){ //random mating model
    for( unsigned int g = 0; g < 2; g++ ) {
      if(options->isAdmixed(g)) {
	unsigned gg = g; // index of which prior (alpha) to use - use alpha[1] if second gamete and no indadmixhiermodel
	if(options->getIndAdmixHierIndicator()) gg = 0;
	for(size_t k = 0; k < K; ++k) { // parameters of conjugate Dirichlet update
	  dirparams[k] = alpha[gg][k] + (double)sumLocusAncestry[k + K*g];
	  if( Loci->isX_data() ) {
	    // if male, paternal elements (0 to K-1) of sumLocusAncestry will be zero
	    dirparams[k] += (double)sumLocusAncestry_X[g*K + k]; 
	  }
	}
	//generate proposal theta from Dirichlet with parameters dirparams
	gendirichlet(K, dirparams, ThetaProposal+g*K );
      } else copy(Theta+g*Populations, Theta+(g+1)*Populations, ThetaProposal+g*Populations);
    } // end loop over gametes

  } else { //assortative mating model
    for(size_t k = 0; k < K; ++k) {
      dirparams[k] = alpha[0][k] + sumLocusAncestry[k] + sumLocusAncestry[k + K];
      if( Loci->isX_data() ) {
      // if male, paternal elements (0 to K-1) of sumLocusAncestry_X will be zero
	dirparams[k] += (double)(sumLocusAncestry_X[k] + sumLocusAncestry_X[K + k]); 
      }
    }
    gendirichlet(K, dirparams, ThetaProposal );
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
      if( myrand() < AccProb ) accept=true;
    } else {
      accept = true;
    }
  }
  
  if(accept) { // set proposed values as new values    
    // if random-mating model and male, thetaX has length K and is set to maternal admixture proportions
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


void Individual::UpdateHMMForwardProbs(unsigned int j, const AdmixOptions* const options, 
				       const double* const theta, const double* const thetaX, 
				       const vector<double> rho, const vector<double> rhoX) {
  //Updates forward probabilities in HMM for chromosome j
  //also sets Diploid flag in Chromosome (last arg of UpdateParameters)
  if( j != X_posn ){// autosome
    Loci->getChromosome(j)->UpdateHMMForwardProbs(theta, GenotypeProbs[j], GenotypesMissing[j], options, rho, true);
  } else {
    if( !SexIsFemale ) { // X chromosome in male individual, haploid
      Loci->getChromosome(j)->UpdateHMMForwardProbs(thetaX, GenotypeProbs[j], GenotypesMissing[j], options, rhoX, false);
    } else { // X chromosome in female individual, diploid
      Loci->getChromosome(j)->UpdateHMMForwardProbs(thetaX, GenotypeProbs[j], GenotypesMissing[j], options, rhoX, true);
    }
  }
  logLikelihood.HMMisOK = false;//because forward probs in HMM have been changed
}

void Individual::SampleRho(const AdmixOptions* const options, double rhoalpha, double rhobeta, 
			   bool updateSumLogRho) {
  vector<unsigned> sumNumArrivals = getSumNumArrivals();
  vector<unsigned> sumNumArrivals_X = getSumNumArrivals_X();
  // rho_X is set to 0.5*rho - equivalent to setting effective length of X chromosome as half length in morgans
  // conjugate gamma update for rho includes arrivals on X chromosome, and 0.5 * length of X chromosome
  double L = Loci->GetLengthOfGenome(); 
  double LX = 0.0;
  if( Loci->isX_data() ) LX = Loci->GetLengthOfXchrm();
  double EffectiveL = 0.0;
  if(options->isRandomMatingModel() ) {
    // SumNumArrivals_X has length 2, and SumNumArrivals_X[0] remains fixed at 0 if male 
    for( unsigned int g = 0; g < 2; g++ )if(options->isAdmixed(g)) {
      do {
	// effective length of genome is L + 0.5*LX if there is an X chrm: i.e. if g=1 or sex is female
	if(g || SexIsFemale) EffectiveL = L + 0.5*LX;
	else EffectiveL = L;
	_rho[g] = gengam( rhoalpha + (double)(sumNumArrivals[g] + sumNumArrivals_X[g]), rhobeta + EffectiveL );
      } while( _rho[g] > TruncationPt || _rho[g] < 1.0 );
      if( Loci->isX_data()  && (g = 1 || SexIsFemale) ) { // no assignment if g = 0 (paternal gamete) and sex is male
	_rho_X[g] = 0.5 * _rho[g];
      }
    } 
  } else {//assortative mating
    // effective length of genome is  2*(L + 0.5*LX) if sex is female, 2*L + 0.5*LX if sex is male
    if( SexIsFemale) EffectiveL =  2.0 * L + LX;
    else EffectiveL = 2.0 * L + 0.5 * LX;
    _rho[0] = gengam( rhoalpha + (double)(sumNumArrivals[0] + sumNumArrivals[1] + sumNumArrivals_X[0] + sumNumArrivals_X[1]), 
			rhobeta + EffectiveL );
    if(Loci->isX_data()) {
      _rho_X[0] = 0.5 * _rho[0];
    }
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
	Outcome->set( myNumber-1, k, gennor( ExpectedY[k][myNumber-1], 1 / sqrt( lambda[k] ) ));
      else{
	if( myrand() * ExpectedY[k][myNumber-1] < 1 )
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
	UpdateHMMForwardProbs(j, options, Theta, ThetaX, _rho, _rho_X);
      }
      //update of score tests for linkage with ancestry requires update of backward probs
      C->UpdateHMMBackwardProbs(Theta, GenotypeProbs[j]);//TODO: pass correct theta for haploid case
      double admixtureCovars[Populations-1];
      for(int t = 0; t < Populations-1; ++t)admixtureCovars[t] = Covariates->get(myNumber-1, Covariates->nCols()-Populations+1+t);
      UpdateScoreTests(options, admixtureCovars, Outcome, OutcomeType, C, DInvLink, dispersion, ExpectedY);
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

void Individual::UpdateScoreForLinkageAffectedsOnly(int locus, int Pops, int k0, bool RandomMatingModel, 
						    const vector<vector<double> > AProbs){
  // values of ancestry risk ratio at which likelihood ratio is evaluated
  double r1 = 0.5;
  double r2 = 2.0;//hard-coding these for now, can make them vary later


  double theta[2];//paternal and maternal admixture proportions

  double Pi[3];//probs of 0,1,2 copies of Pop1 given admixture
  //int offset = 0;
  //if(!RandomMatingModel)offset = Populations;

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
    
    //probs of 0,1,2 copies of Pop1 given admixture
    Pi[2] = theta[0] * theta[1];
    Pi[1] = theta[0] * (1.0 - theta[1]);
    Pi[0] = (1.0 - theta[0]) * (1.0 - theta[1]);
    
    //compute contribution to likelihood ratio
    LikRatio1[locus *Pops + k] += (AProbs[0][k+k0] + sqrt(r1)*AProbs[1][k+k0] + r1 * AProbs[2][k+k0]) / 
      (Pi[0] + sqrt(r1)*Pi[1] + r1*Pi[2]);
    LikRatio2[locus *Pops + k] += (AProbs[0][k+k0] + sqrt(r2)*AProbs[1][k+k0] + r2 * AProbs[2][k+k0]) / 
      (Pi[0] + sqrt(r2)*Pi[1] + r2*Pi[2]);
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
  
  double X[2 * Populations], Xcopy[2*Populations], XX[4*Populations*Populations];
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

//******************** Chib Algorithm ***************************************
// this function does three things:
// 1. sets allele freqs at their prior mean at end of burnin
// 2. calculates log-likelihood and log prior at thetahat, rhohat, allelefreqsMAP at end of burnin
// 3. accumulates posterior ordinate at thetahat, rhohat, allelefreqsMAP at every iteration after burnin 
// should split into three. 
// function should not be called during annealing runs
void Individual::Chib(int iteration, // double *SumLogLikelihood, double *MaxLogLikelihood,
		      const AdmixOptions* const options, const vector<vector<double> > &alpha, 
		      // double globalrho, 
		      double rhoalpha, double rhobeta, double *thetahat,
		      double *thetahatX, vector<double> &rhohat, vector<double> &rhohatX, // LogWriter& Log, 
		      chib *MargLikelihood, AlleleFreqs* A){

  // *** At end of BurnIn ***
  if( iteration == options->getBurnIn() ){
    //set allelefreqsMAP to current values of allelefreqs
    A->setAlleleFreqsMAP();
    //set HapPairProbsMAP to current values of HapPairProbs
    for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
      //if( (*Loci)(j)->GetNumberOfLoci() > 2 )
      (*Loci)(j)->setHaplotypeProbsMAP();
    }
    
    //loglikelihood at estimates
    for(unsigned j = 0; j < Loci->GetNumberOfChromosomes(); ++j)
      SetGenotypeProbs(j, true);//set genotype probs from happairprobsMAP
    // make sure that genotype probs are reset before next update
    MargLikelihood->setLogLikelihood(getLogLikelihood( options, thetahat, thetahatX, rhohat, rhohatX, true));
  }
  
  // *** After BurnIn ***
  if( iteration > options->getBurnIn() ){
    double LogPrior  = LogPriorTheta(thetahat, thetahatX, options, alpha)
      + LogPriorRho(rhohat, rhohatX, options, rhoalpha, rhobeta);

    double LogPosterior = 0.0;
    double LP = 0.0;
    if( Populations > 1 ){
      LP = CalculateLogPosteriorTheta(options, thetahat, thetahatX, alpha);
      logPosterior[0].push_back(LP);
      LogPosterior += LP;
      LP = CalculateLogPosteriorRho(options, rhohat, rhohatX, rhoalpha, rhobeta);
      logPosterior[1].push_back(LP);
      LogPosterior += LP;
    }
    if( A->IsRandom() ){
      double LogPosteriorFreqs = 0.0, LogPriorFreqs = 0.0;
      for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
	for( int k = 0; k < Populations; k++ ){
	  vector<double> args = A->GetPriorAlleleFreqs(j,k);
	  vector<int> counts = A->GetAlleleCounts(j,k);
	  transform(counts.begin(), counts.end(), args.begin(), args.begin(), plus<double>());//PriorAlleleFreqs + AlleleCounts
	  LogPosteriorFreqs += getDirichletLogDensity( args, A->getAlleleFreqsMAP(j, k));//LogPosterior for Allele Freqs
	  LogPriorFreqs += getDirichletLogDensity( A->GetPriorAlleleFreqs(j, k), A->getAlleleFreqsMAP(j,k) );
	}
      }
      LogPosterior += LogPosteriorFreqs;
      logPosterior[2].push_back( LogPosteriorFreqs  );
      LogPrior += LogPriorFreqs;
    }
    //accumulates vector of sampled values of log posterior density at estimates
    MargLikelihood->addLogPrior(LogPrior);
    MargLikelihood->addLogPosteriorObs( LogPosterior );
    // *SumLogLikelihood += logLikelihood; // no point in accumulating this
  }
}

double Individual::LogPriorTheta(const double* const theta, const double* const thetaX,  
			    const AdmixOptions* const options, const vector<vector<double> > &alpha) const {
  //Computes log prior density at supplied parameter values
  // calls getDirichletLogDensity with parameters as std vector, proportions as double 
  int K = Populations;
  double LogPrior=0.0;

   // ** case of xonly data **
   if( options->isXOnlyAnalysis() ){
      LogPrior += getDirichletLogDensity( alpha[0], theta );
   }
   // ** case of some x data **
   else if( Loci->isX_data() ){
      LogPrior += getDirichletLogDensity( alpha[0], theta )//first gamete
	+ getDirichletLogDensity( alpha[0], thetaX )// X chromosome
	+ getDirichletLogDensity( alpha[1], theta+K )//second gamete
	+ getDirichletLogDensity( alpha[1], thetaX + K );//Xchromosome
   }
   else{
      if( Populations > 1 ){
	if( options->isAdmixed(0) ){//gamete 1
            LogPrior += getDirichletLogDensity( alpha[0], theta );
         }
	if( options->isAdmixed(1) ){//gamete 2
            LogPrior += getDirichletLogDensity( alpha[1], theta + K );
         }
      }
   }

   return LogPrior;
}

double Individual::LogPriorTheta_Softmax(const double* const theta, const double* const thetaX,  
			    const AdmixOptions* const options, const vector<vector<double> > &alpha) const {
  // Computes LogPrior density in softmax basis at supplied parameter values
  // calls getDirichletLogDensity_Softmax with parameters as std vector, proportions as array 
  int K = Populations;
  double LogPrior=0.0;
  // ** case of xonly data **
  if( options->isXOnlyAnalysis() ){
    LogPrior += getDirichletLogDensity_Softmax( alpha[0], theta );
  }
  // ** case of some x data **
  else if( Loci->isX_data() ){
    LogPrior += getDirichletLogDensity_Softmax( alpha[0], theta )//first gamete
      + getDirichletLogDensity_Softmax( alpha[0], thetaX )// X chromosome
      + getDirichletLogDensity_Softmax( alpha[1], theta+K )//second gamete
      + getDirichletLogDensity_Softmax( alpha[1], thetaX + K );//Xchromosome
  } else {
    if( Populations > 1 ){
      if( options->isAdmixed(0) ){//gamete 1
	LogPrior += getDirichletLogDensity_Softmax( alpha[0], theta );
      }
      if( options->isAdmixed(1) ){//gamete 2
	LogPrior += getDirichletLogDensity_Softmax( alpha[1], theta + K );
      }
    }
  }
  return LogPrior;
}

//called on test individual after burnin
// should replace with a function to calculate posterior/prior ratio
double Individual::CalculateLogPosteriorTheta(const AdmixOptions* const options, const double* const theta, const double* const thetaX, 
					      const vector<vector<double> > &alpha) const{
  // calculates log full conditional at theta, conditional on realized locus ancestry states and jump indicators
  double LogPosterior = 0.0;
  
  vector<double> alphaparams1(Populations), alphaparams0(Populations);
  if( options->isXOnlyAnalysis() ){
    transform(alpha[0].begin(), alpha[0].end(), SumLocusAncestry_X, alphaparams0.begin(), std::plus<double>());
      LogPosterior += getDirichletLogDensity(alphaparams0, theta);
  }
  else if( Loci->isX_data() ){
    for( unsigned int g = 0; g < 2; g++ ){
      transform(alpha[g].begin(), alpha[g].end(), SumLocusAncestry +g*Populations, alphaparams0.begin(), std::plus<double>());
      LogPosterior += getDirichletLogDensity(alphaparams0, theta + g*Populations);
    }
    for( unsigned int g = 0; g < gametes[X_posn]; g++ ){
      transform(alpha[g].begin(), alpha[g].end(), SumLocusAncestry_X +g*Populations, alphaparams0.begin(), std::plus<double>());
      LogPosterior += getDirichletLogDensity(alphaparams0, thetaX + g*Populations);
    }
  }
  else if( options->isSymmetric() ){//both gametes have same admixture
    vector<double> x(2,0.0);
    transform(alpha[0].begin(), alpha[0].end(), SumLocusAncestry, alphaparams0.begin(), std::plus<double>());
    x[0] += getDirichletLogDensity(alphaparams0, theta);
    x[1] += getDirichletLogDensity(alphaparams0, theta+Populations);
    transform(alpha[1].begin(), alpha[1].end(), SumLocusAncestry+Populations, alphaparams1.begin(), std::plus<double>());
      x[0] += getDirichletLogDensity(alphaparams1, theta+Populations);
      x[1] += getDirichletLogDensity(alphaparams1, theta);
      if( isnan(x[0]) || isinf(x[0]) )
	LogPosterior = x[1] - log(2.0);
      else if( isnan(x[1]) || isinf(x[1])  )
	LogPosterior = x[0] - log(2.0);
      else if( x[0] < x[1] )
	LogPosterior = x[1] + log( 1 + exp( x[0] - x[1] ) ) - log(2.0);
      else
	LogPosterior = x[0] + log( 1 + exp( x[1] - x[0] ) ) - log(2.0);
      //if( isnan(LogPosterior) || isinf(LogPosterior) ){
      //PR(alphaparams0);
      //PR(alphaparams1);
      //PR(AdmixtureHat.GetColumn(0));
      //   PR(AdmixtureHat.GetColumn(1));
      //   PR(x[0]);
      //   PR(x[1]);
      //   exit(0);
      //}
  }
  else{//different admixtures for each gamete
    if(  options->isAdmixed(0) ){//admixed first gamete
      //alphaparams0 = alpha + SumLocusAncestry
      transform(alpha[0].begin(), alpha[0].end(), SumLocusAncestry, alphaparams0.begin(), std::plus<double>());
      LogPosterior += getDirichletLogDensity(alphaparams0, theta);
    }
    if(  options->isAdmixed(1) ){//admixed second gamete
      transform(alpha[1].begin(), alpha[1].end(), SumLocusAncestry+Populations, alphaparams1.begin(), std::plus<double>());
      LogPosterior += getDirichletLogDensity(alphaparams1, theta+Populations);
    }
  }
  return LogPosterior;
}

double Individual::LogPriorRho(const vector<double> rho, const vector<double> rhoX, 
			    const AdmixOptions* const options, double rhoalpha, double rhobeta) const {
//Computes LogPrior at supplied parameter values
// standardizes densities with an integrating constant - can simplify by just calculating posterior/prior ratio 
  double LogPrior=0.0;

   // ** case of xonly data **
   if( options->isXOnlyAnalysis() ){
         LogPrior = getGammaLogDensity( rhoalpha, rhobeta, rho[0] );
         LogPrior -= log( gsl_cdf_gamma_Q(rhobeta, rhoalpha, 1.0) );
   }
   // ** case of some x data **
   else if( Loci->isX_data() ){
     LogPrior = getGammaLogDensity( rhoalpha, rhobeta, rho[0] )//gamete1
       + getGammaLogDensity( rhoalpha, rhobeta, rhoX[0] )      //X chromosome
       + getGammaLogDensity( rhoalpha, rhobeta, rho[1] )       //gamete2
       + getGammaLogDensity( rhoalpha, rhobeta, rhoX[1] );    //X chr
     LogPrior /= gsl_cdf_gamma_Q(rhobeta, rhoalpha, 1.0);
   }
   else{
     if( Populations > 1 ){
       for(int g = 0; g < 2; ++g){//loop over gametes
	 if( options->isAdmixed(g) ){
	   
	   LogPrior = getGammaLogDensity( rhoalpha, rhobeta, rho[g] );
	   LogPrior -= log( gsl_cdf_gamma_Q(rhobeta, rhoalpha, 1.0) );
	 }
       }
      }
   }
   return LogPrior;
}

double Individual::LogPriorRho_LogBasis(const vector<double> rho, const vector<double> rhoX, 
			    const AdmixOptions* const options, double rhoalpha, double rhobeta) const {
// computes LogPrior density in log rho basis at supplied parameter values
// no standardization by integrating constant
  double LogPrior=0.0;
   // ** case of xonly data **
   if( options->isXOnlyAnalysis() ){
         LogPrior = getGammaLogDensity_LogBasis( rhoalpha, rhobeta, rho[0] );
   }
   // ** case of some x data **
   else if( Loci->isX_data() ){
     LogPrior = getGammaLogDensity_LogBasis( rhoalpha, rhobeta, rho[0] )//gamete1
       + getGammaLogDensity_LogBasis( rhoalpha, rhobeta, rhoX[0] )      //X chromosome
       + getGammaLogDensity_LogBasis( rhoalpha, rhobeta, rho[1] )       //gamete2
       + getGammaLogDensity_LogBasis( rhoalpha, rhobeta, rhoX[1] );    //X chr
   }
   else{
     if( Populations > 1 ){
       for(int g = 0; g < 2; ++g){//loop over gametes
	 if( options->isAdmixed(g) ){
	   LogPrior = getGammaLogDensity_LogBasis( rhoalpha, rhobeta, rho[g] );
	 }
       }
      }
   }
   return LogPrior;
}

double Individual::CalculateLogPosteriorRho(const AdmixOptions* const options,  
					    const vector<double> rho, const vector<double> rhoX,
					    double rhoalpha, double rhobeta)const{
  // calculates log full conditional density at sum-intensities rho, conditional on realized number of arrivals
  double LogPosterior = 0.0;
  double L = Loci->GetLengthOfGenome(), L_X = 0.0;
  if( Loci->isX_data() ) L_X = Loci->GetLengthOfXchrm();
  vector<unsigned> SumN = getSumNumArrivals();
  vector<unsigned> SumN_X = getSumNumArrivals_X();

  double IntConst1;
  vector<double> alphaparams1(Populations), alphaparams0(Populations);
  if( options->isXOnlyAnalysis() ){
    LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN_X[0], rhobeta + L_X, rho[0] );
    //if(!options->RhoFlatPrior() && !options->logRhoFlatPrior() )
      IntConst1 = IntegratingConst(rhoalpha+(double)SumN_X[0], rhobeta+L_X, 1.0, options->getTruncPt() );
      //else
      //IntConst1 = gsl_cdf_gamma_Q(rhobeta+L_X, rhoalpha+(double)SumNumArrivals_X[0], 1.0);
    LogPosterior -= log(IntConst1);
  }
  else if( Loci->isX_data() ){
    for( unsigned int g = 0; g < 2; g++ ){
      LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN[g], rhobeta + L, rho[g] );
      //if(!options->RhoFlatPrior() && !options->logRhoFlatPrior() )
	IntConst1 = IntegratingConst(rhoalpha+(double)SumN[g], rhobeta+L, 1.0, options->getTruncPt() );
	//else
	//IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumNumArrivals[g], 1.0);
      LogPosterior -= log(IntConst1);
    }
    for( unsigned int g = 0; g < gametes[X_posn]; g++ ){
      LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN_X[g], rhobeta + L_X, rhoX[g] );
      //if( options->RhoFlatPrior() || options->logRhoFlatPrior() )
      //IntConst1 = IntegratingConst(rhoalpha+(double)SumNumArrivals_X[g], rhobeta+L_X, 1.0, options->getTruncPt() );
      //else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L_X, rhoalpha+(double)SumN_X[g], 1.0);
      LogPosterior -= log(IntConst1);
    }
  }
  else if( options->isSymmetric() ){//both gametes have same admixture
    vector<double> x(2,0.0);
    x[0] += getGammaLogDensity( rhoalpha + (double)SumN[0], rhobeta + L, rho[0] );
    x[1] += getGammaLogDensity( rhoalpha + (double)SumN[1], rhobeta + L, rho[0] );
    x[0] += getGammaLogDensity( rhoalpha + (double)SumN[1], rhobeta + L, rho[1] );
    x[1] += getGammaLogDensity( rhoalpha + (double)SumN[0], rhobeta + L, rho[1] );
    double IntConst1, IntConst2;
    //if( options->RhoFlatPrior() || options->logRhoFlatPrior() ){
    //IntConst1 = IntegratingConst(rhoalpha+(double)SumNumArrivals[0], rhobeta+L, 1.0, options->getTruncPt() );
    //IntConst2 = IntegratingConst(rhoalpha+(double)SumNumArrivals[1], rhobeta+L, 1.0, options->getTruncPt() );
    //}
    //else{
	IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[0], 1.0);
	IntConst2 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[1], 1.0);
	//}
    x[0] -= log(IntConst1);
    x[1] -= log(IntConst2);
      if( isnan(x[0]) || isinf(x[0]) )
	LogPosterior = x[1] - log(2.0);
      else if( isnan(x[1]) || isinf(x[1])  )
	LogPosterior = x[0] - log(2.0);
      else if( x[0] < x[1] )
	LogPosterior = x[1] + log( 1 + exp( x[0] - x[1] ) ) - log(2.0);
      else
	LogPosterior = x[0] + log( 1 + exp( x[1] - x[0] ) ) - log(2.0);
  }
  else{//different admixtures for each gamete
    if(  options->isAdmixed(0) ){//admixed first gamete
      LogPosterior = getGammaLogDensity( rhoalpha + (double)SumN[0], rhobeta + L, rho[0] );
      //if( options->RhoFlatPrior() || options->logRhoFlatPrior() )
	IntConst1 = IntegratingConst(rhoalpha+(double)SumN[0], rhobeta+L, 1.0, options->getTruncPt() );
	//else
	//IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumNumArrivals[0], 1.0);
      LogPosterior -= log( IntConst1 );
    }
    if(  options->isAdmixed(1) ){//admixed second gamete
      LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN[1], rhobeta + L, rho[1] );
      //if( options->RhoFlatPrior() || options->logRhoFlatPrior() )
      //IntConst1 = IntegratingConst(rhoalpha+(double)SumNumArrivals[1], rhobeta+L, 1.0, options->getTruncPt() );
      //else
	IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[1], 1.0);
      LogPosterior -= log( IntConst1 );
    }
  }
  return LogPosterior;
}

double Individual::IntegratingConst( double alpha, double beta, double a, double b )const
{
   double I = gsl_cdf_gamma_P( b*beta, alpha, 1 ) - gsl_cdf_gamma_P( a*beta, alpha, 1);
   return I;
}

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
