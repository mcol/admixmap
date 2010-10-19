/**
 *   ADMIXMAP
 *   AdmixedIndividual.cc
 *   Class to represent an individual in an admixture model and update individual-level parameters
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *
 * This program is free software distributed WITHOUT ANY WARRANTY.
 * You can redistribute it and/or modify it under the terms of the GNU General Public License,
 * version 2 or later, as published by the Free Software Foundation.
 * See the file COPYING for details.
 *
 */
#include "AdmixedIndividual.h"
#include "config.h" // USE_GENOTYPE_PARSER
#include "AdmixOptions.h"
#include "InputAdmixData.h"
#include "GenotypeIterator.h"
#include "bclib/Regression.h"
#include "bclib/misc.h"
#include "bclib/dist.h"
//#include "bclib/linalg.h"
#include <algorithm>
#include <limits>
#include <sstream>


#define PR(x) cout << #x << " = " << x << endl;

using bclib::Rand;
using genepi::cvector;
using genepi::RhoType;

using std::vector;


// Print the first drawing of the admixture proportions
#define DEBUG_INITIAL_ADMIXTURE_PROPS 0


#if 0
    #define DEBUG_TH_PROP(X) X
#else
    #define DEBUG_TH_PROP(X)
#endif



//******** Constructors **********
// AdmixedIndividual::AdmixedIndividual() {//should initialise pointers here
//   dirparams = 0;
//   SumSoftmaxTheta = 0;
//   SumLocusAncestry = 0;
//   SumLocusAncestry_X = 0;
// }

AdmixedIndividual::AdmixedIndividual(int number, const AdmixOptions* const options, const InputAdmixData* const Data,
				     bool undertest=false):Individual(number),IAmUnderTest(undertest){
  GenotypesMissing = new bool*[numChromosomes];

  for( unsigned int j = 0; j < numChromosomes; j++ ){
    GenotypesMissing[j] = new bool[ Loci->GetSizeOfChromosome(j) ];
  }

  //retrieve genotypes
  Data->GetGenotype(number, *Loci, &genotypes, GenotypesMissing);

  #if USE_GENOTYPE_PARSER
    isHaploid = genotypes.at(0).at(0).isHaploid(); //note: assumes at least one autosome before X-chr
  #else
    isHaploid = (genotypes[0][0].size() == 1);	   //note: assumes at least one autosome before X-chr
  #endif

  Individual::Initialise(options, Data);
  int numCompositeLoci = Loci->GetNumberOfCompositeLoci();

  Theta.setDimensions(NumGametes, NumHiddenStates);
  SetUniformAdmixtureProps();

  // loop over composite loci to set possible haplotype pairs compatible with genotype
  for(unsigned j = 0; j < (unsigned)numCompositeLoci; ++j) {
    #if USE_GENOTYPE_PARSER
	const ploidy p = genotypes[j][0].isDiploid() ? diploid : haploid;
    #else
	ploidy p = (genotypes[j][0].size()>1) ? diploid : haploid;
    #endif
    AdmixGenotypeIterator G(genotypes[j], p);
    (*Loci)(j)->HaplotypeSetter.setPossibleHaplotypePairs(&G, PossibleHapPairs[j]);

    // initialise sampledHapPairs with the first of the possible happairs.
    // if only one possible happair or if annealing (which uses hamiltonian sampler), sampling of hap pair will be skipped.
    sampledHapPairs.push_back(PossibleHapPairs[j][0]);
  }
  //Now the PossibleHapPairs have ben determined and missing genotype indicators have been set,
  //the genotypes are deleted as they are no longer needed
  if( options->getHWTestIndicator())SetMissingGenotypes();
  DeleteGenotypes();

  //allocate genotype probs - these are actually the emission probs given each phased hidden state
  GenotypeProbs = new double*[numChromosomes];
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    if(isHaploid || (!SexIsFemale && Loci->isXChromosome(j))){//haploid on this chromosome
      GenotypeProbs[j] = new double[Loci->GetSizeOfChromosome(j)*NumHiddenStates];
    }
    else{
      GenotypeProbs[j] = new double[Loci->GetSizeOfChromosome(j)*NumHiddenStates*NumHiddenStates];
    }
  }

  AncestryProbs = false;
  if (options->getLocusAncestryProbsIndicator()) {
    AncestryProbs = true;
    // vector for the sum over iterations of hidden state copy number probs
    // (of dimension Loci x Populations x 3)
    SumProbs.resize(numCompositeLoci * NumHiddenStates * 3, 0.0);
  }

  if(options->getChibIndicator() || options->getIndAdmixModeFilename())
    thetahat.setDimensions(NumGametes, NumHiddenStates);

  InitialiseSumIntensities(options);
  SumSoftmaxTheta = new double[ NumHiddenStates * NumGametes ];
  fill(SumSoftmaxTheta, SumSoftmaxTheta + NumHiddenStates*NumGametes, 0.0);
  dirparams = new double[NumHiddenStates]; //to hold dirichlet parameters for conjugate updates of theta
  ThetaProposal.setDimensions(NumGametes, NumHiddenStates);
  // SumLocusAncestry is sum of locus ancestry states over loci at which jump indicator xi is 1
  SumLocusAncestry = new int[NumHiddenStates * 2];
  SumLocusAncestry_X = new int[NumHiddenStates * 2];
  if(!options->isGlobalRho())SumNumArrivals.resize(2*numCompositeLoci);

  // ** set up StepSizeTuner object for random walk updates of admixture **
  NumberOfUpdates = 0;
  w = 1;
  step0 = 0.3; // initial sd of random walk proposal distribution
  step = step0;
  ThetaTuner.SetParameters( step0, 0.0001, 10.0, 0.44);

  if (Xdata) {

    psi.resize(NumHiddenStates, 1.0);

    if (!options->isGlobalPsi()) {
      psistep.resize(NumHiddenStates, 1.0);
      SumLogPsi.resize(NumHiddenStates);
      TunePsiSampler.resize(NumHiddenStates);
      NumberOfPsiUpdates = 0;

      // initialise the tuner for psi: we create a stepsize tuner for each
      // element of psi, but we never use the first element as it corresponds
      // to psi[0] which is 1.0 by definition
      for (int i = 1; i < NumHiddenStates; ++i)
        TunePsiSampler[i].SetParameters(1.0, 0.01, 10, 0.44);
    }
  }
}

void AdmixedIndividual::InitialiseSumIntensities(const AdmixOptions* const options){
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
  _rho.assign(NumGametes, 0.0); // set to 0 for unadmixed gametes
  rhohat.assign(NumGametes, 0.0); // set to 0 for unadmixed gametes

  for(unsigned g = 0; g < NumGametes; ++g) {
    if(options->isAdmixed(g)) {
      _rho[g] = init;
      rhohat[g] = init;
    }
  }
  sumlogrho.assign(_rho.size(), 0.0);
}
///sets possible hap pairs for a single SNP
void AdmixedIndividual::SetPossibleHaplotypePairs(const vector<vector<unsigned short> > Genotype, vector<hapPair> &PossibleHapPairs){
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

void AdmixedIndividual::SetMissingGenotypes(){
  //allocates and sets an array of bools indicating whether genotypes at each locus are missing
  //used in HW score test; NB call before genotypes are deleted
  if(genotypes.size()==0)throw string("determining missing genotypes after genotypes have been deleted");
  missingGenotypes = new bool[Loci->GetTotalNumberOfLoci()];
  unsigned index = 0;
  for(unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); ++j)
    for(int k = 0; k < Loci->getNumberOfLoci(j); ++k){
      #if USE_GENOTYPE_PARSER
	missingGenotypes[index++] = genotypes[j][k].isMissing2();
      #else
	missingGenotypes[index++] = (genotypes[j][k][0] == 0);
      #endif
    }
}
//********** Destructor **********
AdmixedIndividual::~AdmixedIndividual() {
  delete[] SumSoftmaxTheta;
  delete[] dirparams;
  delete[] SumLocusAncestry;
  delete[] SumLocusAncestry_X;
  //this might not work, relies on Loci still being in scope in top level
  //GPArray.dealloc(Loci->GetNumberOfCompositeLoci());
}

void AdmixedIndividual::SetStaticMembers( Genome & pLoci, const Options & options ){
  Individual::SetStaticMembers( pLoci, options );
}

// draw initial values for admixture proportions theta from Dirichlet prior
void AdmixedIndividual::drawInitialAdmixtureProps( const AlphaType & alpha ) {
  const size_t K = NumHiddenStates;
  for( unsigned g = 0; g < NumGametes; ++g ) {
    double sum = 0.0;
    for(size_t k = 0; k < K; ++k) {
      sum += alpha[g][k];
    }
    for(size_t k = 0; k < K; ++k) {
      thetahat[g][k] = alpha[g][k] / sum; // set thetahat to prior mean
      dirparams[k] = alpha[g][k];
      //cout << dirparams[k] << " ";
    }
    // draw theta from Dirichlet with parameters dirparams
    Rand::gendirichlet( K, dirparams, Theta[g] );

#if DEBUG_INITIAL_ADMIXTURE_PROPS
    cout << "drawInitialAdmixtureProps " << getMyNumber() << ": ";
    for (size_t k = 0; k < K; ++k)
      cout << Theta[g][k] << " ";
    cout << endl;
#endif
  }
}


//********** set admixture proportions *********
void AdmixedIndividual::setAdmixtureProps(const AdmixtureProportions& rhs) {

  Theta = rhs;
}


// chibindicator is passed to the CompositeLocus object: if set to true,
// CompositeLocus will use HapPairProbsMAP instead of HapPairProbs when the
// allele frequencies are not fixed.
void AdmixedIndividual::SetGenotypeProbs(int j, int jj,
                                         unsigned locus, bool chibindicator) {

  const bool diploid = !isHaploid && (j != (int)X_posn || SexIsFemale);
  const int K = diploid ? NumHiddenStates*NumHiddenStates : NumHiddenStates;

  if (!GenotypesMissing[j][jj]) {
    if (diploid)
      (*Loci)(locus)->GetGenotypeProbs(GenotypeProbs[j] + jj*K,
                                       PossibleHapPairs[locus],
                                       chibindicator);
    else
      (*Loci)(locus)->GetHaploidGenotypeProbs(GenotypeProbs[j] + jj*K,
                                              PossibleHapPairs[locus],
                                              chibindicator);
  }
  else {
    for (int k = 0; k < K; ++k)
      GenotypeProbs[j][jj * K + k] = 1.0;
  }
}


/// Called after energy has been evaluated, before updating model parameters
void AdmixedIndividual::AnnealGenotypeProbs(int j, const double coolness) {

  int locus = Loci->getChromosome(j)->GetLocus(0);

  const bool diploid = !isHaploid && (j != (int)X_posn || SexIsFemale);
  const int K = diploid ? NumHiddenStates*NumHiddenStates : NumHiddenStates;

  // loop over composite loci
  for (unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); ++jj) {

    if (!GenotypesMissing[j][jj]) {
      for (int k = 0; k < K; ++k) // loop over ancestry states
        GenotypeProbs[j][jj*K + k] = pow(GenotypeProbs[j][jj*K + k], coolness);
    }
    locus++;
  }
}

//******************** Accessors ***********************************************************

///returns sum of sumintensities over gametes
double AdmixedIndividual::getSumrho()const {
  double sumrho = 0;
  for( unsigned int g = 0; g < _rho.size(); g++ )
    sumrho += _rho[g];
  return sumrho;
}

const RhoType & AdmixedIndividual::getRho() const {
   return _rho;
}

double AdmixedIndividual::getPsi(int pop) const {
  return psi[pop];
}

const int *AdmixedIndividual::getSumLocusAncestry()const{
  return SumLocusAncestry;
}
const int *AdmixedIndividual::getSumLocusAncestryX()const{
  return SumLocusAncestry_X;
}

/// Return the number of arrivals across the genome
const vector<unsigned> AdmixedIndividual::getSumNumArrivals() const {
  vector<unsigned> SumN(2,0);
  unsigned locus = 0;
  for(unsigned i = 0; i < numChromosomes; ++i){
      for(unsigned j = 0; j < Loci->GetSizeOfChromosome(i); ++j, ++locus)
	if(i != X_posn){//skip X chromosome
	  //accumulate sums over loci
	  SumN[0] += SumNumArrivals[locus*2];//first gamete
	  if(!isHaploidatLocus(locus)) SumN[1] += SumNumArrivals[locus*2+1];//second gamete
	}
  }
  return SumN;
}

const vector<unsigned>AdmixedIndividual::getSumNumArrivals_X()const{
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
void AdmixedIndividual::getSumNumArrivals(vector<unsigned> *sum)const{
  //should check size of argument
  for(unsigned i = 0; i < sum->size(); ++i){//sum over gametes
    (*sum)[i] += SumNumArrivals[2*i] + SumNumArrivals[2*i + 1];
  }
}

//****************** Log-Likelihoods **********************
// public function:
// calls private function to get log-likelihood at current parameter values, and stores it either as loglikelihood.value or as loglikelihood.tempvalue
// store should be false when calculating energy for an annealed run, or when evaluating proposal for global sum-intensities
double AdmixedIndividual::getLogLikelihood( const Options& options, const bool forceUpdate, const bool store) {

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
double AdmixedIndividual::getLogLikelihood(const Options& options,
                                           const AdmixtureProportions& theta,
                                           const RhoType& rho, bool updateHMM) {
  double LogLikelihood = 0.0;
  if(NumHiddenStates == 1) LogLikelihood = getLogLikelihoodOnePop();
  else {
        LogLikelihood = Individual::getLogLikelihood(options, theta, rho, updateHMM);
  }
  return LogLikelihood; // if HMM update not required, can just use stored log-likelihood
}


double AdmixedIndividual::getLogLikelihoodAtPosteriorMeans(const Options& options) {
  // should set allele freqs also to posterior means, and recalculate prob genotypes at these freqs before calling getloglikelihood

  double LogLikelihood = 0.0;
  if(NumHiddenStates == 1)
    LogLikelihood = getLogLikelihoodOnePop();
  else {
    AdmixtureProportions ThetaBar(NumGametes, NumHiddenStates);
    cvector<double> rhobar;

    //getPosteriorMeans(ThetaBar, rhobar, options->getTotalSamples() - options->getBurnIn());
    unsigned size = NumHiddenStates * NumGametes;
    const double scale_factor = options.getTotalSamples() - options.getBurnIn();
    for (unsigned i = 0; i < size; ++i)
      SumSoftmaxTheta[i] /= scale_factor;

    for (unsigned i = 0; i < sumlogrho.size(); ++i)
      rhobar.push_back( exp(sumlogrho[i]/scale_factor) );

    //apply softmax transformation to obtain thetabar
    double *a = new double[NumHiddenStates];
    bool* b = new bool[NumHiddenStates];
    for( unsigned int g = 0; g < NumGametes; g++ ){
      for(int k = 0; k < NumHiddenStates; ++k)
        b[k] = (Theta[g][k] > 0.0) ? true : false; // skip elements set to zero
      bclib::softmax(NumHiddenStates, a, SumSoftmaxTheta+g*NumHiddenStates, b);
      for (int k = 0; k < NumHiddenStates; ++k) {
        if (b[k]) ThetaBar[g][k] = a[k];
        // rescale SumSoftmaxTheta back
        SumSoftmaxTheta[g*NumHiddenStates + k] *= scale_factor;
      }
    }
    delete[] a;
    delete[] b;

    for( unsigned int j = 0; j < numChromosomes; j++ ) {
      UpdateHMMInputs(j, options, ThetaBar, rhobar);
      LogLikelihood += Loci->getChromosome(j)->HMM->getLogLikelihood( !isHaploid && (!Loci->isXChromosome(j) || SexIsFemale) );
    }
  }

  return LogLikelihood;
}

double AdmixedIndividual::getLogLikelihoodOnePop(){ //convenient for a single population as no arguments required
  double LogLikelihood = 0.0;
  double Prob; //one pop so 1x1 array
  for( unsigned int j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ) {// loop over composite loci
    if(!GenotypeIsMissing(j)){
      (*Loci)(j)->GetGenotypeProbs(&Prob,getPossibleHapPairs(j), false);
      LogLikelihood += log( Prob );
    }
  }
  return LogLikelihood;
}

void AdmixedIndividual::getPosteriorMeans(double* ThetaMean, RhoType& rhoMean,
                                          unsigned samples) const {
  unsigned size = NumHiddenStates * NumGametes;
  const double dSamples = samples;
  for (unsigned i = 0; i < size; ++i)
    SumSoftmaxTheta[i] /= dSamples;
  rhoMean.clear();
  for(unsigned i = 0; i < sumlogrho.size(); ++i)
    rhoMean.push_back( exp(sumlogrho[i]/dSamples) );

  if (ThetaMean == 0)
    throw string("ERROR in AdmixedIndividual::getPosteriorMeans: "
                 "ThetaMean is NULL");

  // apply softmax transformation to obtain thetabar
  bool* b = new bool[NumHiddenStates];
  for( unsigned int g = 0; g < NumGametes; g++ ){
    for(int k = 0; k < NumHiddenStates; ++k)
      b[k] = (Theta[g][k] > 0.0) ? true : false; // skip elements set to zero
    bclib::softmax(NumHiddenStates, ThetaMean+g*NumHiddenStates, SumSoftmaxTheta+g*NumHiddenStates, b);
    // rescale SumSoftmaxTheta back
    for(int k = 0; k < NumHiddenStates; ++k)
      SumSoftmaxTheta[g*NumHiddenStates + k] *= dSamples;
  }
  delete[] b;
}

void AdmixedIndividual::WritePosteriorMeans(ostream& os, unsigned samples, bool globalrho)const{
  double* ThetaBar = new double[NumGametes*NumHiddenStates];
  cvector<double> rhobar;

  getPosteriorMeans(ThetaBar, rhobar, samples);

  copy(ThetaBar, ThetaBar + NumHiddenStates * NumGametes, ostream_iterator<double>(os, "\t"));
  if(!globalrho)
    copy(rhobar.begin(), rhobar.end(), ostream_iterator<double>(os, "\t"));

  delete[] ThetaBar;

  if ( AncestryProbs )
    SumProbs /= (double) samples;
}

void AdmixedIndividual::WritePosteriorMeansXChr(ostream& os, unsigned samples) const {

  // We do not explicitly store the X chromosome admixtures, but they are
  // computed from the autosomal admixtures. Therefore, to get the posterior
  // means for the X chromosome admixtures, we first have to get the autosomal
  // posterior means, store them in an AdmixtureProportions object and then
  // ask it to compute the X chromosome admixtures.

  const unsigned size = NumGametes * NumHiddenStates;
  AdmixtureProportions ThetaPosteriorMeans(NumGametes, NumHiddenStates);
  double *ThetaMean = new double[size];
  cvector<double> rhobar;

  // get the autosomal posterior means and store them
  getPosteriorMeans(ThetaMean, rhobar, samples);
  for (unsigned int g = 0; g < NumGametes; ++g)
    for (int k = 0; k < NumHiddenStates; ++k)
      ThetaPosteriorMeans[g][k] = ThetaMean[NumHiddenStates * g + k];

  // get the X chromosome posterior means
  const double *ThetaMeanX = ThetaPosteriorMeans.flatXChromosome(psi);
  copy(ThetaMeanX, ThetaMeanX + size, ostream_iterator<double>(os, "\t"));

  delete[] ThetaMean;
}

void AdmixedIndividual::WritePosteriorMeansLoci(ostream& os)const{
  //getPosteriorMeans(ThetaBar, rhobar, samples);
  const unsigned int size = SumProbs.size();
  for (unsigned int i = 0; i < size - 1; ++i)
    os << SumProbs[i] << ",";
  os << SumProbs[size - 1]; // last element without following comma
}

//************** Updating (Public) **********************************************************
void AdmixedIndividual::ResetSufficientStats(){
  if(NumHiddenStates>1) {
    // ** reset SumLocusAncestry to zero
    for(int j = 0; j < NumHiddenStates *2; ++j) {
      SumLocusAncestry[j] = 0;
      SumLocusAncestry_X[j] = 0;
    }

    //SumNumArrivals is the number of arrivals between each pair of adjacent loci
    fill(SumNumArrivals.begin(), SumNumArrivals.end(), 0);
  }
}


void AdmixedIndividual::SampleJumpIndicators(bool sampleArrivals){
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

/// Use an EM algorithm to search for posterior modes of individual parameters
/// theta and rho
// uses current values of allele freqs
void AdmixedIndividual::FindPosteriorModes(const AdmixOptions& options,
                                           const AlphaType& alpha,
                                           double rhoalpha, double rhobeta,
                                           AlleleFreqs* A, ofstream& modefile) {
  if( A->IsRandom() ) {
    // set genotype probs using HapPairProbsMAP and AlleleProbsMAP
    for(unsigned j = 0; j < Loci->GetNumberOfChromosomes(); ++j){
      unsigned locus = Loci->getChromosome(j)->GetLocus(0);
      for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ) {
  	SetGenotypeProbs(j, jj, locus, true); // setting last arg to true forces use of ...ProbsMAP
	++locus;
      }
    }
    //     cout << "printing genotypeprobs" << endl;
    //     for(unsigned j = 0; j < 5; ++j) { //Loci->GetNumberOfChromosomes(); ++j){
    //       for(unsigned int jj = 0; jj < Loci->GetSizeOfChromosome(j); jj++ ) {
    //   	for(int k = 0; k < NumHiddenStates*NumHiddenStates; ++k) {// loop over ancestry states
    // 	  cout << GenotypeProbs[j][jj*NumHiddenStates*NumHiddenStates+k] << " ";
    // 	}
    // 	cout << endl;

    //       }
    //     }
  }

  const unsigned numEMiters = 10;
  unsigned NumEstepiters = 10;
  double LogUnnormalizedPosterior = - numeric_limits<double>::max( );
  bool isadmixed = options.isAdmixed(0);
  if(NumGametes ==2) isadmixed = isadmixed | options.isAdmixed(1);//indicates if either gamete is admixed
  //use current parameter values as initial values

  if(isadmixed) {
    double *SumLocusAncestryHat = new double[2*NumHiddenStates];
    for(unsigned EMiter = 0; EMiter < numEMiters; ++EMiter) { // begin iteration over E and M steps
      double SumNumArrivalsHat[2] = {0,0};
      fill(SumLocusAncestryHat, SumLocusAncestryHat + 2*NumHiddenStates, 0.0);

      //E-step: fix theta and rho, sample Locus Ancestry and Number of Arrivals
      NumEstepiters *= 2;
      for (unsigned Estepiters = 0; Estepiters < NumEstepiters ; ++Estepiters) {
	if(NumHiddenStates >1){
	  ResetSufficientStats();
	  SampleHiddenStates(options);
	  SampleJumpIndicators((!options.isGlobalRho()));
	}
	vector<unsigned> SumN = getSumNumArrivals();
	vector<unsigned> SumN_X = getSumNumArrivals_X(); // element 0 should be zero in a male
	// accumulate sums
	SumNumArrivalsHat[0] += SumN[0] + SumN_X[0];
	SumNumArrivalsHat[1] += SumN[1] + SumN_X[1];
	for(int i = 0; i < 2*NumHiddenStates; ++i) {
	  SumLocusAncestryHat[i] += SumLocusAncestry[i];
	}
      }
      // set SumLocusAncestry and SumNumArrivals to their averages over current E step
      for(int i = 0; i < 2*NumHiddenStates; ++i) {
	SumLocusAncestryHat[i] /= (double)NumEstepiters;
      }
      if(options.getDisplayLevel() >2)
	cout << "E step\t";
      //     for(unsigned int g = 0; g < NumGametes; ++g) {
      //       for(int k = 0; k < NumHiddenStates; ++k) {
      // 	cout << SumLocusAncestryHat[g*NumHiddenStates + k] << "\t";
      //       }
      //     }
      // cout << "\t" << SumNumArrivals[0] << "\t" << SumNumArrivals[1] << "\t";
      SumNumArrivalsHat[0] /= (double)NumEstepiters;
      SumNumArrivalsHat[1] /= (double)NumEstepiters;

      //M-step: set mode for log rho and for theta in softmax basis - no minus ones in exponents
      // use vector rhohat and array ThetaProposal to store updates before accept/reject
      unsigned gg = 0; // indexes which admixture prior (alpha) to use: use alpha[1] if second gamete and no indadmixhiermodel
      for(unsigned g = 0; g < NumGametes; ++g) {
	if( options.isAdmixed(g) ) {
	  if(!options.isGlobalRho()) {
	    rhohat[g] = ( rhoalpha + SumNumArrivalsHat[g] ) / ( rhobeta + EffectiveL[g] );
	  }
	  if(g==1 && !options.getIndAdmixHierIndicator()) {
	    gg = 1;
	  }
	  double sum = accumulate(alpha[gg].begin(), alpha[gg].end(),0.0, std::plus<double>())
	    + accumulate(SumLocusAncestryHat+g*NumHiddenStates,
			 SumLocusAncestryHat+(g+1)*NumHiddenStates, 0.0, std::plus<double>());
	  for(int k = 0; k < NumHiddenStates; ++k) {
            if (alpha[gg][k] > 0.0)
              ThetaProposal[g][k] = (alpha[gg][k] + SumLocusAncestryHat[g*NumHiddenStates+k]) / sum;
	  }
	} else {//unadmixed gamete - set ThetaProposal to (fixed) values in thetahat, which will be something like 1,0,0
          ThetaProposal[g] = thetahat[g];
	}
      } // end loop over gametes

      // evaluate log unnormalized posterior density
    double logpriorhat =  LogPriorTheta_Softmax(ThetaProposal, options, alpha) +
      LogPriorRho_LogBasis(rhohat, options, rhoalpha, rhobeta);
    loglikhat = getLogLikelihood(options, ThetaProposal, rhohat, false);
    double LogUnnormalizedPosteriorHat  = logpriorhat + loglikhat;

    // accept update only if density increases
    if (LogUnnormalizedPosteriorHat > LogUnnormalizedPosterior) {
      if(isadmixed) {
	setAdmixtureProps(ThetaProposal);
	copy(rhohat.begin(), rhohat.end(), _rho.begin());
      }
      LogUnnormalizedPosterior = LogUnnormalizedPosteriorHat;
    }
    if(options.getDisplayLevel() >2)
      cout << "LogPrior " << logpriorhat << "\tLogLikelihood " << loglikhat << "\tLogUnNormalizedPosterior "
	   << LogUnnormalizedPosterior << endl << flush;
    } //end EM outer loop
    delete[] SumLocusAncestryHat;
  } // end block conditional on isadmixed

  // print values to file
  modefile<<setiosflags(ios::fixed)<<setprecision(3);
  modefile << getMyNumber() << "\t";
  if(!options.isGlobalRho()) {
    for (unsigned i = 0; i < NumGametes; ++i)
      modefile << _rho[i] << "\t ";
   }
  for (unsigned i = 0; i < NumGametes; ++i) {
    for (int k = 0; k < NumHiddenStates; ++k)
      modefile << Theta[i][k] << "\t ";
  }

  if(getIndex()==0 && options.getChibIndicator()){ // copy modes into hat arrays to use in Chib algorithm
    thetahat = Theta;

    // check for zeros where not specified by prior
    for(unsigned i = 0; i < NumGametes; ++i) {

      // index of alpha to use: if indadmixhiermodel use 0 for both gametes,
      // otherwise use 1 for the second gamete
      unsigned gg = options.getIndAdmixHierIndicator() ? 0 : i;
      double sum = 0.0;

      for (int k = 0; k < NumHiddenStates; ++k) {
        // if prior not zero and mode < 0.001, use 0.001
        if (alpha[gg][k] > 0.0 && thetahat[i][k] < 0.001)
          thetahat[i][k] = 0.001;
	sum += thetahat[i][k];
      }

      // re-normalize
      thetahat[i] /= sum;
    }

    copy(_rho.begin(), _rho.end(), rhohat.begin());
  }
  // compute log likelihood at posterior modes
  loglikhat = getLogLikelihood(options, thetahat, rhohat, true); //this value will be used in chib numerator
  logLikelihood.HMMisOK = true;
  logLikelihood.ready = true;
  logLikelihood.value = loglikhat;
}


void AdmixedIndividual::SampleTheta( const int iteration, double *SumLogTheta, const bclib::DataMatrix* const Outcome,
				     const DataType* const OutcomeType, const vector<double> & lambda, const int NumCovariates,
				     bclib::DataMatrix * Covariates, const vector<const double*> & beta, const PopAdmix::PopThetaType & poptheta,
				     const AdmixOptions& options, const AlphaType & alpha,
				     double DInvLink, const double dispersion, CopyNumberAssocTest& ancestryAssocTest,
				     const bool RW, const bool anneal )
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
  int K = NumHiddenStates;

  //calculate Metropolis acceptance probability ratio for proposal theta
  if(!options.getTestForAdmixtureAssociation() && getMyNumber() < Outcome->nCols()){
    RegressionType RegType;
    int NumOutcomes = Outcome->nCols();
    for( int k = 0; k < NumOutcomes; k++ ){
      if(OutcomeType[k] == Binary)RegType = Logistic; else RegType = Linear;
      logpratio += LogAcceptanceRatioForRegressionModel( RegType, options.isRandomMatingModel(), K, NumCovariates,
							 Covariates, beta[k], Outcome->get( getIndex(), k ),
							 poptheta, lambda[k]);
    }
  }

  //Accept or reject proposed value - if conjugate update and no regression model, proposal will be accepted because logpratio = 0
  Accept_Reject_Theta(logpratio, K, options.isRandomMatingModel(), RW );

  // update the value of admixture proportions used in the regression model
  if( options.getNumberOfOutcomes() > 0 )
    UpdateAdmixtureForRegression(K, NumCovariates, poptheta, options.isRandomMatingModel(), Covariates);

  if(!anneal && iteration > options.getBurnIn()){ // accumulate sums in softmax basis for calculation of posterior means
    double* a = new double[NumHiddenStates];
    for( unsigned int g = 0; g < NumGametes; g++ ){
      Theta[g].inv_softmax_gt0(a);
      transform(a, a+NumHiddenStates, SumSoftmaxTheta+g*NumHiddenStates, SumSoftmaxTheta+g*NumHiddenStates, std::plus<double>());
    }
    delete[] a;
  }
  if(!IAmUnderTest){
    for( int k = 0; k < K; k++ ){
      SumLogTheta[ k ] += log( Theta[0][k] );
      if(NumGametes==2 )
        SumLogTheta[ k ] += log( Theta[1][k] );
    }

    //increment B using new Admixture Props
    //Xcov is a vector of admixture props as covariates as in UpdateScoreForAncestry
    if(iteration >= options.getBurnIn() && options.getTestForLinkageWithAncestry()){
      double* admixtureCovars = new double[NumHiddenStates-1];
      for(int t = 0; t < NumHiddenStates-1; ++t)admixtureCovars[t] = Covariates->get(getIndex(), Covariates->nCols()-NumHiddenStates+1+t);
       ancestryAssocTest.UpdateB(DInvLink, dispersion, admixtureCovars);
      delete[] admixtureCovars;
    }
  }
}

void AdmixedIndividual::SamplePsi(const AdmixOptions& options,
                                  const cvector<double>& priormean,
                                  const cvector<double>& priorprec,
                                  bool updateSumLogPsi) {

  // In the update of the odds ratio vector psi, we skip element 0 as it
  // is 1.0 by definition, and we do a random-walk on the other elements
  for (int el = 1; el < NumHiddenStates; ++el) {

    // propose log psi from normal distribution with SD step
    const double logpsi = log(psi[el]);
    const double logpsiprop = Rand::gennor(logpsi, psistep[el]);
    const double psiprop = exp(logpsiprop);

    // force update, don't store the result
    double LogLikelihood = getLogLikelihoodXChr(options, true, false);

    // HMM probs overwritten by next indiv
    HMMIsBad(true);  // XXX needed?

    // get log likelihood at proposed values

    // set the proposed value for psi
    const double storepsi = psi[el];
    psi[el] = psiprop;

    // force update, don't store the result
    double LogLikelihoodAtProposal = getLogLikelihoodXChr(options, true, false);

    // HMM probs overwritten by next indiv
    HMMIsBad(true);  // XXX needed?

    // gaussian prior: -0.5 * tau * (logpsi - mu)^2
    const double LogPriorRatio = -0.5 * priorprec[el] * (logpsiprop - logpsi)
                                      * (logpsiprop + logpsi - 2 * priormean[el]);
    const double LogLikelihoodRatio = LogLikelihoodAtProposal - LogLikelihood;
    const double LogAccProbRatio = LogLikelihoodRatio + LogPriorRatio;

    // generic Metropolis step
    const bool accept = (LogAccProbRatio >= 0) ||
                        (log(Rand::myrand()) < LogAccProbRatio);

    if (!accept) {
      psi[el] = storepsi;
    }

    // update sampler object every w updates
    if( !(NumberOfPsiUpdates % w) )
      psistep[el] = TunePsiSampler[el].UpdateStepSize(exp(LogAccProbRatio));

    // Accumulate sum of log psi after burnin
    if (updateSumLogPsi)
      SumLogPsi[el] += logpsi;
  }
}

// ****** End Public Interface *******

double AdmixedIndividual::ProposeThetaWithRandomWalk( const AdmixOptions& options, const AlphaType & alpha ) {
  double LogLikelihoodRatio = 0.0;
  double LogPriorRatio = 0.0;

  double* a = new double[NumHiddenStates]; // should be at class scope
  bool* b = new bool[NumHiddenStates];

  //generate proposals
  for( unsigned int g = 0; g < NumGametes; g++ ){
    if(options.isAdmixed(g)){
      // inverse softmax transformation from proportions to numbers on real line that sum to 0
      for(int k = 0; k < NumHiddenStates; ++k) {
	if (Theta[g][k] > 0.0) {
	  b[k] = true; //to skip elements set to zero
	} else b[k] = false;
      }
      Theta[g].inv_softmax_gt0(a);

      DEBUG_TH_PROP(
	  fprintf( stderr, "DBG-TH-PR-1: %d %zu", getMyNumber(), g );
	  for ( int k = 0 ; k < NumHiddenStates ; ++k )
	      fprintf( stderr, " %.9lf", Theta[g][k] );
	  putc( '\n', stderr );
	  fprintf( stderr, "DBG-TH-PR-2: %d %zu", getMyNumber(), g );
	  for ( int k = 0 ; k < NumHiddenStates ; ++k )
	      fprintf( stderr, " %.9lf", a[k] );
	  putc( '\n', stderr );

	  fprintf( stderr, "PRNG POISON TEST: %.12lf\n", RNG_UNIFORM() );
      ) // end DEBUG_TH_PROP()

      //random walk step - on all elements of array a
      for(int k = 0; k < NumHiddenStates; ++k) {
	if( b[k] )
	    a[k] = Rand::gennor(a[k], step);
      }

      DEBUG_TH_PROP(
	  fprintf( stderr, "DBG-TH-PR-3: %d %zu", getMyNumber(), g );
	  for ( int k = 0 ; k < NumHiddenStates ; ++k )
	      fprintf( stderr, " %.9lf", a[k] );
	  putc( '\n', stderr );
      ) // end DEBUG_TH_PROP()

      //reverse transformation from numbers on real line to proportions
      bclib::softmax(NumHiddenStates, a, a, b);
      for (int i = 0; i < NumHiddenStates; ++i)
        if (b[i]) ThetaProposal[g][i] = a[i];

      DEBUG_TH_PROP(
	  fprintf( stderr, "DBG-TH-PR-4: %d %zu", getMyNumber(), g );
	  for ( int k = 0 ; k < NumHiddenStates ; ++k )
	      fprintf( stderr, " %.9lf", ThetaProposal[g][k] );
	  putc( '\n', stderr );
      ) // end DEBUG_TH_PROP()

      //compute contribution of this gamete to log prior ratio
      for(int k = 0; k < NumHiddenStates; ++k) {
	if( b[k] ) {
	  // prior densities must be evaluated in softmax basis
          LogPriorRatio += alpha[g][k]*(log(ThetaProposal[g][k]) -
                                        log(Theta[g][k]));
	}
      }
      //     //compute contribution of this gamete to log prior ratio
      //       LogPriorRatio += getDirichletLogDensity_Softmax(alpha[g], ThetaProposal[g) -
      // 	getDirichletLogDensity_Softmax(alpha[g], Theta[g]);
    }
    else
      ThetaProposal[g] = Theta[g];
  }// end loop over gametes

  delete[] a;
  delete[] b;

  DEBUG_TH_PROP( dumpTheta( "New proposed theta" ); )

  //get log likelihood at current parameter values - do not force update, store result of update
  LogLikelihoodRatio -= getLogLikelihood(options, false, true);
  DEBUG_TH_PROP( fprintf( stderr, "Accept-LLR: %.9lf", getLogLikelihood(options, false, true) ); )

  //get log likelihood at proposal theta and current rho - force update
  // store result in loglikelihood.tempvalue, and accumulate loglikelihood ratio
  logLikelihood.tempvalue = getLogLikelihood(options, ThetaProposal, _rho, true);
  LogLikelihoodRatio += logLikelihood.tempvalue;
  DEBUG_TH_PROP( fprintf( stderr, "  Propose-LLR: %.9lf  Ratio: %.9lf  PriorRat: %.9lf\n",
			logLikelihood.tempvalue, LogLikelihoodRatio, LogPriorRatio ); )
  return LogLikelihoodRatio + LogPriorRatio;// log ratio of full conditionals
}

// Proposes new values for individual admixture proportions
// as conjugate Dirichlet posterior conditional on prior parameter vector alpha and
// multinomial likelihood given by sampled values of ancestry at loci where jump indicator xi is 1 (SumLocusAncestry)
// proposes new values for both gametes if random mating model
void AdmixedIndividual::ProposeTheta(const AdmixOptions& options, const AlphaType & alpha,
			      int* sumLocusAncestry, int* sumLocusAncestry_X){
  size_t K = NumHiddenStates;
  if( isHaploid || options.isRandomMatingModel() ){ //random mating model
    for( unsigned int g = 0; g < NumGametes; g++ ) {
      if(options.isAdmixed(g)) {
	unsigned gg = g; // index of which prior (alpha) to use - use alpha[1] if second gamete and no indadmixhiermodel
	if(options.getIndAdmixHierIndicator()) gg = 0;
	for(size_t k = 0; k < K; ++k) { // parameters of conjugate Dirichlet update
	  dirparams[k] = alpha[gg][k] + (double)(sumLocusAncestry[k + K*g] + sumLocusAncestry_X[g*K + k]);
	  // if male, paternal elements (0 to K-1) of sumLocusAncestry_X will be zero
	}
	//generate proposal theta from Dirichlet with parameters dirparams
	Rand::gendirichlet(K, dirparams, ThetaProposal[g]);
      } else { // gamete unadmixed
        ThetaProposal[g] = Theta[g];
      }
    } // end loop over gametes

  } else { //assortative mating model
    for(size_t k = 0; k < K; ++k) {
      dirparams[k] = alpha[0][k] + double(sumLocusAncestry[k] + sumLocusAncestry_X[k] +
					  sumLocusAncestry[K + k] + sumLocusAncestry_X[K + k]);
    }
    Rand::gendirichlet(K, dirparams, ThetaProposal[0]);
  }
}

double AdmixedIndividual::LogAcceptanceRatioForRegressionModel( RegressionType RegType, bool RandomMatingModel,
							 int NumHiddenStates, int NumCovariates,
							 const bclib::DataMatrix * Covariates, const double* beta,
							 const double Outcome, const PopAdmix::PopThetaType & poptheta, const double lambda) {
  // returns log of ratio of likelihoods of new and old values of population admixture
  // in regression models.  individual admixture theta is standardized about the mean poptheta calculated during burn-in.
  double logprobratio = 0.0, XBeta = 0.0, currentXBeta = 0.0;
  vector<double> avgtheta(NumHiddenStates);avgtheta[0] = 0.0;
  vector<double> currentavgtheta(NumHiddenStates);currentavgtheta[0] = 0.0;
  if( RandomMatingModel && NumGametes==2)
    for(int k = 1;k < NumHiddenStates; ++k){
      avgtheta[k] = (ThetaProposal[0][k] + ThetaProposal[1][k])/ 2.0 - poptheta[k];
      currentavgtheta[k] = (Theta[0][k] + Theta[1][k])/ 2.0 - poptheta[k];
    }
  else
    for(int k = 1;k < NumHiddenStates; ++k){
      avgtheta[k] = ThetaProposal[0][k] - poptheta[k];
      currentavgtheta[k] = Theta[0][k] - poptheta[k];
    }

  for( int jj = 0; jj < NumCovariates - NumHiddenStates + 1; jj++ ){
    XBeta += Covariates->get( getIndex(), jj ) * beta[jj];
    currentXBeta += Covariates->get( getIndex(), jj ) * beta[jj];
  }
  for( int k = 1; k < NumHiddenStates; k++ ){
    XBeta += avgtheta[ k ] * beta[NumCovariates - NumHiddenStates + k ];
    currentXBeta += currentavgtheta[ k ] * beta[NumCovariates - NumHiddenStates + k];
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
void AdmixedIndividual::UpdateAdmixtureForRegression( int NumHiddenStates, int NumCovariates,
                                               const PopAdmix::PopThetaType & poptheta, bool RandomMatingModel,
					       bclib::DataMatrix *Covariates) {
  vector<double> avgtheta(NumHiddenStates);
  if(RandomMatingModel && NumGametes==2)//average over gametes
    for(int k = 0; k < NumHiddenStates; ++k) avgtheta[k] = (Theta[0][k] + Theta[1][k]) / 2.0;
  else
    for(int k = 0; k < NumHiddenStates; ++k) avgtheta[k] = Theta[0][k];
  for( int k = 1; k < NumHiddenStates ; k++ )
    Covariates->set( getIndex(), NumCovariates - NumHiddenStates + k, avgtheta[ k ] - poptheta[ k ] );
}

void AdmixedIndividual::Accept_Reject_Theta( double logpratio, /*bool xdata, */ int NumHiddenStates, bool RandomMatingModel, bool RW ) {
  // Metropolis update for admixture proportions theta, taking log of acceptance probability ratio as argument
  // uses log ratio because this is easier than ratio to calculate for linear regression model
  // if conjugate update and no regression model, logpratio remains set to 0, so all proposals are accepted
  bool test = true;
  bool accept = false;
  double AccProb = 1.0;
  // loop over populations: if any element of proposed Dirichlet parameter vector is too small, reject update without test step
  for( int k = 0; k < NumHiddenStates; k++ ) {
    if ( Theta[0][k] > 0.0 && ThetaProposal[0][k] < 0.0001 ) {
      test = false;
    }
    else if ( NumGametes==2 && RandomMatingModel && Theta[1][k] > 0.0 && ThetaProposal[1][k] < 0.0001 ) {
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

  DEBUG_TH_PROP( fprintf( stderr, "Accept: %s\n", accept ? "yes" : "no" ); )

  if(accept) { // set proposed values as new values
    setAdmixtureProps(ThetaProposal);
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

void AdmixedIndividual::resetStepSizeApproximator(int k) {
  ThetaTuner.resetStepSizeApproximator(k);
}

void AdmixedIndividual::UpdateHMMInputs(unsigned int j, const Options& options,
                                        const AdmixtureProportions& theta,
                                        const RhoType& rho) {
  //Updates inputs to HMM for chromosome j
  //also sets Diploid flag in Chromosome (last arg of SetStateArrivalProbs)
  const bool diploid = !isHaploid && (j!=X_posn || SexIsFemale);
  const bool isRandomMating = options.isRandomMatingModel();
  Chromosome* C = Loci->getChromosome(j);

  C->HMM->SetGenotypeProbs(GenotypeProbs[j], GenotypesMissing[j]);
  //   if(!SexIsFemale) cout << "chr " << j << " " << X_posn << " ";
  //   if(!diploid) cout << "haploid GenotypeProbs set for chromosome " << j << "\n";

  if(!options.isGlobalRho()){
    //set locus correlation, f, if individual- or gamete-specific rho
    C->SetLocusCorrelation(rho, isRandomMating);
  }

  // set StateArrivalProbs in HMM
  // in the haploid case in a random mating model, pass the pointer to the
  // maternal admixture proportions
  const int maternalShift = (diploid || !isRandomMating) ? 0 : NumHiddenStates;
  if (C->isXChromosome())
    C->HMM->SetStateArrivalProbs(theta.flatXChromosome(psi) + maternalShift,
                                 isRandomMating, diploid);
  else
    C->HMM->SetStateArrivalProbs(theta.flat() + maternalShift,
                                 isRandomMating, diploid);

  logLikelihood.HMMisOK = false;//because forward probs in HMM have been changed
}

void AdmixedIndividual::SampleRho(const AdmixOptions& options, double rhoalpha, double rhobeta,
			   bool updateSumLogRho) {
  vector<unsigned> sumNumArrivals = getSumNumArrivals();
  vector<unsigned> sumNumArrivals_X = getSumNumArrivals_X();
  // rho_X is set to 0.5*rho - equivalent to setting effective length of X chromosome as half length in Morgans
  // conjugate gamma update for rho includes arrivals on X chromosome, and 0.5 * length of X chromosome
  if(isHaploid || options.isRandomMatingModel() ) {
    // SumNumArrivals_X has length 2, and SumNumArrivals_X[0] remains fixed at 0 if male
    for( unsigned int g = 0; g < NumGametes; g++ )if(options.isAdmixed(g)) {
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

//********************** Score Tests ***************************
void AdmixedIndividual::UpdateScores(const AdmixOptions& options, bclib::DataMatrix *Outcome, bclib::DataMatrix *Covariates,
				     const vector<bclib::Regression*> & R, AffectedsOnlyTest& affectedsOnlyTest, CopyNumberAssocTest& ancestryAssocTest){
//merge with updatescoretests
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    Chromosome* C = Loci->getChromosome(j);
    // update of forward probs here is unnecessary if SampleTheta was called and proposal was accepted
      //Update Forward/Backward probs in HMM
      if( !logLikelihood.HMMisOK ) {
	UpdateHMMInputs(j, options, Theta, _rho);
      }
      //update of score tests for linkage with ancestry requires update of backward probs
      double* admixtureCovars = 0;
      if(options.getTestForLinkageWithAncestry()) {
	admixtureCovars = new double[NumHiddenStates-1];
	for(int t = 0; t < NumHiddenStates-1; ++t)admixtureCovars[t] = Covariates->get(getIndex(), Covariates->nCols()-NumHiddenStates+1+t);
      }
      UpdateScoreTests(options, admixtureCovars, Outcome, C, R, affectedsOnlyTest, ancestryAssocTest);
      if(options.getTestForLinkageWithAncestry()) {
	delete[] admixtureCovars;
      }
  } //end chromosome loop
}

void AdmixedIndividual::UpdateScoreTests(const AdmixOptions& options, const double* admixtureCovars, bclib::DataMatrix *Outcome,
					 Chromosome* chrm, const vector<bclib::Regression*> R, AffectedsOnlyTest& affectedsOnlyTest, CopyNumberAssocTest& ancestryAssocTest){
  bool IamAffected = false;
  try {
    if( options.getTestForAffectedsOnly()){
      //determine which regression is logistic, in case of 2 outcomes
      unsigned col = 0;
      if(options.getNumberOfOutcomes() >1 && R[0]->getRegressionType()!=Logistic )col = 1;
      //check if this individual is affected
      if(options.getNumberOfOutcomes() == 0 || Outcome->get(getIndex(), col) == 1) IamAffected = true;
    }

    //we don't bother computing scores for the first population when there are two
    int KK = NumHiddenStates,k0 = 0;
    if(NumHiddenStates == 2) {KK = 1;k0 = 1;}

    int locus;
    for( unsigned int jj = 0; jj < chrm->GetSize(); jj++ ){
      locus = chrm->GetLocus(jj);
      //retrieve AncestryProbs from HMM
      // DDF says: should use a const-reference here:
      vector<vector<double> > AProbs =
	chrm->getHiddenStateCopyNumberProbs(!isHaploid && (!chrm->isXChromosome() || SexIsFemale),  jj );

      // diploid is true if (not isHaploid) and (female or not X chr)
      bool diploid = !isHaploid && (SexIsFemale ||
                                    (Loci->GetChrNumOfLocus(locus) != X_posn));

      //Update affecteds only scores
      if(IamAffected){
        // if random mating model and X chromosome in male, should pass the
        // maternal gamete X chromosome admixture
	if(options.isRandomMatingModel() && !SexIsFemale && (Loci->GetChrNumOfLocus(locus) == X_posn)) {
          affectedsOnlyTest.Update(locus, k0,
                                   Theta.flatXChromosome(psi) + NumHiddenStates,
                                   options.isRandomMatingModel(),
                                   diploid, AProbs);
        }
        // if X chromosome
        else if (Loci->GetChrNumOfLocus(locus) == X_posn) {
          affectedsOnlyTest.Update(locus, k0, Theta.flatXChromosome(psi),
                                   options.isRandomMatingModel(),
                                   diploid, AProbs);
	} else { // just pass Theta
          affectedsOnlyTest.Update(locus, k0, Theta.flat(), options.isRandomMatingModel(),
                                   diploid, AProbs);
	}
      }

      //update ancestry score tests
      if( options.getTestForLinkageWithAncestry() ){
	ancestryAssocTest.Update(locus, admixtureCovars, R[0]->getDispersion(),
				 Outcome->get(getIndex(), 0) - R[0]->getExpectedOutcome(getIndex()),
				 R[0]->DerivativeInverseLinkFunction(getIndex()),
                                 diploid, AProbs);
      }

      if ( AncestryProbs ) {
	// accumulate hidden state probabilities at given locus
	for( int k = 0; k < NumHiddenStates; k++ ){
	  for( unsigned a = 0; a < 3; ++a) {
	    SumProbs[locus*NumHiddenStates*3 + k*3 + a] += AProbs[a][k];
	  }
	}
      }
      ++locus;
    }//end within-chromosome loop
  } catch (string msg) {
    cout << msg;
  }
}

// this function does two things:
// 1. calculates log-likelihood at thetahat, rhohat, allelefreqsMAP
// 2. calculates log prior at same values
void AdmixedIndividual::setChibNumerator(const AdmixOptions& options, const AlphaType &alpha,
				  double rhoalpha, double rhobeta, chib *MargLikelihood, AlleleFreqs* A) {

  // 1. pass value of log-likelihood at MAP parameter values, calculated after finding posterior modes, to chib
  MargLikelihood->setLogLikelihood(loglikhat);

  // 2. calculate log prior at MAP parameter values
  double LogPrior = LogPriorTheta_Softmax(thetahat, options, alpha)
    + LogPriorRho_LogBasis(rhohat, options, rhoalpha, rhobeta);
  if( A->IsRandom() ) {
    double LogPriorFreqs = 0.0;
    for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
      for( int k = 0; k < NumHiddenStates; k++ ){
	LogPriorFreqs += bclib::getDirichletLogDensity( A->GetPriorAlleleFreqs(j, k), A->getAlleleFreqsMAP(j,k) );
      }
    }
    LogPrior += LogPriorFreqs;
  }
  MargLikelihood->setLogPrior(LogPrior);
}

void AdmixedIndividual::updateChib(const AdmixOptions& options, const AlphaType &alpha,
				   double rhoalpha, double rhobeta,
				   // double *thetahat, vector<double> &rhohat,
				   chib *MargLikelihood, AlleleFreqs* A){
  // *** After BurnIn *** - accumulate samples of posterior ordinates for theta, rho, allelefreqs separately
  double LogPosterior = 0.0;
  double LP = 0.0;
  if( NumHiddenStates > 1 ){
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
      for( int k = 0; k < NumHiddenStates; k++ ){
	vector<double> args = A->GetPriorAlleleFreqs(j,k);
	vector<int> counts = A->GetAlleleCounts(j,k);
	//print_vector(args);
	transform(counts.begin(), counts.end(), args.begin(), args.begin(), plus<double>());//PriorAlleleFreqs + AlleleCounts
	//print_vector(args);
	LogPosteriorFreqs += bclib::getDirichletLogDensity( args, A->getAlleleFreqsMAP(j, k) );//LogPosterior for Allele Freqs
      }
    }
    LogPosterior += LogPosteriorFreqs;
    logPosterior[2].push_back( LogPosteriorFreqs  );
  }
  MargLikelihood->addLogPosteriorObs( LogPosterior );
}

// TODO: fix these two functions to work with assortative mating
double AdmixedIndividual::LogPriorTheta_Softmax(const AdmixtureProportions& theta,
                                                const AdmixOptions& options,
                                                const AlphaType& alpha) const {
  double LogPrior=0.0;
  for(unsigned g = 0; g < NumGametes; ++g) { //loop over gametes
    if( options.isAdmixed(g) ){
      LogPrior += bclib::getDirichletLogDensity_Softmax(alpha[g].getVector_unsafe(), theta[g]);
    }
  }
  return LogPrior;
}

double AdmixedIndividual::LogPosteriorTheta_Softmax(const AdmixOptions& options,
                                                    const AdmixtureProportions& theta,
                                                    const AlphaType& alpha) const {

  // calculates log full conditional at theta, conditional on realized locus ancestry states and jump indicators
  double LogPosterior = 0.0;
  vector<double> alphaparams(NumHiddenStates); // , alphaparams1(NumHiddenStates);  // to be set to alpha + SumLocusAncestry
  for(unsigned g = 0; g < NumGametes; ++g) { //loop over gametes
    if( options.isAdmixed(g) ) {
      transform(alpha[g].begin(), alpha[g].end(), SumLocusAncestry + g*NumHiddenStates,
		alphaparams.begin(), std::plus<double>());
      LogPosterior += bclib::getDirichletLogDensity_Softmax(alphaparams, theta[g]);
    }
  }
  //   if(  options->isAdmixed(1) ){//admixed second gamete
  //     transform(alpha[1].begin(), alpha[1].end(), SumLocusAncestry+NumHiddenStates, alphaparams1.begin(), std::plus<double>());
  //     LogPosterior += getDirichletLogDensity_Softmax(alphaparams1, theta+NumHiddenStates);
  //   }
  return LogPosterior;
}

double AdmixedIndividual::LogPriorRho_LogBasis( const RhoType & rho, const AdmixOptions& options, double rhoalpha,
					       double rhobeta) const {
  // computes log prior density in log rho basis at supplied parameter values
  double LogPrior=0.0;
  if( options.isRandomMatingModel() ) {
    for(unsigned g = 0; g < NumGametes; ++g) { //loop over gametes
      if( options.isAdmixed(g) ){
	LogPrior += bclib::getGammaLogDensity_LogBasis( rhoalpha, rhobeta, rho[g] );
      }
    }
  } else { // assortative mating: rho assumed same on both gametes
    LogPrior = bclib::getGammaLogDensity_LogBasis( rhoalpha, rhobeta, rho[0] );
  }
  return LogPrior;
}

double AdmixedIndividual::LogPosteriorRho_LogBasis(const AdmixOptions & options, const RhoType & rho,
					    double rhoalpha, double rhobeta)const{
  // calculates log full conditional density at sum-intensities rho, conditional on realized number of arrivals
  double LogPosterior = 0.0;
  vector<unsigned> SumN = getSumNumArrivals();
  vector<unsigned> SumN_X = getSumNumArrivals_X();
  if(options.isRandomMatingModel() ) { // SumNumArrivals_X has length 2, and SumNumArrivals_X[0] remains fixed at 0 if male
    for( unsigned int g = 0; g < NumGametes; g++ ) {
      if(options.isAdmixed(g)) {
	LogPosterior += bclib::getGammaLogDensity_LogBasis( rhoalpha + (double)(SumN[g] + SumN_X[g]), rhobeta + EffectiveL[g], rho[g] );
      }
    }
  } else {//assortative mating, rho assumed same on both gametes
    if(options.isAdmixed(0)) {
      LogPosterior+= bclib::getGammaLogDensity_LogBasis( rhoalpha + (double)(SumN[0] + SumN[1] + SumN_X[0] + SumN_X[1]),
					 rhobeta + EffectiveL[0] + EffectiveL[1], rho[0] );
    }
  }
  return LogPosterior;
}

// these three functions can be used to monitor stability of estimates of components of posterior ordinate
double AdmixedIndividual::getLogPosteriorTheta()const{
  if(NumHiddenStates > 1){
    vector<double>::const_iterator max = max_element(logPosterior[0].begin(), logPosterior[0].end());
    return bclib::AverageOfLogs(logPosterior[0], *max);
  }
  else return 0.0;
}
double AdmixedIndividual::getLogPosteriorRho()const{
  if(NumHiddenStates > 1){
    vector<double>::const_iterator max = max_element(logPosterior[1].begin(), logPosterior[1].end());
    return bclib::AverageOfLogs(logPosterior[1], *max);
  }
  else return 0.0;
}
double AdmixedIndividual::getLogPosteriorAlleleFreqs()const{
  vector<double>::const_iterator max = max_element(logPosterior[2].begin(), logPosterior[2].end());
  return bclib::AverageOfLogs(logPosterior[2], *max);
}



// ====== DEBUGGING METHODS (overridden from PedBase) ======
#if PEDBASE_DEBUG_METHODS
    void AdmixedIndividual::dumpTheta( const char * prefix ) const
	{
	const size_t K = NumHiddenStates;
	for ( unsigned int g = 0 ; g < NumGametes ; ++g )
	    {
	    fprintf( stderr, "%s(%d) theta[%d]: ", prefix, getMyNumber(), g );
	    for ( size_t x = 0 ; x < K ; ++x )
		fprintf( stderr, " %18.15lf", Theta[g][x] );
	    fprintf( stderr, "\n" );

	    fprintf( stderr, "	theta-prop: " );
	    for ( size_t x = 0 ; x < K ; ++x )
		fprintf( stderr, " %18.15lf", ThetaProposal[g][x] );
	    fprintf( stderr, "\n" );
	    }
	}
#endif
