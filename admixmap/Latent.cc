/** 
 *   ADMIXMAP
 *   Latent.cc 
 *   Class to hold and update population admixture and sumintensities parameters and their priors
 *   Copyright (c) 2002 - 2006 LSHTM
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 */
#include "Latent.h"
#include "Chromosome.h"
#include "functions.h"
#include <algorithm>
#include <numeric>
#include "gsl/gsl_math.h"

using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

Latent::Latent( AdmixOptions* op, const Genome* const loci)
{
  options = 0;
  rhoalpha = 0.0;
  rhobeta = 0.0;
  SumLogRho = 0.0;
  options = op;
  Loci = loci;
  poptheta = 0;
  globaltheta = 0;
  globalthetaproposal = 0;
  rho.push_back(0.0);
}

void Latent::Initialise(int Numindividuals, const std::string* const PopulationLabels, LogWriter &Log){
  Log.setDisplayMode(On);
  K = options->getPopulations();

  //ergodic average of population admixture, which is used to centre 
  // the values of individual admixture in the regression model  
  poptheta = new double[ K ];
  for( int i = 0; i < K; i++ ) poptheta[i] = 0.0;

  // ** Initialise population admixture distribution Dirichlet parameters alpha **
  alpha = options->getInitAlpha();
  SumAlpha.resize( K );
  if(!options->getIndAdmixHierIndicator())  copy(alpha[0].begin(), alpha[0].end(), SumAlpha.begin());

  if(K > 1){
    // ** set up sampler for alpha **
    // should be able to pass initial step size to the sampler
    unsigned obs = Numindividuals;
    if( options->isRandomMatingModel() ){
      obs *= 2;//for 2 gametes per individual
    } 
    PopAdmixSampler.SetSize( obs, K );

    //initialise global admixture proportions
    if(options->getHapMixModelIndicator()){
      globaltheta = new double[K];
      globalthetaproposal = new double[K];
      fill(globaltheta, globaltheta+K, 1.0/(double)K);
      ThetaTuner.SetParameters(1.0 /*<-initial stepsize on softmax scale*/, 0.00, 10.0, 0.44);
    }

    // ** get prior on sum-of-intensities parameter rho or on rate parameter of its population distribution
    rhoalpha = options->getRhoalpha();
    if(/*options->getHapMixModelIndicator() ||*/ (options->getIndAdmixHierIndicator() && !options->isGlobalRho() )){
      // get prior on rate parameter beta and initialize it at prior mean
      rhobeta0 = options->getRhobetaShape();
      rhobeta1 = options->getRhobetaRate();
      rhobeta = rhobeta0 / rhobeta1;
      if(options->getHapMixModelIndicator()){
	//initialise rho vector
	for(unsigned j = 0; j < Loci->GetNumberOfCompositeLoci()-1; ++j)
	  rho.push_back(rhoalpha / rhobeta);
      }
    }
    else{
      rhobeta = options->getRhobeta();
      if( options->isGlobalRho()){
	// set up sampler for global variable
	rho[0] = rhoalpha / rhobeta ;//initialise global sumintensities parameter at prior mean for globalrho
	// ** set up TuneRW object for global rho updates **
	NumberOfUpdates = 0;
	w = 1;
	step0 = 1.0; // sd of proposal distribution for log rho
	//need to choose sensible value for this initial RW sd
	step = step0;
	TuneRhoSampler.SetParameters( step0, 0.01, 10, 0.44);
      }
    }

    // ** Open paramfile **
    if ( options->getIndAdmixHierIndicator()){
      Log.setDisplayMode(Quiet);
      if( strlen( options->getParameterFilename() ) ){
	outputstream.open( options->getParameterFilename(), ios::out );
	if( !outputstream )
	  {
	    Log.setDisplayMode(On);
	    Log << "ERROR: Couldn't open paramfile\n";
	    exit( 1 );
	  }
	else{
	  Log << "Writing population-level parameters to " << options->getParameterFilename() << "\n";
	  InitializeOutputFile(PopulationLabels);
	}
      }
      else{
	Log << "No paramfile given\n";
      }
    }
  }//end if Populations > 1
}

Latent::~Latent()
{
  delete[] poptheta;
  delete[] globaltheta;
  delete[] globalthetaproposal;
}

void Latent::UpdateGlobalTheta(int iteration, IndividualCollection* individuals, Chromosome** C){
  if(options->getHapMixModelIndicator()){
    if(!(iteration%2))ConjugateUpdateGlobalTheta(individuals->getSumLocusAncestry(options->getPopulations()));
    else
      UpdateGlobalThetaWithRandomWalk(individuals, C); 
  }
  individuals->setAdmixtureProps(globaltheta, options->getPopulations());//shouldn't be necessary
}

void Latent::UpdatePopAdmixParams(int iteration, const IndividualCollection* const individuals, LogWriter &Log, bool anneal=false)
 {
   
   if( options->getPopulations() > 1 && individuals->getSize() > 1 &&
       options->getIndAdmixHierIndicator() ){

   // ** Sample for population admixture distribution Dirichlet parameters, alpha **
   // For a model in which the distribution of individual admixture in the population is a mixture
   // of components, we will have one Dirichlet parameter vector for each component, 
   // updated only from those individuals who belong to the component
   
     //sample alpha conditional on individual admixture proportions
     PopAdmixSampler.Sample( individuals->getSumLogTheta(), &alpha[0] );
     copy(alpha[0].begin(), alpha[0].end(), alpha[1].begin()); // alpha[1] = alpha[0]

  }
   // ** accumulate sum of Dirichlet parameter vector over iterations  **
   transform(alpha[0].begin(), alpha[0].end(), SumAlpha.begin(), SumAlpha.begin(), std::plus<double>());//SumAlpha += alpha[0];
   
   if( iteration < options->getBurnIn() && options->getPopulations() > 1) {
     // accumulate ergodic average of population admixture, which is used to centre 
     // the values of individual admixture in the regression model
     double sum = accumulate(SumAlpha.begin(), SumAlpha.end(), 0.0);
     if(options->getNumberOfOutcomes() > 0)for( int j = 0; j < options->getPopulations(); j++ )poptheta[j] = SumAlpha[j] / sum;
   }
   
   if( iteration == options->getBurnIn() && options->getPopulations() > 1) {
     if(options->getNumberOfOutcomes() > 0){
       Log.setDisplayMode(Off);
       Log << "Individual admixture centred in regression model around: ";
       for(int i = 0; i < options->getPopulations(); ++i)Log << poptheta[i] << "\t";
       Log << "\n";
     }
     fill(SumAlpha.begin(), SumAlpha.end(), 0.0);
   }
   
   if( !anneal && iteration > options->getBurnIn() && options->getPopulations() > 1 ){
     // accumulate sum of log of sumintensities after burnin.
     if(options->isGlobalRho()) SumLogRho += log(rho[0]);
     else SumLogRho += log(rhoalpha) - log(rhobeta);

   }
}

void Latent::UpdateSumIntensities(const IndividualCollection* const IC, Chromosome **C) {
  if( options->isGlobalRho() || options->getHapMixModelIndicator()) { // update rho with random walk MH
    vector<double> rhoprop(rho.size(), 0.0);
    double LogLikelihood = 0.0;
    double LogLikelihoodAtProposal = 0.0;
    double LogLikelihoodRatio = 0.0;
    double LogPriorRatio = 0.0;
    double LogAccProbRatio = 0.0;
    bool accept = false;

    NumberOfUpdates++;
    for(unsigned j = 0; j < rho.size(); ++j)
      rhoprop[j] = exp(gennor(log(rho[j]), step)); // propose log rho from normal distribution with SD step
    
    //get log likelihood at current parameter values, annealed if this is an annealing run
    for(int i = 0; i < IC->getSize(); ++i) {
      Individual* ind = IC->getIndividual(i);
      ind->HMMIsBad(true);//to force HMM update
      LogLikelihood += ind->getLogLikelihood(options, C, false, true); // don't force update, store result if updated
      ind->HMMIsBad(true); // HMM probs overwritten by next indiv, but stored loglikelihood still ok
   }
     // set ancestry correlations using proposed value of sum-intensities
    // value for X chromosome set to half the autosomal value 
    for( unsigned int j = 0; j < Loci->GetNumberOfChromosomes(); j++ ) {
      C[j]->SetLociCorr(rhoprop);
    }
    //get log HMM likelihood at proposal rho and current admixture proportions
    for(int i = 0; i < IC->getSize(); ++i) {
      Individual* ind = IC->getIndividual(i);
      LogLikelihoodAtProposal += ind->getLogLikelihood(options, C, true, false); // force update, do not store result 
      ind->HMMIsBad(true); // set HMM probs as bad but stored log-likelihood is still ok
      // line above should not be needed for a forced update with result not stored
    }
    LogLikelihoodRatio = LogLikelihoodAtProposal - LogLikelihood;

    //compute prior ratio
    LogPriorRatio = 0.0;
    for(unsigned j = 0; j < rho.size(); ++j)
      LogPriorRatio += getGammaLogDensity(rhoalpha, rhobeta, rhoprop[j]) - getGammaLogDensity(rhoalpha, rhobeta, rho[j]);
    LogAccProbRatio = LogLikelihoodRatio + LogPriorRatio; 

    // generic Metropolis step
    if( LogAccProbRatio < 0 ) {
      if( log(myrand()) < LogAccProbRatio ) accept = true;
    } else accept = true;  
    
    if(accept) {
      copy(rhoprop.begin(), rhoprop.end(), rho.begin());
      for(int i = 0; i < IC->getSize(); ++i){
	Individual* ind = IC->getIndividual(i);
	ind->storeLogLikelihood(false); // store log-likelihoods calculated at rhoprop, but do not set HMM probs as OK 
      }
    } else { 
      // restore ancestry correlations in Chromosomes using original value of sum-intensities
      for( unsigned int j = 0; j < Loci->GetNumberOfChromosomes(); j++ )
      C[j]->SetLociCorr(rho);
    } // stored loglikelihoods are still ok

    //update sampler object every w updates
    if( !( NumberOfUpdates % w ) ){
      step = TuneRhoSampler.UpdateStepSize( exp(LogAccProbRatio) );  
    }

  }//end if global rho model

  else { //non global rho model
    if(IC->getSize()>1 && options->getIndAdmixHierIndicator() ) { // >1 individual and hierarchical model
      // update scale parameter of gamma distribution of sumintensities in population 
      if( options->isRandomMatingModel() )
	rhobeta = gengam( 2*rhoalpha * IC->getSize() + rhobeta0, IC->GetSumrho() + rhobeta1 );
      else
	rhobeta = gengam( rhoalpha* IC->getSize() + rhobeta0, IC->GetSumrho() + rhobeta1 );
    } // otherwise do not update rhobeta
  }
}

void Latent::ConjugateUpdateGlobalTheta(const vector<int> sumLocusAncestry){
  //sumlocusancestry is summed over all individuals
  size_t K = options->getPopulations();
  double dirparams[K];
  for(size_t k = 0; k < K; ++k) {
    dirparams[k] = alpha[0][k] + sumLocusAncestry[k] + sumLocusAncestry[k + K];
  }
  gendirichlet(K, dirparams, globaltheta );
}

void Latent::UpdateGlobalThetaWithRandomWalk(IndividualCollection* IC, Chromosome** C) {
  double LogLikelihoodRatio = 0.0;
  double LogPriorRatio = 0.0;
  double logpratio = 0.0;
  
  //generate proposals
  // inverse softmax transformation from proportions to numbers on real line that sum to 0
  bool* b = new bool[K];
  double* a = new double[K]; // should be at class scope
  for(int k = 0; k < K; ++k) {
    if(globaltheta[k] > 0.0) {
      b[k] = true; //to skip elements set to zero
    } else b[k] = false;
  }
  inv_softmax(K, globaltheta, a, b);
  //random walk step - on all elements of array a
  for(int k = 0; k < K; ++k) {
    if( b[k] ) a[k] = gennor(a[k], thetastep);  
  }
  //reverse transformation from numbers on real line to proportions 
  softmax(K, globalthetaproposal, a, b);
  //compute contribution of this gamete to log prior ratio
  for(int k = 0; k < K; ++k) {
    if( b[k] ) { 
      // prior densities must be evaluated in softmax basis
      LogPriorRatio += alpha[0][k]*( log(globalthetaproposal[k]) - log(globaltheta[k]) ); 
    }
  }
  delete[] a;
  delete[] b; 

  vector<double> dummyrho(1);dummyrho[0] = rho[0];
  for(int i = 0; i < IC->getSize(); ++i){
    //get log likelihood at current parameter values - do not force update, store result of update
    LogLikelihoodRatio -= IC->getIndividual(i)->getLogLikelihood(options, C, false, true); 
    
    //get log likelihood at proposal theta and current rho - force update 
    LogLikelihoodRatio += IC->getIndividual(i)->getLogLikelihood(options, C, globaltheta, globaltheta, dummyrho, dummyrho, true);
  }
  IC->HMMIsBad(true);

  logpratio = LogLikelihoodRatio + LogPriorRatio;// log ratio of full conditionals
  Accept_Reject_Theta(logpratio, K);
}

void Latent::Accept_Reject_Theta( double logpratio, int Populations) {
  // Metropolis update for admixture proportions theta, taking log of acceptance probability ratio as argument
  bool test = true;
  bool accept = false;
  double AccProb = 1.0; 
  // loop over populations: if any element of proposed parameter vector is too small, reject update without test step
  for( int k = 0; k < Populations; k++ ) {
    if( globaltheta[ k ] > 0.0 && globalthetaproposal[ k ] < 0.0001 ) {
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
    copy(globalthetaproposal, globalthetaproposal+K, globaltheta);
    //possibly need to reset loglikelihood in Individuals
  } 
  
  //update step size in tuner object every w updates
  if( !( NumberOfUpdates % w ) ) {
    thetastep = ThetaTuner.UpdateStepSize( AccProb );
  }
}

// void Latent::UpdateSumIntensities(){


// }

void Latent::InitializeOutputFile(const std::string* const PopulationLabels)
{
  // Header line of paramfile
  //Pop. Admixture
  for( int i = 0; i < options->getPopulations(); i++ ) {
    outputstream << "\""<<PopulationLabels[i] << "\" ";
  }
  //SumIntensities
  if( options->isGlobalRho() ) outputstream << "sumIntensities\t";
  else outputstream << "sumIntensities.mean\t";
  outputstream << endl;
}

void Latent::OutputErgodicAvg( int samples, std::ofstream *avgstream)
{
  for( int j = 0; j < options->getPopulations(); j++ ){
    avgstream->width(9);
    *avgstream << setprecision(6) << SumAlpha[j] / samples << "\t";
  }
  avgstream->width(9);
  *avgstream << setprecision(6) << exp(SumLogRho / samples) << "\t";
}
void Latent::OutputParams(ostream* out){
  for( int j = 0; j < options->getPopulations(); j++ ){
    out->width(9);
    if(options->getHapMixModelIndicator())
      (*out) << setprecision(6) << globaltheta[ j ] << "\t";
    else
      (*out) << setprecision(6) << alpha[0][ j ] << "\t";
  }
  
  out->width(9);
  if( !options->isGlobalRho() )
    (*out) << setprecision(6) << rhoalpha / rhobeta  << "\t";
  else
    (*out) << setprecision(6) << rho[0] << "\t";
}

void Latent::OutputParams(int iteration, LogWriter &Log){
  //output initial values to logfile
  if( iteration == -1 )
    {
      Log.setDisplayMode(Off);
      Log.setPrecision(6);
      for( int j = 0; j < options->getPopulations(); j++ ){
	//Log->width(9);
	Log << alpha[0][j];
      }
      
      if( !options->isGlobalRho() )
	Log << rhoalpha / rhobeta;
      else
	Log << rho[0];
    }
  //output to screen
  if( options->getDisplayLevel() > 2 )
    {
      OutputParams(&cout);
    }
  //Output to paramfile after BurnIn
  if( iteration > options->getBurnIn() ){
    OutputParams(&outputstream);
    outputstream << endl;
  }
}

const vector<double > &Latent::getalpha0()const{
  return alpha[0];
}
const std::vector<vector<double> > &Latent::getalpha()const{
  return alpha;
}

double Latent::getrhoalpha()const{
  return rhoalpha;
}
double Latent::getrhobeta()const{
  return rhobeta;
}
double Latent::getglobalrho()const{
  return rho[0];
}
const vector<double> &Latent::getrho()const{
  return rho;
}
double Latent::getSumLogRho()const{
  return SumLogRho;
}
const double *Latent::getpoptheta()const{
  return poptheta;
}

void Latent::printAcceptanceRates(LogWriter &Log) {
  if(options->getHapMixModelIndicator()){
    Log << "Expected acceptance rate in global admixture sampler: "
	<< ThetaTuner.getExpectedAcceptanceRate()
	<< "\nwith final step size of "
	<< thetastep << "\n";
  }
  else{
    Log << "Expected acceptance rate in population admixture sampler: "
	<< PopAdmixSampler.getExpectedAcceptanceRate()
	<< "\nwith final step size of "
	<< PopAdmixSampler.getStepSize() << "\n";
  }
}
double Latent::getRhoSamplerAccRate()const{
  return TuneRhoSampler.getExpectedAcceptanceRate();
}

double Latent::getRhoSamplerStepsize()const{
  return step;
}
