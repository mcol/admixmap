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
#include "gsl/gsl_specfunc.h"

using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

Latent::Latent( AdmixOptions* op, Genome* loci)
{
  options = 0;
  rhoalpha = 0.0;
  rhobeta = 0.0;
  options = op;
  Loci = loci;
  poptheta = 0;
  globaltheta = 0;
  //globalthetaproposal = 0;
  rho.push_back(0.0);
  RhoSampler = 0;
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
      //globalthetaproposal = new double[K];
      fill(globaltheta, globaltheta+K, 1.0/(double)K);
      //ThetaTuner.SetParameters(1.0 /*<-initial stepsize on softmax scale*/, 0.00, 10.0, 0.44);
    }
    
    SumLogRho.push_back(0.0);
    // ** get prior on sum-of-intensities parameter rho or on rate parameter of its population distribution
    rhoalpha = options->getRhoalpha();
    if(options->getHapMixModelIndicator() || (options->getIndAdmixHierIndicator() && !options->isGlobalRho() )){
       // get prior on rate parameter beta and initialize it at prior mean
      rhobeta0 = options->getRhobetaShape();
      rhobeta1 = options->getRhobetaRate();
      rhobeta = rhobeta0 / rhobeta1;
      rho[0] = rhoalpha * rhobeta1 / (rhobeta0 - 1.0);
      if(options->getHapMixModelIndicator()){
	//initialise rho vector
	for(unsigned j = 0; j < Loci->GetNumberOfCompositeLoci()-1; ++j){
	  rho.push_back(rho[0]);
	  SumLogRho.push_back(0.0);
	}
	RhoArgs.NumPops = K;
	RhoArgs.NumLoci = Loci->GetNumberOfCompositeLoci();
	RhoArgs.Distances = Loci->GetDistances();
	RhoSampler = new HamiltonianMonteCarlo[Loci->GetNumberOfChromosomes()];
	for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c)
	  RhoSampler[c].SetDimensions(Loci->GetSizeOfChromosome(c)-1, 0.0001/*initial stepsize*/, 0.000/*min stepsize*/, 
				      1.0/*max stepsize*/, 10/*num leapfrogs*/,  0.8/*target accept rate*/, RhoEnergy, RhoGradient);
    
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
  //delete[] globalthetaproposal;
  delete[] RhoSampler;
}

// void Latent::UpdateGlobalTheta(int iteration, IndividualCollection* individuals){
//   if(options->getHapMixModelIndicator()){
//     if(!(iteration%2))ConjugateUpdateGlobalTheta(individuals->getSumLocusAncestry(options->getPopulations()));
//     else
//       UpdateGlobalThetaWithRandomWalk(individuals); 
//   }
//   individuals->setAdmixtureProps(globaltheta, options->getPopulations());//shouldn't be necessary
// }

void Latent::UpdatePopAdmixParams(int iteration, const IndividualCollection* const individuals, LogWriter &Log)
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
  
}

void Latent::UpdateGlobalSumIntensities(const IndividualCollection* const IC, bool sumlogrho) {
  if( options->isGlobalRho() ) { // update rho with random walk MH
    double rhoprop = rho[0];
    double LogLikelihood = 0.0;
    double LogLikelihoodAtProposal = 0.0;
    double LogLikelihoodRatio = 0.0;
    double LogPriorRatio = 0.0;
    double LogAccProbRatio = 0.0;
    bool accept = false;

    NumberOfUpdates++;
    rhoprop = exp(gennor(log(rho[0]), step)); // propose log rho from normal distribution with SD step
    
    //get log likelihood at current parameter values, annealed if this is an annealing run
    for(int i = 0; i < IC->getSize(); ++i) {
      Individual* ind = IC->getIndividual(i);
      ind->HMMIsBad(true);//to force HMM update
      LogLikelihood += ind->getLogLikelihood(options, false, true); // don't force update, store result if updated
      ind->HMMIsBad(true); // HMM probs overwritten by next indiv, but stored loglikelihood still ok
   }
     // set ancestry correlations using proposed value of sum-intensities
    // value for X chromosome set to half the autosomal value 
    Loci->SetLociCorr(rhoprop);

    //get log HMM likelihood at proposal rho and current admixture proportions
    for(int i = 0; i < IC->getSize(); ++i) {
      Individual* ind = IC->getIndividual(i);
      LogLikelihoodAtProposal += ind->getLogLikelihood(options, true, false); // force update, do not store result 
      ind->HMMIsBad(true); // set HMM probs as bad but stored log-likelihood is still ok
      // line above should not be needed for a forced update with result not stored
    }
    LogLikelihoodRatio = LogLikelihoodAtProposal - LogLikelihood;

    //compute prior ratio
    LogPriorRatio = getGammaLogDensity(rhoalpha, rhobeta, rhoprop) - getGammaLogDensity(rhoalpha, rhobeta, rho[0]);
    LogAccProbRatio = LogLikelihoodRatio + LogPriorRatio; 

    // generic Metropolis step
    if( LogAccProbRatio < 0 ) {
      if( log(myrand()) < LogAccProbRatio ) accept = true;
    } else accept = true;  
    
    if(accept) {
      rho[0] = rhoprop;
      for(int i = 0; i < IC->getSize(); ++i){
	Individual* ind = IC->getIndividual(i);
	ind->storeLogLikelihood(false); // store log-likelihoods calculated at rhoprop, but do not set HMM probs as OK 
      }
    } else { 
      // restore ancestry correlations in Chromosomes using original value of sum-intensities
      Loci->SetLociCorr(rho);
    } // stored loglikelihoods are still ok

    //update sampler object every w updates
    if( !( NumberOfUpdates % w ) ){
      step = TuneRhoSampler.UpdateStepSize( exp(LogAccProbRatio) );  
    }
    if(sumlogrho )SumLogRho[0] += log(rho[0]);// accumulate sum of log of sumintensities after burnin.
  }//end if global rho model

  else if(!options->getHapMixModelIndicator()){ //individual- or gamete-specific rho model
    if(IC->getSize()>1 && options->getIndAdmixHierIndicator() ) { // >1 individual and hierarchical model
      // update scale parameter of gamma distribution of sumintensities in population 
      if( options->isRandomMatingModel() )
	rhobeta = gengam( 2*rhoalpha * IC->getSize() + rhobeta0, IC->GetSumrho() + rhobeta1 );
      else
	rhobeta = gengam( rhoalpha* IC->getSize() + rhobeta0, IC->GetSumrho() + rhobeta1 );
    } // otherwise do not update rhobeta
    if(sumlogrho )SumLogRho[0] += log(rhoalpha) - log(rhobeta);// accumulate sum of log of mean of sumintensities after burnin.
  }
}

// void Latent::ConjugateUpdateGlobalTheta(const vector<int> sumLocusAncestry){
//   //sumlocusancestry is summed over all individuals
//   size_t K = options->getPopulations();
//   double dirparams[K];
//   for(size_t k = 0; k < K; ++k) {
//     dirparams[k] = alpha[0][k] + sumLocusAncestry[k] + sumLocusAncestry[k + K];
//   }
//   gendirichlet(K, dirparams, globaltheta );
// }

// void Latent::UpdateGlobalThetaWithRandomWalk(IndividualCollection* IC) {
//   double LogLikelihoodRatio = 0.0;
//   double LogPriorRatio = 0.0;
//   double logpratio = 0.0;
  
//   //generate proposals
//   // inverse softmax transformation from proportions to numbers on real line that sum to 0
//   bool* b = new bool[K];
//   double* a = new double[K]; // should be at class scope
//   for(int k = 0; k < K; ++k) {
//     if(globaltheta[k] > 0.0) {
//       b[k] = true; //to skip elements set to zero
//     } else b[k] = false;
//   }
//   inv_softmax(K, globaltheta, a, b);
//   //random walk step - on all elements of array a
//   for(int k = 0; k < K; ++k) {
//     if( b[k] ) a[k] = gennor(a[k], thetastep);  
//   }
//   //reverse transformation from numbers on real line to proportions 
//   softmax(K, globalthetaproposal, a, b);
//   //compute contribution of this gamete to log prior ratio
//   for(int k = 0; k < K; ++k) {
//     if( b[k] ) { 
//       // prior densities must be evaluated in softmax basis
//       LogPriorRatio += alpha[0][k]*( log(globalthetaproposal[k]) - log(globaltheta[k]) ); 
//     }
//   }
//   delete[] a;
//   delete[] b; 

//   vector<double> dummyrho(1);dummyrho[0] = rho[0];
//   for(int i = 0; i < IC->getSize(); ++i){
//     //get log likelihood at current parameter values - do not force update, store result of update
//     LogLikelihoodRatio -= IC->getIndividual(i)->getLogLikelihood(options, false, true); 
    
//     //get log likelihood at proposal theta and current rho - force update 
//     LogLikelihoodRatio += IC->getIndividual(i)->getLogLikelihood(options, globaltheta, globaltheta, dummyrho, dummyrho, true);
//   }
//   IC->HMMIsBad(true);

//   logpratio = LogLikelihoodRatio + LogPriorRatio;// log ratio of full conditionals
//   Accept_Reject_Theta(logpratio, K);
// }

// void Latent::Accept_Reject_Theta( double logpratio, int Populations) {
//   // Metropolis update for admixture proportions theta, taking log of acceptance probability ratio as argument
//   bool test = true;
//   bool accept = false;
//   double AccProb = 1.0; 
//   // loop over populations: if any element of proposed parameter vector is too small, reject update without test step
//   for( int k = 0; k < Populations; k++ ) {
//     if( globaltheta[ k ] > 0.0 && globalthetaproposal[ k ] < 0.0001 ) {
//       test = false;
//     }
//   }

//   if(test) { // generic Metropolis step
//     if( logpratio < 0 ) {
//       AccProb = exp(logpratio); 
//       if( myrand() < AccProb ) accept=true;
//     } else {
//       accept = true;
//     }
//   }
  
//   if(accept) { // set proposed values as new values    
//     copy(globalthetaproposal, globalthetaproposal+K, globaltheta);
//     //possibly need to reset loglikelihood in Individuals
//   } 
  
//   //update step size in tuner object every w updates
//   if( !( NumberOfUpdates % w ) ) {
//     thetastep = ThetaTuner.UpdateStepSize( AccProb );
//   }
// }

//conjugate update of locus-specific sumintensities, conditional on observed numbers of arrivals
void Latent::SampleSumIntensities(const vector<unsigned> &SumNumArrivals, unsigned NumIndividuals, bool sumlogrho) {
  double sum = 0.0;
  for(unsigned j = 1; j < rho.size(); ++j){
    double EffectiveL = Loci->GetDistance(j) * 2 * NumIndividuals;//length of interval * # gametes
    rho[j] = gengam( rhoalpha + (double)(SumNumArrivals[j]), rhobeta + EffectiveL );
    sum += rho[j];
  }

  //sample rate parameter of gamma prior on rho
  rhobeta = gengam( rhoalpha * (double)(rho.size()-1) + rhobeta0, sum + rhobeta1 );
  //set locus correlation
  Loci->SetLociCorr(rho);

  //accumulate sums of log of rho
  if(sumlogrho)
    transform(rho.begin(), rho.end(), SumLogRho.begin(), SumLogRho.begin(), std::plus<double>());
}

void Latent::SampleSumIntensities(const int* SumAncestry, bool sumlogrho){
  RhoArgs.SumAncestry = SumAncestry;
  RhoArgs.theta = globaltheta;
  vector<double>::iterator p = rho.begin()+1;
  double sum = 0.0;
  for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
    Chromosome* C = Loci->getChromosome(c); 
    //take logs of rho
    for(unsigned i = 0; i < C->GetSize()-1; ++i)
      *(p + i) = log(*(p +i));
    RhoArgs.NumLoci = C->GetSize();
    RhoArgs.Distances = C->GetDistances();

    RhoSampler[c].Sample(&(*p), &RhoArgs);// *p is the element in rho corresponding to second locus on current chromosome
    //take exponents of logrho
    for(unsigned i = 0; i < C->GetSize()-1; ++i){
      *(p + i) = exp(*(p + i));
      sum += *(p+i);
    }
    p += C->GetSize();
    RhoArgs.SumAncestry += C->GetSize()*(options->getPopulations()+1);
   }
  //sample rate parameter of gamma prior on rho
  rhobeta = gengam( rhoalpha * (double)(rho.size()-1) + rhobeta0, sum + rhobeta1 );
  //set locus correlation
  Loci->SetLociCorr(rho);

  //accumulate sums of log of rho
  if(sumlogrho)
    for(vector<double>::iterator i = SumLogRho.begin(), p=rho.begin(); i < SumLogRho.end(); ++i, ++p)*i += log(*p); 
}

double Latent::RhoEnergy(const double* const x, const void* const vargs){
  //x is the log of rho, with length NumLoci - 1
  const RhoArguments* args = (const RhoArguments*)vargs;
  unsigned K = args->NumPops;
  unsigned L = args->NumLoci;
  const int* n = args->SumAncestry;
  const double* d = args->Distances;
  const double* theta = args->theta;
  double rho, f;
  double E = 0.0;
  gsl_sf_result result;
  int status = 0;

  for(unsigned j = 1; j < L; ++j){
    rho = exp(x[j-1]);
    status = gsl_sf_exp_e(d[j]*rho, &result);
    if(status)throw("exp error in RhoEnergy");
    f = result.val;
    status = gsl_sf_lngamma_e(1.0-f, &result);
    if(status)throw("log error in RhoEnergy");
    E += n[j*(K+1)] * result.val;
    for(unsigned k = 0; k < K; ++k)
      E += n[j*(K+1) + k+1] * log(f + theta[k]*(1.0 - f));
  }
  return -E; 
}

void Latent::RhoGradient( const double* const x, const void* const vargs, double* g ){
  const RhoArguments* args = (const RhoArguments*)vargs;
  unsigned K = args->NumPops;
  unsigned L = args->NumLoci;
  const int* n = args->SumAncestry;
  const double* d = args->Distances;
  const double* theta = args->theta;
  double rho, f;
  gsl_sf_result result;
  int status = 0;

  for(unsigned j = 1; j < L; ++j){
    rho = exp(x[j-1]);
    status = gsl_sf_exp_e(d[j]*rho, &result);
    if(status)throw("exp error in RhoGradient");
    f = result.val;
    g[j-1] = n[j*(K+1)] / (1.0 - f);
    for(unsigned k = 0; k < K; ++k)
      g[j-1] -= (n[j*(K+1) + k+1]*(1.0 - theta[k])) / (f + theta[k]*(1.0 - f));
    g[j-1] *= rho*d[j]*f;//jacobian
  }
}

void Latent::InitializeOutputFile(const std::string* const PopulationLabels)
{
  // Header line of paramfile
  //Pop. Admixture
  if(!options->getHapMixModelIndicator())
    for( int i = 0; i < options->getPopulations(); i++ ) {
      outputstream << "\""<<PopulationLabels[i] << "\"\t";
    }
  //SumIntensities
  if(options->getHapMixModelIndicator())
    outputstream << "SumIntensities.Mean\tSumIntensities.Variance";
  else{
    if( options->isGlobalRho() ) outputstream << "sumIntensities\t";
    else outputstream << "sumIntensities.mean\t";
  }
  outputstream << endl;
  
}

void Latent::OutputErgodicAvg( int samples, std::ofstream *avgstream)
{
  if(options->getPopulations()>1){
    if(!options->getHapMixModelIndicator())
      for( int j = 0; j < options->getPopulations(); j++ ){
	avgstream->width(9);
	*avgstream << setprecision(6) << SumAlpha[j] / samples << "\t";
      }
    avgstream->width(9);
    *avgstream << setprecision(6) << exp(SumLogRho[0] / samples) << "\t";
  }
}

//output to given output stream
void Latent::OutputParams(ostream* out){
  //pop admixture params
  if(!options->getHapMixModelIndicator())
    for( int j = 0; j < options->getPopulations(); j++ ){
      out->width(9);
      (*out) << setprecision(6) << alpha[0][ j ] << "\t";
    }
  //sumintensities
  out->width(9);
  if(options->getHapMixModelIndicator()){
    double sum = accumulate(rho.begin()+1, rho.end(), 0.0, std::plus<double>());
    double var = 0.0;
    for(unsigned j = 1; j < rho.size(); ++j)var += rho[j]*rho[j];
    var = var - (sum*sum) / (double)rho.size();
    (*out) << setiosflags(ios::fixed) << setprecision(6) << sum / rho.size() << "\t" << var /(rho.size()-1) << "\t";
  }
  else{
    if( options->isGlobalRho() )
      (*out) << setprecision(6) << rho[0] << "\t";
    else
      (*out) << setprecision(6) << rhoalpha / rhobeta  << "\t";
  }
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
const vector<double> &Latent::getSumLogRho()const{
  return SumLogRho;
}
const double *Latent::getpoptheta()const{
  return poptheta;
}

void Latent::printAcceptanceRates(LogWriter &Log) {
  if(!options->getHapMixModelIndicator()){
//     Log << "Expected acceptance rate in global admixture sampler: "
// 	<< ThetaTuner.getExpectedAcceptanceRate()
// 	<< "\nwith final step size of "
// 	<< thetastep << "\n";
//   }
//   else{
    Log << "Expected acceptance rate in population admixture sampler: "
	<< PopAdmixSampler.getExpectedAcceptanceRate()
	<< "\nwith final step size of "
	<< PopAdmixSampler.getStepSize() << "\n";
  }
}
double Latent::getRhoSamplerAccRate()const{
  if(options->getHapMixModelIndicator())
    return RhoSampler[0].getAcceptanceRate();
  else
    return TuneRhoSampler.getExpectedAcceptanceRate();
}

double Latent::getRhoSamplerStepsize()const{
  if(options->getHapMixModelIndicator())
    return RhoSampler[0].getStepsize();
  else
    return step;
}
