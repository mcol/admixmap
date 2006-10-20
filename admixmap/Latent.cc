/** 
 *   ADMIXMAP
 *   Latent.cc 
 *   Class to hold and update population admixture and sumintensities parameters and their priors
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 */
#include "Latent.h"
#include "Chromosome.h"
#include "utils/misc.h"
#include "utils/dist.h"//for log gamma density
#include <algorithm>
#include <numeric>
#include "gsl/gsl_math.h"
#include "gsl/gsl_specfunc.h"
#include "Comms.h"
#ifdef PARALLEL
#include <mpe.h>//for MPI event logging
#endif
using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

Latent::Latent( AdmixOptions* op, Genome* loci)
{
  options = 0;
  //rhopriormean = 0.0;
  //rhopriorvar = 0.0;
  options = op;
  Loci = loci;
  poptheta = 0;
  globaltheta = 0;
  //globalthetaproposal = 0;
  rho.push_back(0.0);
  HapMixLambdaSampler = 0;

}

void Latent::Initialise(int Numindividuals, const Vector_s& PopulationLabels, LogWriter &Log){
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
    if(!options->getHapMixModelIndicator())PopAdmixSampler.SetSize( obs, K );
    
    //initialise global admixture proportions
    if(options->getHapMixModelIndicator()){
      globaltheta = new double[K];
      //globalthetaproposal = new double[K];
      fill(globaltheta, globaltheta+K, 1.0/(double)K);
      //ThetaTuner.SetParameters(1.0 /*<-initial stepsize on softmax scale*/, 0.00, 10.0, 0.44);
    }

    if(Comms::isMaster())SumLogRho.push_back(0.0);
    // ** get prior on sum-of-intensities parameter rho or on rate parameter of its population distribution

    if(options->getHapMixModelIndicator()){
      //set prior of rho prior rate
	const vector<double>& priorparams = options->getHapMixLambdaPrior();

	LambdaArgs.h = rhoalpha = priorparams[0]; // h
	LambdaArgs.beta_shape = priorparams[1];
	LambdaArgs.beta_rate = priorparams[2];
	LambdaArgs.beta = rhobeta = priorparams[1]/priorparams[2];//initialise beta at prior mean

      unsigned numIntervals = Loci->GetNumberOfCompositeLoci()-Loci->GetNumberOfChromosomes();
      if(Comms::isMaster()){
	LambdaArgs.NumPops = K;
	LambdaArgs.NumIntervals = numIntervals;
	
	//set up Hamiltonian sampler for rho
 	HapMixLambdaSampler = new HamiltonianMonteCarlo[numIntervals];
	const vector<float>& lambdasamplerparams = options->getrhoSamplerParams();
	size_t size = lambdasamplerparams.size();
	float initial_stepsize = size? lambdasamplerparams[0] : 0.06;
	float min_stepsize = size? lambdasamplerparams[1] : 0.0001;
	float max_stepsize = size? lambdasamplerparams[2] : 1.0;
	float target_acceptrate = size? lambdasamplerparams[3] : 0.9;
	int num_leapfrog_steps = size? (int)lambdasamplerparams[4] : 20;
	NumberOfUpdates = 0;
	w = 1;
	
//Hamiltonian sampler
  	for(unsigned j = 0; j < numIntervals; ++j){
  	  HapMixLambdaSampler[j].SetDimensions(1, initial_stepsize, min_stepsize, max_stepsize, num_leapfrog_steps, 
  				      target_acceptrate, LambdaEnergy, LambdaGradient);
 	}
	
      }//end sampler initialisation
      //initialise rho vector
      int locus = 0;
      const double initial_lambda = options->getInitialHapMixLambda();
      rho.clear();
      for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
	  ++locus;//skip first locus on each chromosome
	  for(unsigned i = 1; i < Loci->GetSizeOfChromosome(c); ++i){
	      if(initial_lambda >0.0)
		  rho.push_back(initial_lambda);
	      else
		  rho.push_back(Rand::gengam(rhoalpha*Loci->GetDistance(locus), rhobeta));
	      if(Comms::isMaster()){
		  SumLogRho.push_back(0.0);
	      }
	      ++locus;
	  }
      }

      rhoproposal.resize(numIntervals);
    }//end if hapmixmodel
    else{
      rhoalpha = options->getRhoalpha();
      if( (options->getIndAdmixHierIndicator() && !options->isGlobalRho() )){
	// get prior on rate parameter beta and initialize it at prior mean
	rhobeta0 = options->getRhobetaShape();
	rhobeta1 = options->getRhobetaRate();
	rhobeta = rhobeta0 / rhobeta1;
	double initial_rho = rhoalpha * rhobeta1 / (rhobeta0 - 1.0);
	rho[0] = initial_rho;
	
      }
      else{//global rho or non-hierarchical model on rho
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
	const vector<float>& rhosamplerparams = options->getrhoSamplerParams();
	size_t size = rhosamplerparams.size();
	float initial_stepsize = size? rhosamplerparams[0] : step0;
	float min_stepsize = size? rhosamplerparams[1] : 0.01;
	float max_stepsize = size? rhosamplerparams[2] : 10;
	float target_acceptrate = size? rhosamplerparams[3] : 0.44;
	TuneRhoSampler.SetParameters( initial_stepsize, min_stepsize, max_stepsize, target_acceptrate);
	}
      }
    }

    // ** Open paramfile **
    if ( Comms::isMaster() && options->getIndAdmixHierIndicator()){
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

void Latent::resetStepSizeApproximator(int k) {
    TuneRhoSampler.resetStepSizeApproximator(k);
}


Latent::~Latent()
{
  delete[] poptheta;
  delete[] globaltheta;
  //delete[] globalthetaproposal;
  if(options->getHapMixModelIndicator()) {
    delete[] HapMixLambdaSampler;
  }
}

/** Samples for population admixture distribution Dirichlet parameters, alpha **
 For a model in which the distribution of individual admixture in the population is a mixture
 of components, we have one Dirichlet parameter vector for each component, 
 updated only from those individuals who belong to the component
*/
void Latent::UpdatePopAdmixParams(int iteration, const IndividualCollection* const individuals, LogWriter &Log)
 {
   if( options->getPopulations() > 1 && individuals->getSize() > 1 &&
       options->getIndAdmixHierIndicator() ){
     const double* sumlogtheta = individuals->getSumLogTheta();

     if(Comms::isMaster()){
       //sample alpha conditional on individual admixture proportions
       //cout << "alpha " << alpha[0][0] << " " << alpha[0][1] <<  " sumlogtheta " 
       //	    << individuals->getSumLogTheta()[0] << " " <<  individuals->getSumLogTheta()[1] << endl;
       try{
	 PopAdmixSampler.Sample( sumlogtheta, &alpha[0], options->PopAdmixturePropsAreEqual() );
       }
       catch(string s){
	 throw string("Error encountered while sampling population admixture parameters:\n" +s);
       }
     }
#ifdef PARALLEL
     Comms::BroadcastVector(alpha[0]);
#endif
     copy(alpha[0].begin(), alpha[0].end(), alpha[1].begin()); // alpha[1] = alpha[0]

  }
   if(Comms::isMaster()){
     // ** accumulate sum of Dirichlet parameter vector over iterations  **
     transform(alpha[0].begin(), alpha[0].end(), SumAlpha.begin(), SumAlpha.begin(), std::plus<double>());//SumAlpha += alpha[0];
     
     if( iteration == options->getBurnIn() && options->getPopulations() > 1) {
       if(options->getNumberOfOutcomes() > 0){
	 Log << Off << "Individual admixture centred in regression model around: ";
	 for(int i = 0; i < options->getPopulations(); ++i)Log << poptheta[i] << "\t";
	 Log << "\n";
       }
       fill(SumAlpha.begin(), SumAlpha.end(), 0.0);
     }
   }
   
   if( iteration < options->getBurnIn() && options->getPopulations() > 1) {
     if(Comms::isMaster()){
       // accumulate ergodic average of population admixture, which is used to centre 
       // the values of individual admixture in the regression model
       double sum = accumulate(SumAlpha.begin(), SumAlpha.end(), 0.0);
       if(options->getNumberOfOutcomes() > 0)for( int j = 0; j < options->getPopulations(); j++ )poptheta[j] = SumAlpha[j] / sum;
     }
#ifdef PARALLEL
     Comms::Broadcast(poptheta, options->getPopulations());
#endif
   }
   
 
}

///updates global sumintensities in a globalrho model, using random-walk Metropolis-Hastings
void Latent::UpdateGlobalSumIntensities(const IndividualCollection* const IC, bool sumlogrho) {
  if( options->isGlobalRho() ) {
    double LogLikelihood = 0.0;
    double LogLikelihoodAtProposal = 0.0;
    double LogLikelihoodRatio = 0.0;
    int accept = 0;
    double rhoprop = rho[0];
    double logrhoprop = 0.0;
    double logrho0 = log(rho[0]);
    const int NumWorkers = Comms::getNumWorkers();

    if(Comms::isMaster()){
     
      NumberOfUpdates++;
      logrhoprop = Rand::gennor(logrho0, step);
      rhoprop = exp(logrhoprop); // propose log rho from normal distribution with SD step
    }
#ifdef PARALLEL
    //broadcast proposal
    Comms::Broadcast(&rhoprop);
#endif
    
    if(Comms::isWorker()){
       //get log likelihood at current parameter values, annealed if this is an annealing run
      for(int i = Comms::getWorkerRank(); i < IC->getSize(); i += NumWorkers) {
	Individual* ind = IC->getIndividual(i);
	ind->HMMIsBad(true);//to force HMM update
	LogLikelihood += ind->getLogLikelihood(options, false, true); // don't force update, store result if updated
	ind->HMMIsBad(true); // HMM probs overwritten by next indiv, but stored loglikelihood still ok
      }
      // set ancestry correlations using proposed value of sum-intensities
      // value for X chromosome set to half the autosomal value 
      Loci->SetLocusCorrelation(rhoprop);
      
      //get log HMM likelihood at proposal rho and current admixture proportions
      for(int i = Comms::getWorkerRank(); i < IC->getSize(); i += NumWorkers) {
	Individual* ind = IC->getIndividual(i);
	LogLikelihoodAtProposal += ind->getLogLikelihood(options, true, false); // force update, do not store result 
	ind->HMMIsBad(true); // set HMM probs as bad but stored log-likelihood is still ok
	// line above should not be needed for a forced update with result not stored
      }
      LogLikelihoodRatio = LogLikelihoodAtProposal - LogLikelihood;
    }//end worker block
  
#ifdef PARALLEL
    //reduce logL ratio
    Comms::Reduce(&LogLikelihoodRatio);
#endif
    
    if(Comms::isMaster()){
      //compute log ratio of prior densities in log rho basis
      double LogPriorRatio = rhoalpha * (logrhoprop - logrho0) - rhobeta * (rhoprop - rho[0]); 
      // getGammaLogDensity(rhoalpha, rhobeta, rhoprop) - getGammaLogDensity(rhoalpha, rhobeta, rho[0]);
      double LogAccProbRatio = LogLikelihoodRatio + LogPriorRatio; 
      
      // generic Metropolis step
      if( LogAccProbRatio < 0 ) {
	if( log(Rand::myrand()) < LogAccProbRatio ) accept = true;
      } else accept = 1;  

      if(accept){
	rho[0] = rhoprop;
	logrho0 = logrhoprop;
      }
      //update sampler object every w updates
      if( !( NumberOfUpdates % w ) ){
	step = TuneRhoSampler.UpdateStepSize( exp(LogAccProbRatio) );  
      }
      if(sumlogrho )SumLogRho[0] += logrho0;// accumulate sum of log of sumintensities after burnin.
    }
#ifdef PARALLEL
    //tell workers whether proposal has been accepted
    MPE_Log_event(17, 0, "Bcastrho");
    Comms::Broadcast(&accept);
    MPE_Log_event(18, 0, "Bcasted");
#endif

    if(Comms::isWorker()){      
      if(accept) {
	for(int i = Comms::getWorkerRank(); i < IC->getSize(); i += NumWorkers){
	  Individual* ind = IC->getIndividual(i);
	  ind->storeLogLikelihood(false); // store log-likelihoods calculated at rhoprop, but do not set HMM probs as OK 
	}
      } else { 
	// restore ancestry correlations in Chromosomes using original value of sum-intensities
	Loci->SetLocusCorrelation(rho);
      } // stored loglikelihoods are still ok
    }
  }//end if global rho model
  
  else if(!options->getHapMixModelIndicator()){ //individual- or gamete-specific rho model
    if(IC->getSize()>1 && options->getIndAdmixHierIndicator() ) { // >1 individual and hierarchical model
      double sumrho = IC->GetSumrho();    
      if(Comms::isMaster()){
	
	// update scale parameter of gamma distribution of sumintensities in population 
	if( options->isRandomMatingModel() )
	  rhobeta = Rand::gengam( 2*rhoalpha * IC->getSize() + rhobeta0, sumrho + rhobeta1 );
	else
	  rhobeta = Rand::gengam( rhoalpha* IC->getSize() + rhobeta0, sumrho + rhobeta1 );
      }
#ifdef PARALLEL
    Comms::Broadcast(&rhobeta);
#endif
    } // otherwise do not update rhobeta

    // accumulate sum of log of mean of sumintensities after burnin.
    if(sumlogrho && Comms::isMaster())SumLogRho[0] += log(rhoalpha) - log(rhobeta);
  }
}

///samples locus-specific number of arrivals in a hapmixmodel, using Hamiltonian Monte Carlo
void Latent::SampleHapMixLambda(const int* SumAncestry, bool accumulateLogs){
  //SumAncestry is a (NumberOfCompositeLoci) * 2 array of counts of gametes with ancestry states that are
  // unequal (element 0) and equal 
  //sumlogrho indicates whether to accumulate sums of log rho

  //LambdaArgs.theta = globaltheta;
  //double sum = 0.0;
  int locus = 0;
  int index = 0;
  
  if(Comms::isMaster()){
#ifdef PARALLEL
    MPE_Log_event(9, 0, "sampleLambda");
#endif
    LambdaPriorArgs.sumlambda = 0.0;
    try{
      vector<double>::iterator lambda_iter = rho.begin();
      vector<double>::iterator sumloglambda_iter = SumLogRho.begin();
      for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
	++locus;//skip first locus on each chromosome
	for(unsigned i = 1; i < Loci->GetSizeOfChromosome(c); ++i){
	  double loglambda = log(*lambda_iter);//sampler is on log scale
	  
	  LambdaArgs.Distance = Loci->GetDistance(locus);//distance between this locus and last
	  LambdaArgs.SumAncestry = SumAncestry + locus*2;//sums of ancestry for this locus
	  
	  //sample new value
	  HapMixLambdaSampler[index++].Sample(&loglambda, &LambdaArgs);
	  *lambda_iter = exp(loglambda);

	  //accumulate sums of log of rho
	  if(accumulateLogs)
	    *(sumloglambda_iter++) += loglambda;
	  //accumulate sums, used to sample beta
	  LambdaPriorArgs.sumlambda += *lambda_iter;

	  ++locus;
	  ++lambda_iter;
	}
      }
    }
    catch(string s){
      stringstream err;err << "Error encountered while sampling lambda " << locus << ":\n" + s;
      throw(err.str());
    }
#ifdef PARALLEL
    MPE_Log_event(10, 0, "sampledLambda");
#endif
  }
//broadcast rho to all processes, in Comm (excludes freqsampler)
//no need to broadcast globaltheta if it is kept fixed
#ifdef PARALLEL
  MPE_Log_event(17, 0, "Bcastlambda");
  Comms::BroadcastVector(rho);
  MPE_Log_event(18, 0, "Bcasted");
#endif    
  //set locus correlation (workers only)
  if(Comms::isWorker())  
    {
      Loci->SetLocusCorrelation(rho);
      for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
	//set global state arrival probs in hapmixmodel
	//TODO: can skip this if xonly analysis with no females
	//NB: assumes always diploid in hapmixmodel
	Loci->getChromosome(c)->SetStateArrivalProbs(options->isRandomMatingModel());
      }
    }
  if(Comms::isMaster()){
      SampleHapmixLambdaPriorParameters();
  }
}

void Latent::SampleHapmixLambdaPriorParameters(){
   try{
       const double pshape = LambdaArgs.beta_shape + LambdaArgs.h*Loci->GetLengthOfGenome();
       const double prate = LambdaArgs.beta_rate + LambdaPriorArgs.sumlambda;
       rhobeta = Rand::gengam(pshape, prate);
       LambdaArgs.beta = rhobeta;

   }
   catch(string s){
     string err = "Error encountered while sampling lambda prior params:\n" + s;
     throw(err);
   }
  
}

///energy function for sampling locus-specific lambda, expected #arrivals
///conditional on locus ancestry and gamma prior params 
double Latent::LambdaEnergy(const double* const x, const void* const vargs){
  //x is log rho
  const LambdaArguments* args = (const LambdaArguments*)vargs;
  unsigned K = args->NumPops;
  const double d = args->Distance;
  double theta = 1.0 / (double)K;
  double E = 0.0;
  try {
    double rho = myexp(*x);
    double f = myexp(-rho);
    int sumequal = args->SumAncestry[1], sumnotequal = args->SumAncestry[0];
    double probequal = theta + f*(1.0 - theta);
    E -= sumnotequal * log(1.0-probequal);//constant term in log(1-theta) omitted
    E -= sumequal * log(probequal);

    // log unnormalized gamma prior on log rho
    E -= args->h*d * (*x) - args->beta * rho;
  } catch(string s){
    throw string("Error in LambdaEnergy: " + s);
  }
  return E; 
}

///gradient function for sampling locus-specific lambda, expected #arrivals
///conditional on locus ancestry and gamma prior params
void Latent::LambdaGradient( const double* const x, const void* const vargs, double* g ){
  const LambdaArguments* args = (const LambdaArguments*)vargs;
  unsigned K = args->NumPops;
  const double d = args->Distance;
  double theta = 1.0 / (double)K;
  g[0] = 0.0;
  try {
    double rho = myexp(*x);
    double f = myexp(-rho);
    //first compute dE / df
    int sumequal = args->SumAncestry[1], sumnotequal = args->SumAncestry[0];
    double probequal = theta + f*(1.0 - theta);

    g[0] +=  sumnotequal / (1.0 - probequal) - sumequal / probequal;//dprobequal / df
    g[0] *= -rho*d*f*(1.0-theta); // dprobequal / df * df / dx 
    
    //derivative of minus log prior wrt log rho
    g[0] -= args->h*d  - args->beta * rho;

  } catch(string s) {
    throw string("Error in LambdaGradient: " + s);
  }
}

void Latent::InitializeOutputFile(const Vector_s& PopulationLabels) {
  if(Comms::isMaster()){
    // Header line of paramfile
    //Pop. Admixture
    if(!options->getHapMixModelIndicator())
      for( int i = 0; i < options->getPopulations(); i++ ) {
	outputstream << "\""<<PopulationLabels[i] << "\"\t";
      }
    //SumIntensities
    if(options->getHapMixModelIndicator())
      outputstream << "LambdaIntensities.Mean\tLambda.Variance\th\tbeta";
    else{
      if( options->isGlobalRho() ) outputstream << "sumIntensities\t";
      else outputstream << "sumIntensities.mean\t";
    }
    outputstream << endl;
  }
}

void Latent::OutputErgodicAvg( int samples, std::ofstream *avgstream)
{
  if(!options->getHapMixModelIndicator()){
    if(options->getPopulations()>1){
      for( int j = 0; j < options->getPopulations(); j++ ){
	avgstream->width(9);
	*avgstream << setprecision(6) << SumAlpha[j] / samples << "\t";
      }
      avgstream->width(9);
      *avgstream << setprecision(6) << exp(SumLogRho[0] / samples) << "\t";
    }
  }
  else{
    double sum = 0.0;
    for(vector<double>::const_iterator i = SumLogRho.begin(); i < SumLogRho.end(); ++i)sum += exp(*i / samples);
    avgstream->width(9);
    *avgstream << setprecision(6) << sum / (Loci->GetNumberOfCompositeLoci() - Loci->GetNumberOfChromosomes()) << "\t";
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
    double sum = 0.0;
    double var = 0.0;

    //unsigned j = 1;
    for(vector<double>::const_iterator lambda_iter = rho.begin(); lambda_iter != rho.end(); ++lambda_iter){
	double r = *lambda_iter;// / Loci->GetDistance(j++);
	sum += r;
	var += r * r;
    }
    double size = (double)rho.size();
    var = var - (sum*sum) / size;
    (*out) << setiosflags(ios::fixed) << setprecision(6) << sum / size << "\t" << var /size << "\t"
	   << LambdaArgs.h << "\t" << LambdaArgs.beta << "\t" ;
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
const double *Latent::getGlobalTheta()const{
  return globaltheta;
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
  if(options->getHapMixModelIndicator()){
    double av = 0;//average acceptance rate
    for(unsigned j = 0; j < rho.size(); ++j){
	av += HapMixLambdaSampler[j].getAcceptanceRate();
    }
    Log << "Average Expected acceptance rate in sumintensities samplers: "
	<< av / (double)(Loci->GetNumberOfCompositeLoci()-Loci->GetNumberOfChromosomes());

    av = 0;//average stepsize
    for(unsigned j = 0; j < rho.size(); ++j){
      av += HapMixLambdaSampler[j].getStepsize();
      if(isnan(HapMixLambdaSampler[j].getStepsize()))cout << j << " " << HapMixLambdaSampler[j].getStepsize() << " " << HapMixLambdaSampler[j].getAcceptanceRate() << endl;  
    }
    Log << "\nwith average final step size of "
	<< av / (double)(Loci->GetNumberOfCompositeLoci()-Loci->GetNumberOfChromosomes())
	<< "\n";

  }
  else if( options->isGlobalRho() ){
    Log << "Expected acceptance rate in global sumintensities sampler: "
	<< TuneRhoSampler.getExpectedAcceptanceRate()
	<< "\nwith final step size of "
	<< step	<< "\n";
  }

}
