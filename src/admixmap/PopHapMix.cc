/** 
 *   HAPMIXMAP
 *   PopHapMix.cc 
 *   Class to hold and update population parameters in a hapmixmodel
 *   Copyright (c) 2006 - 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 */
#include "PopHapMix.h"
#include "Chromosome.h"
#include "utils/misc.h"
#include "utils/dist.h"//for log gamma density
#include <algorithm>
#include <numeric>
#include <sstream>
#include "gsl/gsl_math.h"
#include "gsl/gsl_specfunc.h"
#include "Comms.h"
#include "MixturePropsWrapper.hh"
#include "EventLogger.hh"//for MPE event logging
using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

PopHapMix::PopHapMix( HapMixOptions* op, Genome* loci)
{
  options = op;
  Loci = loci;
  MixtureProps = 0;
  //MixturePropsproposal = 0;
  HapMixLambdaSampler = 0;
  MixturePropsPrior = 0;
  dirparams = 0;
}

void PopHapMix::Initialise(const string& distanceUnit, LogWriter& Log, const vector<string>& BlockStateLabels){
  Log.setDisplayMode(On);
  K = options->getNumberOfBlockStates();

  if(K > 1){ 
    if ( Comms::isMaster()){
      
      InitialiseMixtureProportions(Log);
      
      InitialiseArrivalRates(Log);
      
      // ** Open paramfile **
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
          InitializeOutputFile(BlockStateLabels, distanceUnit);
        }
      }
      else{
        Log << "No paramfile given\n";
      }
    }
    
    //broadcast lambda to all processes, in Comm (excludes freqsampler)
    //no need to broadcast MixtureProps if it is kept fixed
#ifdef PARALLEL
    EventLogger::LogEvent(17, 0, "Bcastlambda");
    Comms::BroadcastVector(lambda);
    EventLogger::LogEvent(18, 0, "Bcasted");
    //TODO: broadcast Mixture Params
#endif    

  }//end if K > 1
}

///initialise global admixture proportions
void PopHapMix::InitialiseMixtureProportions(LogWriter& Log){
  const unsigned L = Loci->GetNumberOfCompositeLoci();

  MixtureProps = new double[K*L];
  MixturePropsPrior = new double[K];
  dirparams = new double[K];
  
  /*
    set prior on mixture props.
    Mixture proportions have a Dirichlet prior with parameters MixturePropsDispersion / K 
    ( so that the proportions are equal). The prior dispersion (MixturePropsPriorDispersion) has a 
    Gamma prior, which defaults to Gamma(1,1).
  */
  
  //defaults 
  double MPDShape = (double)K;// dispersion prior shape
  double MPDRate = 1.0;//dispersion prior rate
  double DirParamInit = 1.0;//initial value of Dirichlet Params
  //check for user-specified prior
  const vector<double>& userprior = options->getMixturePropsPrior();
  if(userprior.size()){//user has specified a prior
    MPDShape = userprior[0];
    MPDRate = userprior[1];
    DirParamInit = (MPDShape / MPDRate) / (double)K;
  }
  
  //set initial values of Dirichlet prior params
  fill(MixturePropsPrior, MixturePropsPrior + K, DirParamInit);
  
  //write prior parameters to screen and log
  Log << Quiet << "Dirichlet prior on mixture proportions with parameters: " << MixturePropsPrior[0];
  for(unsigned k = 1; k < K; ++k) Log << ", " << MixturePropsPrior[k];
  Log << "\n";
  Log << "Gamma(" << MPDShape << ", " << MPDRate << ") prior on mixture proportion prior dispersion\n" ;
  
  //set initial values of Mixture Props as uniform
  //TODO: initialise to prior means
  fill(MixtureProps, MixtureProps + K*L, 1.0/(double)K);
  
  //setup sampler for mixure props dispersion
  MixturePropsDispersionSampler.SetSize(L, K, 0.01/*<-initial stepsize*/, 20/*<-num leapfrog steps*/);
  MixturePropsDispersionSampler.SetPriorEta(MPDShape, MPDRate);
  
  //MixturePropsproposal = new double[K];
  //ThetaTuner.SetParameters(1.0 /*<-initial stepsize on softmax scale*/, 0.00, 10.0, 0.44);
  
}

void PopHapMix::InitialiseArrivalRates(LogWriter& Log){
 const unsigned L = Loci->GetNumberOfCompositeLoci();
  //if(Comms::isMaster())SumLogLambda.push_back(0.0);
  //set priors
  const vector<double>& priorparams = options->getHapMixLambdaPrior();
    
  LambdaArgs.beta_shape = priorparams[2];
  LambdaArgs.beta_rate = priorparams[3];
  LambdaArgs.beta = priorparams[2]/priorparams[3];//initialise beta at prior mean
  unsigned numIntervals = L-Loci->GetNumberOfChromosomes(); 
  lambda.assign(numIntervals, 0);
  //TODO: if distance is too large, leave lambda at zero and do not sample (as if a new chromosome)
   
  LambdaArgs.NumBlockStates = K;
  LambdaArgs.NumIntervals = numIntervals;
  
  //set up Hamiltonian sampler for lambda
  HapMixLambdaSampler = new HamiltonianMonteCarlo[numIntervals];
  const vector<float>& lambdasamplerparams = options->getLambdaSamplerParams();
  size_t size = lambdasamplerparams.size();
  float initial_stepsize = size? lambdasamplerparams[0] : 0.06;
  float min_stepsize = size? lambdasamplerparams[1] : 0.0001;
  float max_stepsize = size? lambdasamplerparams[2] : 1.0;
  float target_acceptrate = size? lambdasamplerparams[3] : 0.9;
  int num_leapfrog_steps = size? (int)lambdasamplerparams[4] : 20;
  
  //set constant args for sampler for h
  hargs.shape = priorparams[0];
  hargs.rate = priorparams[1];
  hargs.Dlogbeta = Loci->GetLengthOfGenome()*log(LambdaArgs.beta);
  hargs.sum_lngamma_hd = 0.0;
  hargs.NumIntervals = numIntervals;
  hargs.distances = new double[numIntervals];
  LambdaArgs.h = hargs.shape / hargs.rate; // initialise h at prior mean 
  
  //write prior to screen and log
  Log << Quiet << "Gamma(h, beta) distribution on number of arrivals per Mb. \nGamma( " << hargs.shape << ", " << hargs.rate << " ) prior on h and "
      << "Gamma( " << LambdaArgs.beta_shape << ", " << LambdaArgs.beta_rate << " ) prior on beta\n";
  
  //Hamiltonian sampler
  for(unsigned j = 0; j < numIntervals; ++j){
    HapMixLambdaSampler[j].SetDimensions(1, initial_stepsize, min_stepsize, max_stepsize, num_leapfrog_steps, 
                                         target_acceptrate, LambdaEnergy, LambdaGradient);
  }
  //Random Walk sampler for log h
  //hTuner.SetParameters( 0.01, 0.0001, 1000.0, 0.44);
  
  //Adaptive rejection sampler for h
  hARS.Initialise(false, true, 1000000.0/*<- upper bound*/, 10.0/*<-lower bound*/, hlogf, hdlogf);
  
  //end sampler initialisation
  //initialise lambda vector
  int locus = 0;//indexes loci
  int d = 0; //indexes intervals
  const char* initfilename = options->getInitialHapMixLambdaFilename();
  const bool useinitfile = (strlen(initfilename) > 0);
  ifstream initfile;
  if(useinitfile){
    Log << Quiet << "Reading initial values of arrival rates from " << initfilename << "\n";
    initfile.open(initfilename);
    //read initial values of h, beta
    initfile >> LambdaArgs.h >> LambdaArgs.beta;
    //copy initial values of lambda straight from file into vector
    istream_iterator<double>firstlambda(initfile);
    //copy(firstlambda, firstlambda + numIntervals, lambda.begin());
    copy(firstlambda, istream_iterator<double>(), lambda.begin());//does not check number of elements
    initfile.close();
  }
  
  for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
    ++locus;//skip first locus on each chromosome
    for(unsigned i = 1; i < Loci->GetSizeOfChromosome(c); ++i){
      hargs.sum_lngamma_hd += lngamma(LambdaArgs.h*Loci->GetDistance(locus));
      hargs.distances[d] = Loci->GetDistance(locus);
      if(!useinitfile)//initialise at prior mean
        lambda[d]= Rand::gengam(LambdaArgs.h*Loci->GetDistance(locus), LambdaArgs.beta);
      if(Comms::isMaster()){
        SumLogLambda.push_back(0.0);
      }
      ++d;
      ++locus;
    }//end locus loop
  }//end chromosome loop
  
}

PopHapMix::~PopHapMix()
{
  delete[] MixtureProps;
  //delete[] MixturePropsproposal;
  delete[] HapMixLambdaSampler;
  delete[] MixturePropsPrior;
  delete[] dirparams;
}


/**
   samples locus-specific number of arrivals in a hapmixmodel, using Hamiltonian Monte Carlo.
   StateArrivalCounts is a (NumberOfCompositeLoci) * (NumBlockStates +1) array of counts of gametes with ancestry states that are
   accumulatelogs indicates whether to accumulate sums of log lambda.
*/
void PopHapMix::SampleHapMixLambda(const int* ConcordanceCounts, bool accumulateLogs){
  //double sum = 0.0;
  int locus = 0;
  int interval = 0;
  if(Comms::isMaster()){
    EventLogger::LogEvent(9, 0, "sampleLambda");
    LambdaPriorArgs.sumlambda = 0.0;
    hargs.sum_dloglambda = 0.0;
    try{
      vector<double>::iterator lambda_iter = lambda.begin();
      vector<double>::iterator sumloglambda_iter = SumLogLambda.begin();
 
      for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
	++locus;//skip first locus on each chromosome
	for(unsigned i = 1; i < Loci->GetSizeOfChromosome(c); ++i){
	  double loglambda = log(*lambda_iter);//sampler is on log scale
          LambdaArgs.theta = MixtureProps + locus*K;
	  LambdaArgs.Distance = hargs.distances[interval];//distance between this locus and last
	  //LambdaArgs.NumConcordant = StateArrivalCounts[locus*2+1];
	  //LambdaArgs.NumDiscordant = StateArrivalCounts[locus*2];
	  LambdaArgs.ConcordanceCounts = ConcordanceCounts + locus*2*K;
	  //sample new value
	  HapMixLambdaSampler[interval].Sample(&loglambda, &LambdaArgs);
	  *lambda_iter = exp(loglambda);

	  //accumulate sums of log of lambda
	  if(accumulateLogs)
	    *(sumloglambda_iter++) += loglambda;
	  //accumulate sums, used to sample beta
	  LambdaPriorArgs.sumlambda += *lambda_iter;
	  hargs.sum_dloglambda += hargs.distances[interval]*loglambda;

	  ++interval;
	  ++locus;
	  ++lambda_iter;
	}
      }
    }
    catch(string s){
      stringstream err;err << "Error encountered while sampling lambda " << interval << ":\n" + s;
      throw(err.str());
    }
    EventLogger::LogEvent(10, 0, "sampledLambda");
  }

//broadcast lambda to all processes, in Comm (excludes freqsampler)
//no need to broadcast MixtureProps if it is kept fixed
#ifdef PARALLEL
  EventLogger::LogEvent(17, 0, "Bcastlambda");
  Comms::BroadcastVector(lambda);
  EventLogger::LogEvent(18, 0, "Bcasted");
#endif    

  if(Comms::isMaster()){
      //SamplehRandomWalk();
      Sampleh_ARS();
      SampleRateParameter();
  }
}

void PopHapMix::Sampleh_RandomWalk(){
   try{
//sample h component of lambda shape parameter, using random walk
     const double h = LambdaArgs.h;
     const double logh = log(h);
     const double logproposal = Rand::gennor(logh, hTuner.getStepSize());
     const double proposal = exp(logproposal);
     
     double logLikelihoodRatio = (proposal - h)*(hargs.Dlogbeta + hargs.sum_dloglambda) + hargs.sum_lngamma_hd;
     double sum = 0.0;
     int locus = 0;
     for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
       ++locus;//skip first locus on each chromosome
       for(unsigned i = 1; i < Loci->GetSizeOfChromosome(c); ++i){
	       sum += lngamma(proposal*Loci->GetDistance(locus));
       }
       ++locus;
     }
     logLikelihoodRatio -= sum;
     double logPriorRatio = (hargs.shape)*(logproposal - logh) - hargs.rate*(proposal - h);
     //accept/reject step
     const double LogAccProbRatio = logLikelihoodRatio + logPriorRatio;
     bool accept = false;
     if( (LogAccProbRatio) < 0 ) {
       if( log(Rand::myrand()) < LogAccProbRatio ) accept = true;
     } else accept = true;  
     
     if(accept){
       LambdaArgs.h = proposal;
       hargs.sum_lngamma_hd = sum;//to save recalculating every time
     }
     //update sampler object every w updates
     //if( !( NumberOfUpdates % w ) ){
     hTuner.UpdateStepSize( exp(LogAccProbRatio) );  
     //}
   }
   catch(string s){
     string err = "Error encountered while sampling h in lambda prior:\n" + s;
     throw(err);
   }
}

void PopHapMix::Sampleh_ARS(){
  try{
    LambdaArgs.h = hARS.Sample(&hargs, hd2logf);
  }
   catch(string s){
     string err = "Error encountered while sampling h in lambda prior:\n" + s;
     throw(err);
   }
}

void PopHapMix::SampleRateParameter(){
  try{
//sample rate parameter with conjugate update
       const double pshape = LambdaArgs.beta_shape + LambdaArgs.h*Loci->GetLengthOfGenome();
       const double prate = LambdaArgs.beta_rate + LambdaPriorArgs.sumlambda;
       LambdaArgs.beta = Rand::gengam(pshape, prate);
       hargs.Dlogbeta = Loci->GetLengthOfGenome()*log(LambdaArgs.beta);
   }
   catch(string s){
     string err = "Error encountered while sampling lambda rate parameter:\n" + s;
     throw(err);
   }

}

///energy function for sampling locus-specific lambda, expected #arrivals
///conditional on locus ancestry and gamma prior params 
double PopHapMix::LambdaEnergy(const double* const x, const void* const vargs){
  //x is log lambda
  const LambdaArguments* args = (const LambdaArguments*)vargs;
  const double d = args->Distance;
  unsigned K = args->NumBlockStates;
  const double* theta = args->theta;
  const int* n = args->ConcordanceCounts;

  double E = 0.0;
  try {
    double lambda = myexp(*x);
    double f = myexp(-lambda);

    for(unsigned k = 0; k < K; ++k){
      E -= n[k]  * log(     theta[k]*(1.0 - f) )//discordant
	 + n[k+K]* log( f + theta[k]*(1.0 - f) );//concordant
    }

    // log unnormalized gamma prior on log lambda
    E -= args->h*d * (*x) - args->beta * lambda;
  } catch(string s){
    throw string("Error in LambdaEnergy: " + s);
  }
  return E; 
}

///gradient function for sampling locus-specific lambda, expected #arrivals
///conditional on locus ancestry and gamma prior params
void PopHapMix::LambdaGradient( const double* const x, const void* const vargs, double* g ){
  const LambdaArguments* args = (const LambdaArguments*)vargs;
  const double d = args->Distance;
  unsigned K = args->NumBlockStates;
  const double* theta = args->theta;
  const int* n = args->ConcordanceCounts;

  g[0] = 0.0;
  try {
    double lambda = myexp(*x);
    double f = myexp(-lambda);

    //first compute dE / df
    for(unsigned k = 0; k < K; ++k){
      g[0] += n[k]  / (1.0 - f);//discordant
      g[0] -= n[k+K]*(1.0 - theta[k]) / (f + theta[k]*(1.0 - f));//concordant
    }

    g[0] *= -lambda*f; // df / dlambda * dlambda / dx 
    
    //derivative of minus log prior wrt log lambda
    g[0] -= args->h*d  - args->beta * lambda;

  } catch(string s) {
    throw string("Error in LambdaGradient: " + s);
  }
}

//log posterior for h
double PopHapMix::hlogf(double h, const void* const vargs){
  const h_args* const args = (const h_args* const)vargs;
  double sum = 0.0;
  double loglikelihood = 0.0;
  double logprior = 0.0;
  try{
    for(unsigned j = 0; j < args->NumIntervals; ++j){
      sum += lngamma(h*args->distances[j]);
    }
    loglikelihood = h*(args->Dlogbeta + args->sum_dloglambda) - sum;
    logprior = (args->shape-1.0)*log(h) - args->rate*h;
  }
  catch(std::string s){
    throw(s);
  }
  return loglikelihood + logprior;
}

//derivative of log posterior
double PopHapMix::hdlogf(double h, const void* const vargs){
  const h_args* const args = (const h_args* const)vargs;
  double sum = 0.0;
  try{
    for(unsigned j = 0; j < args->NumIntervals; ++j){
      sum += args->distances[j]*digamma(h*args->distances[j]);
    }
  }
  catch(std::string s){
    throw(s);
  }
  return (args->Dlogbeta + args->sum_dloglambda) - sum + (args->shape-1.0)/h - args->rate;
}

//second derivative of log posterior
double PopHapMix::hd2logf(double h, const void* const vargs){
  const h_args* const args = (const h_args* const)vargs;
  double sum = 0.0;
  try{
    for(unsigned j = 0; j < args->NumIntervals; ++j){
      sum += args->distances[j]*args->distances[j]*trigamma(h*args->distances[j]);
    }
  }
  catch(std::string s){
    throw(s);
  }
  return  -sum - (args->shape-1.0)/(h*h);
}

void PopHapMix::InitializeOutputFile(const vector<string>& BlockStateLabels, const string& distanceUnit ) {
  if(Comms::isMaster()){
    // Header line of paramfile
    for(unsigned k = 0; k < K; ++k)
      outputstream << BlockStateLabels[k] << "\t";

    outputstream << "Dispersion\t"
		 << "ExpectedArrivals.Mean\tExpectedArrivals.Variance\t"
		 << "shapeParam\trateParam\t"
		 << "Exp.Arrivals.per"<< distanceUnit << endl;
  }
}

//output sample mean of ergodic average of lambda
void PopHapMix::OutputErgodicAvg( int samples, std::ofstream& avgstream)
{
  double sum = 0.0;
  for(vector<double>::const_iterator i = SumLogLambda.begin(); i < SumLogLambda.end(); ++i)sum += exp(*i / samples);
  avgstream.width(9);
  avgstream << setprecision(6) << sum / (double)lambda.size() << "\t";

}

///output to given output stream
void PopHapMix::OutputParams(ostream& out){
  out.width(9);

  //mixture props dispersion
  double sum = 0.0;
  for(unsigned k = 0; k < K; ++k)
    sum += MixturePropsPrior[k];

  out << sum << "\t";

  //sample mean and variance of lambdas
  sum = 0.0;
  double var = 0.0;
  
  //unsigned j = 1;
  for(vector<double>::const_iterator lambda_iter = lambda.begin(); lambda_iter != lambda.end(); ++lambda_iter){
    double r = *lambda_iter;// / Loci->GetDistance(j++);
    sum += r;
    var += r * r;
  }
  double size = (double)lambda.size();
  var = var - (sum*sum) / size;

  out << setiosflags(ios::fixed) << setprecision(6) << sum / size << "\t" << var /size << "\t"
    //lambda prior parameters
	 << LambdaArgs.h << "\t" << LambdaArgs.beta << "\t"
    //expected number of arrivals per unit distance
	 << LambdaArgs.h / LambdaArgs.beta << "\t" ;

}

void PopHapMix::OutputMixtureProps(ostream& out)const{
  vector<double> sum(K, 0.0);
  const unsigned L = Loci->GetNumberOfCompositeLoci();
  for(unsigned locus = 0; locus < L; ++locus){
    for(unsigned k = 0; k < K; ++k)
      sum[k] += MixtureProps[locus*K +k];
  }
  for(unsigned k = 0; k < K; ++k)
    out << sum[k] / (double)L << "\t";
}

void PopHapMix::OutputParams(int iteration, LogWriter &){
  //output to screen
  if( options->getDisplayLevel() > 2 )
    {
      OutputParams(cout);
    }
  //Output to paramfile after BurnIn
  if( iteration > options->getBurnIn() ){
    OutputMixtureProps(outputstream);
    OutputParams(outputstream);
    outputstream << endl;
  }
}

const double *PopHapMix::getGlobalTheta()const{
  return MixtureProps;
}

void PopHapMix::printAcceptanceRates(LogWriter &Log) {
  double av = 0;//average acceptance rate
  for(unsigned j = 0; j < lambda.size(); ++j){
    av += HapMixLambdaSampler[j].getAcceptanceRate();
  }
  Log << "Average Expected acceptance rate in samplers for numbers of arrivals: "
      << av / (double)(lambda.size());
  
  av = 0;//average stepsize
  for(unsigned j = 0; j < lambda.size(); ++j){
    av += HapMixLambdaSampler[j].getStepsize();
    //if(isnan(HapMixLambdaSampler[j].getStepsize()))cout << j << " " << HapMixLambdaSampler[j].getStepsize() << " " << HapMixLambdaSampler[j].getAcceptanceRate() << endl;  
  }
  Log << "\nwith average final step size of "
      << av / (double)(lambda.size())
      << "\n";
//   Log << "Expected acceptance rate in h sampler: " << hTuner.getExpectedAcceptanceRate()
//       << "\nwith final step size of " << hTuner.getStepSize() << "\n";

  Log << "Acceptance rate in sampler for mixture proportion dispersion: "
      << MixturePropsDispersionSampler.getExpectedAcceptanceRate()
      << "\nwith final step size of " << MixturePropsDispersionSampler.getStepSize() << "\n";
  
}
///Outputs h, beta and current values of lambda to file
void PopHapMix::OutputLambda(const char* filename)const{
  if(strlen(filename)){
    ofstream outfile(filename);
    //     int locus = 0;
    //     vector<double>::const_iterator r = rho.begin();
    //     for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
    // 	++locus;//skip first locus on each chromosome
    // 	for(unsigned i = 1; i < Loci->GetSizeOfChromosome(c); ++i){
    // 	    outfile << *r << " " << Loci->GetDistance(locus)<< endl;
    // 	    ++r;
    // 	    ++locus;
    // 	}
    //     }
    outfile << LambdaArgs.h << " " << LambdaArgs.beta << " ";
    copy(lambda.begin(), lambda.end(), ostream_iterator<double>(outfile, " "));

    
    outfile.close();
  }
}
void PopHapMix::OutputLambdaPosteriorMeans(const char* filename, int samples)const{
  ofstream outfile(filename);
  for(vector<double>::const_iterator i = SumLogLambda.begin(); i < SumLogLambda.end(); ++i)
      outfile << exp(*i / samples) << endl;
    outfile.close();
}

///sample mixture proportions with conjugate update
void PopHapMix::SampleMixtureProportions(const int* SumArrivalCounts){
  if(Comms::isMaster()){
    const unsigned L = Loci->GetNumberOfCompositeLoci();
    //TODO: ?? make SumLogTheta a class variable
    double SumLogTheta[K];

    for(size_t k = 0; k < K; ++k)
      SumLogTheta[k] = 0.0;
    
    for(unsigned j  = 0; j < L; ++j){
      for(size_t k = 0; k < K; ++k){
	dirparams[k] = MixturePropsPrior[k] + SumArrivalCounts[j*K + k];
      }
      Rand::gendirichlet(K, dirparams, MixtureProps+j*K );
      for(size_t k = 0; k < K; ++k){
	SumLogTheta[k] += log(MixtureProps[j*K +k]);
      }
    }

    //sample prior dispersion
    MixturePropsDispersionSampler.SampleEta(SumLogTheta, MixturePropsPrior);
  }

#ifdef PARALLEL
  //TODO: broadcast mixture props
#endif
}

///set global state arrival probs in hapmixmodel
//TODO: can skip this if xonly analysis with no females
//NB: assumes always diploid in hapmixmodel
void PopHapMix::SetHMMStateArrivalProbs(){
  //set locus correlation (workers only)
  if(Comms::isWorker()){
    Loci->SetLocusCorrelation(lambda);
    for( unsigned int j = 0; j < Loci->GetNumberOfChromosomes(); j++ ){
      //TODO: can probably replace isRandomMating with false
      Chromosome* C = Loci->getChromosome(j);
      MixturePropsWrapper MPW(MixtureProps + (C->GetLocus(0) * K), K); 
      C->SetHMMTheta(MPW, options->isRandomMatingModel(), true);
      C->SetStateArrivalProbs(options->isRandomMatingModel(), true);
    }
  }
}

