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
#include "bclib/misc.h"
#include "bclib/dist.h"//for log gamma density
#include <algorithm>
#include <numeric>
#include <sstream>
#include <iomanip>
#include "gsl/gsl_math.h"
#include "gsl/gsl_specfunc.h"
#include "bclib/GSLExceptions.h" //for catching gsl exceptions
#include "bclib/Exceptions.h"//for throwing InfiniteGradient to Hamiltonian
#include "bclib/RObjectWriter.h"

using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;
#define SQ( x ) x*x
#define REC(x) 1.0 / x

PopHapMix::PopHapMix( const HapMixOptions& op, HapMixGenome& loci):
  K(op.getNumberOfBlockStates()), options(op), Loci(loci){

  MixtureProps = 0;
  ArrivalRateSampler = 0;
  MixturePropsPrior = 0;
  dirparams = 0;
  SumTheta = 0;
  SumThetaSq = 0;
  eta = 0.0;
  fixRateParameter = false;
}

void PopHapMix::Initialise(const string& distanceUnit, LogWriter& Log){
  InitialiseMixtureProportions(Log);
  
  InitialiseArrivalRates(Log);
  
  // ** Open paramfile **
  Log.setDisplayMode(Quiet);
  if( strlen( options.getParameterFilename() ) ){
    outputstream.open( options.getParameterFilename(), ios::out );
    if( !outputstream )
      {
	Log << On << "ERROR: Couldn't open paramfile\n";
	exit( 1 );
      }
    else{
      Log << Quiet << "Writing population-level parameters to " << options.getParameterFilename() << "\n";
      InitializeOutputFile(distanceUnit);
    }
  }
  else{
    Log << Quiet << "No paramfile given\n";
  }
  
}

///initialise global admixture proportions
void PopHapMix::InitialiseMixtureProportions(LogWriter& Log){
  const unsigned L = Loci.GetNumberOfCompositeLoci();
  const unsigned KL = K*L;
  MixtureProps = new double[KL];

  //set initial values
  const string initfilename = options.getInitialMixturePropsFilename();
  //get initial values from file
  if(initfilename.size()){
    ReadInitialMixturePropsFromFile(initfilename.c_str(), Log);
  }
  //set initial values of Mixture Props as uniform  
  else{
    //TODO: initialise to prior means
    fill(MixtureProps, MixtureProps + K*L, 1.0/(double)K);
  }

  if(!options.getFixedMixtureProps()){//we are sampling mixture props
    /*
    Mixture proportions have a Dirichlet prior with parameters MixturePropsPrecision / K 
    ( so that the proportions are equal).
    */

    //allocate arrays required for sampling of mixture props
    MixturePropsPrior = new double[K];
    dirparams = new double[K];
    SumTheta = new double[K];
    SumThetaSq = new double[K];

    //set initial value of precision
    double InitialPrecision = options.getMixturePropsPrecision(); //user-specified
    if(InitialPrecision <= 0.0)InitialPrecision = (double)K;// default

    double MPDShape = (double)K;// precision prior shape
    double MPDRate = 1.0;//precision prior rate
    if(!options.getFixedMixturePropsPrecision()){//we are sampling precision
      /*
	The prior precision (MixturePropsPriorPrecision) has a 
	Gamma prior, which defaults to Gamma(1,1).
      */

      //check for user-specified prior
      const vector<double>& userprior = options.getMixturePropsPrecisionPrior();
      if(userprior.size()){//user has specified a prior
	MPDShape = userprior[0];
	MPDRate = userprior[1];
	//InitialPrecision = MPDShape / MPDRate;
      }

      //setup sampler for mixure props precision
      MixturePropsPrecisionSampler.SetSize(L, K, 0.01/*<-initial stepsize*/, 20/*<-num leapfrog steps*/);
      MixturePropsPrecisionSampler.SetPriorEta(MPDShape, MPDRate);
    }

    //set initial values of Dirichlet prior params
    double DirParamInit = InitialPrecision / (double)K;
    fill(MixturePropsPrior, MixturePropsPrior + K, DirParamInit);

    //write prior parameters to screen and log
    Log << Quiet << "Dirichlet prior on mixture proportions" ;
    if(!options.getFixedMixturePropsPrecision()){
      Log << "\nGamma(" << MPDShape << ", " << MPDRate << ") prior on mixture proportion prior precision\n" ;
    }
    else
      Log << " with precision: " << InitialPrecision << "\n";
    

  }//end if not fixed proportions
  else
    Log << Quiet << "Model with fixed mixture proportions\n";
  
}

void PopHapMix::ReadInitialMixturePropsFromFile(const char* initfilename, LogWriter& Log){
  ifstream initfile(initfilename);
  
  if(initfile.is_open()){
    Log << Quiet << "Reading initial values of mixture proportions from " << initfilename << "\n";
    const unsigned KL = K*Loci.GetNumberOfCompositeLoci();
   //read until vector is filled or end-of-file
    for(unsigned index = 0; index < KL; ++index){
      if(!(initfile >>MixtureProps[index])){
      //reached end of file before vector full -> file too small
      Log << On << "\nERROR: Too few entries in initialmixturepropsfile. Expected " << KL
	  << ", found " << index << "\n";
	exit(1);
      }
      if(MixtureProps[index] < 0.0 || MixtureProps[index] > 1.0)
	throw DataOutOfRangeException("mixture proportions", "between 0 and 1", "initialmixturepropsfile");
    }

    //see if anything more than whitespace left in file
    string test;
    initfile >> test;
    if(test.find_first_not_of(" \t\n\r") != string::npos){
      Log << On << "WERROR: Too many entries in initialmixturepropsfile. Expected " << KL << "\n";
      exit(1);
    }
    
    initfile.close();
  }
  else{
    string err("ERROR: cannot open initialmixturepropsfile: ");
    err.append(initfilename);
    throw err;
  }
}

void PopHapMix::InitialiseArrivalRates(LogWriter& Log){
 const unsigned L = Loci.GetNumberOfCompositeLoci();

  //set priors
  const vector<double>& priorparams = options.getLambdaPrior();
  if(priorparams.size()==3){
    fixRateParameter = true;
    LambdaArgs.beta = priorparams[2];
  }
  else{
    fixRateParameter = false;
    LambdaArgs.beta_shape = priorparams[2];
    LambdaArgs.beta_rate = priorparams[3];
    LambdaArgs.beta = priorparams[2]/priorparams[3];//initialise beta at prior mean
  }
  unsigned numIntervals = L-Loci.GetNumberOfChromosomes(); 
  lambda.assign(numIntervals, 0);
  //TODO: if distance is too large, leave lambda at zero and do not sample (as if a new chromosome)
   
  LambdaArgs.NumBlockStates = K;
  LambdaArgs.NumIntervals = numIntervals;
  
  //set up Hamiltonian sampler for lambda
  ArrivalRateSampler = new HamiltonianMonteCarlo[numIntervals];
  const vector<float>& lambdasamplerparams = options.getLambdaSamplerParams();
  size_t size = lambdasamplerparams.size();
  float initial_stepsize = size? lambdasamplerparams[0] : 0.06;
  float min_stepsize = size? lambdasamplerparams[1] : 0.0001;
  float max_stepsize = size? lambdasamplerparams[2] : 1.0;
  float target_acceptrate = size? lambdasamplerparams[3] : 0.9;
  int num_leapfrog_steps = size? (int)lambdasamplerparams[4] : 20;
  
  //set constant args for sampler for h
  hargs.shape = priorparams[0];
  hargs.rate = priorparams[1];
  hargs.Dlogbeta = Loci.GetLengthOfGenome()*log(LambdaArgs.beta);
  hargs.sum_lngamma_hd = 0.0;
  hargs.NumIntervals = numIntervals;
  hargs.distances = new double[numIntervals];
  LambdaArgs.h = hargs.shape / hargs.rate; // initialise h at prior mean 
  
  //write prior to screen and log
  Log << Quiet << "Gamma(h, ";
  if(fixRateParameter)
    Log << LambdaArgs.beta;
  else
    Log << "beta";
  Log << ") distribution on number of arrivals per Mb. \nGamma( " << hargs.shape << ", " << hargs.rate << " ) prior on h";
  if(!fixRateParameter)
    Log << " and Gamma( " << LambdaArgs.beta_shape << ", " << LambdaArgs.beta_rate << " ) prior on beta";
  Log << "\n";
  
  //Hamiltonian sampler
  for(unsigned j = 0; j < numIntervals; ++j){
    ArrivalRateSampler[j].SetDimensions(1, initial_stepsize, min_stepsize, max_stepsize, num_leapfrog_steps, 
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
  const string initfilename = options.getInitialArrivalRateFilename();
  const bool useinitfile = (initfilename.size() > 0);

  if(useinitfile){
    ReadInitialArrivalRatesFromFile(initfilename.c_str(), Log);
  }
  
  for(unsigned c = 0; c < Loci.GetNumberOfChromosomes(); ++c){
    ++locus;//skip first locus on each chromosome
    for(unsigned i = 1; i < Loci.GetSizeOfChromosome(c); ++i){
      hargs.sum_lngamma_hd += lngamma(LambdaArgs.h*Loci.GetDistance(locus));
      hargs.distances[d] = Loci.GetDistance(locus);
      if(!useinitfile)//initialise at prior mean
        lambda[d]= Rand::gengam(LambdaArgs.h*Loci.GetDistance(locus), LambdaArgs.beta);

      SumLogLambda.push_back(0.0);

      ++d;
      ++locus;
    }//end locus loop
  }//end chromosome loop
  
}

void PopHapMix::ReadInitialArrivalRatesFromFile(const char* initfilename, LogWriter& Log){
  ifstream initfile(initfilename);
  
  if(initfile.is_open()){
    Log << Quiet << "Reading initial values of arrival rates from " << initfilename << "\n";
    
    //read initial values of h, beta
    initfile >> LambdaArgs.h >> LambdaArgs.beta;
    if(LambdaArgs.h <= 0.0)
      throw DataOutOfRangeException("Arrival Rate mean", ">0", "initialarrivalratefile");
    if(LambdaArgs.beta <= 0.0)
      throw DataOutOfRangeException("Arrival Rate rate parameter", ">0", "initialarrivalratefile");
    
   //read until vector is filled or end-of-file
    for(vector<double>::iterator ar = lambda.begin(); ar != lambda.end(); ++ar){
      if(!(initfile >>*ar)){
      //reached end of file before vector full -> file too small
	Log << On << "\nERROR: Too few entries in initialarrivalratefile. Expected " 
	    << (unsigned)(lambda.size() + 2) << ", found " 
	    << (unsigned)(lambda.size() - (lambda.end() - ar)+2) << "\n";
	exit(1);
      }
      if(*ar <= 0.0)
	throw DataOutOfRangeException("Arrival Rate", ">0", "initialarrivalratefile");
    }

    //see if anything more than whitespace left in file
    string test;
    initfile >> test;
    if(test.find_first_not_of(" \t\n\r") != string::npos){
      Log << On << "ERROR: Too many entries in initialarrivalratefile\n";
      exit(1);
    }
    initfile.close(); 
  }
  else{
    string err("ERROR: cannot open initialarrivalratefile: ");
    err.append(initfilename);
    throw err;
  }
}

PopHapMix::~PopHapMix(){
  delete[] MixtureProps;
  delete[] ArrivalRateSampler;
  delete[] MixturePropsPrior;
  delete[] dirparams;
  delete[] SumTheta;
  delete[] SumThetaSq;
}


/**
   samples locus-specific number of arrivals in a hapmixmodel, using Hamiltonian Monte Carlo.
   StateArrivalCounts is a (NumberOfCompositeLoci) * (NumBlockStates +1) array of counts of gametes with ancestry states that are
   accumulatelogs indicates whether to accumulate sums of log lambda.
*/
void PopHapMix::SampleArrivalRate(const int* ConcordanceCounts, bool accumulateLogs){
  //double sum = 0.0;
  int locus = 0;
  int interval = 0;
  LambdaPriorArgs.sumlambda = 0.0;
  hargs.sum_dloglambda = 0.0;
  try{
    vector<double>::iterator lambda_iter = lambda.begin();
    vector<double>::iterator sumloglambda_iter = SumLogLambda.begin();
    
    for(unsigned c = 0; c < Loci.GetNumberOfChromosomes(); ++c){
      ++locus;//skip first locus on each chromosome
      for(unsigned i = 1; i < Loci.GetSizeOfChromosome(c); ++i){
	double loglambda = log(*lambda_iter);//sampler is on log scale
	LambdaArgs.theta = MixtureProps + locus*K;
	LambdaArgs.Distance = hargs.distances[interval];//distance between this locus and last
	//LambdaArgs.NumConcordant = StateArrivalCounts[locus*2+1];
	//LambdaArgs.NumDiscordant = StateArrivalCounts[locus*2];
	LambdaArgs.ConcordanceCounts = ConcordanceCounts + locus*2*K;
	//sample new value
	ArrivalRateSampler[interval].Sample(&loglambda, &LambdaArgs);
	*lambda_iter = exp(loglambda);
	
	// 	  if(*lambda_iter > 50.0)
	// 	    {cout << interval << " " << *lambda_iter << endl;
	// 	      system("pause");
	// 	    }
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
  
  //SamplehRandomWalk();
  Sampleh_ARS();
  if(!fixRateParameter)
    SampleRateParameter();
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
     for(unsigned c = 0; c < Loci.GetNumberOfChromosomes(); ++c){
       ++locus;//skip first locus on each chromosome
       for(unsigned i = 1; i < Loci.GetSizeOfChromosome(c); ++i){
	       sum += lngamma(proposal*Loci.GetDistance(locus));
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
    const double pshape = LambdaArgs.beta_shape + LambdaArgs.h*Loci.GetLengthOfGenome();
    const double prate = LambdaArgs.beta_rate + LambdaPriorArgs.sumlambda;
    LambdaArgs.beta = Rand::gengam(pshape, prate);
    hargs.Dlogbeta = Loci.GetLengthOfGenome()*log(LambdaArgs.beta);
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
      E -= n[k]  * log(              theta[k]*(1.0 - f) )//discordant, omitting constant term
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
  double f = 0.0;
  double lambda = 0.0;
  try {
    lambda = myexp(*x);
    f = myexp(-lambda);
  } catch(string s) {
    throw string("Error in LambdaGradient: " + s);
  }
  catch(overflow e){
    throw InfiniteGradient("LambdaGradient", "PopHapMix.cc");
  }

  //first compute dE / df
  for(unsigned k = 0; k < K; ++k){
    g[0] += n[k]  / (1.0 - f);//discordant
    g[0] -= n[k+K]*(1.0 - theta[k]) / (f + theta[k]*(1.0 - f));//concordant
  }
  
  g[0] *= -lambda*f; // df / dlambda * dlambda / dx 
  
  //derivative of minus log prior wrt log lambda
  g[0] -= args->h*d  - args->beta * lambda;
  
  //     if(!gsl_finite(g[0]))
  //       throw InfiniteGradient("LambdaGradient", "PopHapMix.cc");
}

//log posterior for h
double PopHapMix::hlogf(double h, const void* const vargs){
  const h_args* const args = (h_args*)vargs;
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
  const h_args* const args = (h_args*)vargs;
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
  const h_args* const args = (h_args*)vargs;
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

void PopHapMix::InitializeOutputFile(const string& distanceUnit ) {
  // Header line of paramfile
  if(!options.getFixedMixturePropsPrecision())
    outputstream << "MixtureProps.Precision\t";

  if(!options.getFixedMixtureProps())
    outputstream << "MixtureProps.Sample.Precision\t";

  outputstream << "Arrivals.per"<< distanceUnit << ".shapeParam\t";

  if(!fixRateParameter)
    outputstream << "Arrivals.per"<< distanceUnit << ".rateParam\t";

  outputstream << "Arrivals.per"<< distanceUnit << ".Mean"
	       << "\tMean.Block.Length.kb"
	       << endl;
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
  out << setiosflags(ios::fixed) << setprecision(6);

  //mixture props precision
  if(!options.getFixedMixturePropsPrecision()){
    double sum = 0.0;
    for(unsigned k = 0; k < K; ++k)
      sum += MixturePropsPrior[k];
    out << sum << "\t";
  }
  
  if(!options.getFixedMixtureProps())
    //observed precision
    out << eta  << "\t";
 
  out
    //lambda prior parameters
    << LambdaArgs.h << "\t";
  if(!fixRateParameter) 
    out<< LambdaArgs.beta << "\t";
  //expected number of arrivals per unit distance
  out << LambdaArgs.h / LambdaArgs.beta << "\t" ;
}

//output mixture proportions averaged over loci
void PopHapMix::OutputAverageMixtureProps(ostream& out)const{
  vector<double> sum(K, 0.0);
  const unsigned L = Loci.GetNumberOfCompositeLoci();
  for(unsigned locus = 0; locus < L; ++locus){
    for(unsigned k = 0; k < K; ++k)
      sum[k] += MixtureProps[locus*K +k];
  }
  for(unsigned k = 0; k < K; ++k)
    out << sum[k] / (double)L << "\t";
}

void PopHapMix::OutputParams(int iteration, LogWriter &){
  //output to screen
  if( options.getDisplayLevel() > 2 )
    {
      OutputParams(cout);
    }
  //Output to paramfile after BurnIn
  if( iteration > options.getBurnIn() ){
    // OutputAverageMixtureProps(outputstream);
    OutputParams(outputstream);
    //write Haplotype block length to file
    outputstream << (8000 * LambdaArgs.beta) / (LambdaArgs.h * 7)
		 << endl;
  }
}

const double *PopHapMix::getGlobalMixtureProps()const{
  return MixtureProps;
}

void PopHapMix::printAcceptanceRates(LogWriter &Log) {
  double av = 0;//average acceptance rate
  for(unsigned j = 0; j < lambda.size(); ++j){
    av += ArrivalRateSampler[j].getAcceptanceRate();
  }
  Log << "Average Expected acceptance rate in samplers for numbers of arrivals: "
      << av / (double)(lambda.size());
  
  av = 0;//average stepsize
  for(unsigned j = 0; j < lambda.size(); ++j){
    av += ArrivalRateSampler[j].getStepsize();
    //if(isnan(ArrivalRateSampler[j].getStepsize()))cout << j << " " << ArrivalRateSampler[j].getStepsize() << " " << ArrivalRateSampler[j].getAcceptanceRate() << endl;  
  }
  Log << "\nwith average final step size of "
      << av / (double)(lambda.size())
      << "\n";
//   Log << "Expected acceptance rate in h sampler: " << hTuner.getExpectedAcceptanceRate()
//       << "\nwith final step size of " << hTuner.getStepSize() << "\n";

  if(!options.getFixedMixturePropsPrecision()){
    Log << "Acceptance rate in sampler for mixture proportion precision: "
	<< MixturePropsPrecisionSampler.getExpectedAcceptanceRate()
	<< "\nwith final step size of " << MixturePropsPrecisionSampler.getStepSize() << "\n";
  }
  
}
///Outputs h, beta and current values of lambda to file
void PopHapMix::OutputArrivalRates(const char* filename)const{
  if(strlen(filename)){
    ofstream outfile(filename);
    if(outfile.is_open()){
      //write prior params first
      outfile << LambdaArgs.h << " " << LambdaArgs.beta << " ";
      //then copy all lambdas, space-separated
      copy(lambda.begin(), lambda.end(), ostream_iterator<double>(outfile, " "));
      
      outfile.close();
    }
    else
      cerr << "Error: cannot open finalarrivalratefile: " << filename << endl;
  }
}

///outputs mixture props to file
void PopHapMix::OutputMixtureProps(const char* filename)const{
  if(strlen(filename)){
    ofstream outfile(filename);
    if(outfile.is_open()){
      //copy to output stream , space-separated
      copy(MixtureProps, MixtureProps + K*Loci.GetNumberOfCompositeLoci(), ostream_iterator<double>(outfile, " "));
      outfile.close();
    }
    else
      cerr << "Error: cannot open finalmixturepropsfile: " << filename << endl;
  }
}

///output posterior means of arrival rates per unit distance, as R object
void PopHapMix::OutputArrivalRatePosteriorMeans(const char* filename, int samples, const string& distanceUnit)const{
  if(filename == 0 || strlen(filename)==0)
    return;
  RObjectWriter outfile(filename);

  unsigned d = 0;//to index distances
  vector<double>::const_iterator i = SumLogLambda.begin();
  for(unsigned c = 0; c < Loci.GetNumberOfChromosomes(); ++c){
    //write NA for first position
    outfile << "NA";
    for(unsigned j  = 1; j < Loci.GetSizeOfChromosome(c); ++j){//step through rest of loci on chromosome
      //distances are stored in hargs
      outfile << exp(*i / samples) /hargs.distances[d];
      ++i;
      ++d;
      }
  }
  vector<vector< string> > dimnames(2);
  dimnames[1].push_back("ArrivalRatePer" + distanceUnit);
  
  vector<int> dims(2);
  dims[0] = Loci.GetNumberOfCompositeLoci();
  dims[1] = 1;
  
  outfile.close(dims, dimnames);

}

///sample mixture proportions with conjugate update
void PopHapMix::SampleMixtureProportions(const int* SumArrivalCounts){
  if(!options.getFixedMixtureProps()){
    const unsigned L = Loci.GetNumberOfCompositeLoci();
    //TODO: ?? make SumLogTheta a class variable
    double SumLogTheta[K];

    for(size_t k = 0; k < K; ++k){
      SumLogTheta[k] = 0.0;
      SumTheta[k] = 0.0;
      SumThetaSq[k] = 0.0;
    }
    
    for(unsigned j  = 0; j < L; ++j){
      //set Dirichlet paramaters
      for(size_t k = 0; k < K; ++k){
	dirparams[k] = MixturePropsPrior[k] + SumArrivalCounts[j*K + k];
      }
      //sample from Dirichlet with those parameters
      Rand::gendirichlet(K, dirparams, MixtureProps+j*K );

      //accumulate sums, sums-of-squares and sums-of-logs
      for(size_t k = 0; k < K; ++k){
	SumTheta[k] += MixtureProps[j*K +k];
	SumThetaSq[k] += SQ( MixtureProps[j*K +k] );
	SumLogTheta[k] += log(MixtureProps[j*K +k]);
      }
    }

    if(!options.getFixedMixturePropsPrecision()){
      //sample prior precision
      MixturePropsPrecisionSampler.SampleEta(SumLogTheta, MixturePropsPrior);
    }

    //calculate observed variance of mixture props
    double SumThetaVar = 0.0;
    for(size_t k = 0; k < K; ++k){
      SumThetaVar += SumThetaSq[k] / (double) L  - SQ( (SumTheta[k] / (double) L) );
    }

    eta = ( (K - 1.0) / ((double)K * SumThetaVar )) - 1.0 ;

  }
}

/**
   set global state arrival probs in hapmixmodel.
    NB: assumes always diploid
*/
//TODO: can skip this if xonly analysis with no females
void PopHapMix::SetHMMStateArrivalProbs(bool diploid){
  //set locus correlation (workers only)
  Loci.SetLocusCorrelation(lambda);
  for( unsigned int j = 0; j < Loci.GetNumberOfChromosomes(); j++ ){
    Loci.getChromosome(j)->HMM->SetStateArrivalProbs(MixtureProps, 0, diploid);
  }
}

