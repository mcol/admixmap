/** 
 *   ADMIXMAP
 *   PopHapMix.cc 
 *   Class to hold and update population parameters in a hapmixmodel
 *   Copyright (c) 2006 David O'Donnell and Paul McKeigue
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
#ifdef PARALLEL
#include <mpe.h>//for MPI event logging
#endif
using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

PopHapMix::PopHapMix( AdmixOptions* op, Genome* loci)
{
  options = op;
  Loci = loci;
  globaltheta = 0;
  //globalthetaproposal = 0;
  HapMixLambdaSampler = 0;

}

void PopHapMix::Initialise(int , const Vector_s& PopulationLabels, LogWriter &Log){
  Log.setDisplayMode(On);
  K = options->getPopulations();

  if(K > 1){
    //initialise global admixture proportions
    globaltheta = new double[K];
    //globalthetaproposal = new double[K];
    fill(globaltheta, globaltheta+K, 1.0/(double)K);
    //ThetaTuner.SetParameters(1.0 /*<-initial stepsize on softmax scale*/, 0.00, 10.0, 0.44);
    
    //if(Comms::isMaster())SumLogLambda.push_back(0.0);
    //set priors
    const vector<double>& priorparams = options->getHapMixLambdaPrior();
    
    LambdaArgs.beta_shape = priorparams[2];
    LambdaArgs.beta_rate = priorparams[3];
    LambdaArgs.beta = priorparams[2]/priorparams[3];//initialise beta at prior mean
    
    unsigned numIntervals = Loci->GetNumberOfCompositeLoci()-Loci->GetNumberOfChromosomes();
    if(Comms::isMaster()){
      LambdaArgs.NumPops = K;
      LambdaArgs.NumIntervals = numIntervals;
      
      //set up Hamiltonian sampler for lambda
      HapMixLambdaSampler = new HamiltonianMonteCarlo[numIntervals];
      const vector<float>& lambdasamplerparams = options->getrhoSamplerParams();
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
      
    }//end sampler initialisation
    //initialise lambda vector
    int locus = 0;
    int d = 0;
    const char* initfilename = options->getInitialHapMixLambdaFilename();
    const bool useinitfile = (strlen(initfilename) > 0);
    double initvalue;
    ifstream initfile;
    if(useinitfile)initfile.open(initfilename);

    for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
      ++locus;//skip first locus on each chromosome
      for(unsigned i = 1; i < Loci->GetSizeOfChromosome(c); ++i){
	hargs.sum_lngamma_hd += lngamma(LambdaArgs.h*Loci->GetDistance(locus));
	hargs.distances[d++] = Loci->GetDistance(locus);
	if(useinitfile){
	    initfile >> initvalue;
	    lambda.push_back(initvalue);
	}
	else
	    lambda.push_back(Rand::gengam(LambdaArgs.h*Loci->GetDistance(locus), LambdaArgs.beta));
	if(Comms::isMaster()){
	  SumLogLambda.push_back(0.0);
	}
	++locus;
      }
    }
    if(useinitfile)initfile.close();
    
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
  }//end if K > 1
}

PopHapMix::~PopHapMix()
{
  delete[] globaltheta;
  //delete[] globalthetaproposal;
  delete[] HapMixLambdaSampler;
}


///samples locus-specific number of arrivals in a hapmixmodel, using Hamiltonian Monte Carlo
void PopHapMix::SampleHapMixLambda(const int* SumAncestry, bool accumulateLogs){
  //SumAncestry is a (NumberOfCompositeLoci) * 2 array of counts of gametes with ancestry states that are
  // unequal (element 0) and equal (element 1)
  //accumulatelogs indicates whether to accumulate sums of log lambda

  //LambdaArgs.theta = globaltheta;
  //double sum = 0.0;
  int locus = 0;
  int interval = 0;
  if(Comms::isMaster()){
#ifdef PARALLEL
    MPE_Log_event(9, 0, "sampleLambda");
#endif
    LambdaPriorArgs.sumlambda = 0.0;
    hargs.sum_dloglambda = 0.0;
    try{
      vector<double>::iterator lambda_iter = lambda.begin();
      vector<double>::iterator sumloglambda_iter = SumLogLambda.begin();
      for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
	++locus;//skip first locus on each chromosome
	for(unsigned i = 1; i < Loci->GetSizeOfChromosome(c); ++i){
	  double loglambda = log(*lambda_iter);//sampler is on log scale
	  
	  LambdaArgs.Distance = hargs.distances[interval];//distance between this locus and last
	  LambdaArgs.NumConcordant = SumAncestry[locus*2+1];
	  LambdaArgs.NumDiscordant = SumAncestry[locus*2];

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
#ifdef PARALLEL
    MPE_Log_event(10, 0, "sampledLambda");
#endif
  }
//broadcast lambda to all processes, in Comm (excludes freqsampler)
//no need to broadcast globaltheta if it is kept fixed
#ifdef PARALLEL
  MPE_Log_event(17, 0, "Bcastlambda");
  Comms::BroadcastVector(lambda);
  MPE_Log_event(18, 0, "Bcasted");
#endif    
  //set locus correlation (workers only)
  if(Comms::isWorker())  
    {
      Loci->SetLocusCorrelation(lambda);
      for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
	//set global state arrival probs in hapmixmodel
	//TODO: can skip this if xonly analysis with no females
	Loci->getChromosome(c)->SetStateArrivalProbs(options->isRandomMatingModel(), true);
      }
    }
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
  unsigned K = args->NumPops;
  const double d = args->Distance;
  double theta = 1.0 / (double)K;
  double E = 0.0;
  try {
    double lambda = myexp(*x);
    double f = myexp(-lambda);
    double probequal = theta + f*(1.0 - theta);
    E -= args->NumDiscordant * log(1.0-probequal)//constant term in log(1-theta) omitted
	+ args->NumConcordant * log(probequal);

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
  unsigned K = args->NumPops;
  const double d = args->Distance;
  double theta = 1.0 / (double)K;
  g[0] = 0.0;
  try {
    double lambda = myexp(*x);
    double f = myexp(-lambda);
    //first compute dE / dprobequal
    double probequal = theta + f*(1.0 - theta);

    g[0] +=  args->NumDiscordant / (1.0 - probequal) - args->NumConcordant / probequal;//dE / dprobequal
    g[0] *= -lambda*f*(1.0-theta); // dprobequal / df * df / dlambda * dlambda / dx 
    
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

void PopHapMix::InitializeOutputFile(const Vector_s& ) {
  if(Comms::isMaster()){
    // Header line of paramfile
    outputstream << "Lambda.Mean\tLambda.Variance\th\tbeta";
    outputstream << endl;
  }
}

void PopHapMix::OutputErgodicAvg( int samples, std::ofstream *avgstream)
{
  double sum = 0.0;
  for(vector<double>::const_iterator i = SumLogLambda.begin(); i < SumLogLambda.end(); ++i)sum += exp(*i / samples);
  avgstream->width(9);
  *avgstream << setprecision(6) << sum / (double)lambda.size() << "\t";

}

//output to given output stream
void PopHapMix::OutputParams(ostream* out){
  //lambda
  out->width(9);
  double sum = 0.0;
  double var = 0.0;
  
  //unsigned j = 1;
  for(vector<double>::const_iterator lambda_iter = lambda.begin(); lambda_iter != lambda.end(); ++lambda_iter){
    double r = *lambda_iter;// / Loci->GetDistance(j++);
    sum += r;
    var += r * r;
  }
  double size = (double)lambda.size();
  var = var - (sum*sum) / size;
  (*out) << setiosflags(ios::fixed) << setprecision(6) << sum / size << "\t" << var /size << "\t"
	 << LambdaArgs.h << "\t" << LambdaArgs.beta << "\t" ;
}

void PopHapMix::OutputParams(int iteration, LogWriter &){
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

const double *PopHapMix::getGlobalTheta()const{
  return globaltheta;
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
  
}
void PopHapMix::OutputLambda(const char* filename)const{
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
     for(vector<double>::const_iterator i = lambda.begin(); i != lambda.end(); ++i){
 	outfile << *i << "\t";
     }

    outfile.close();
}
void PopHapMix::OutputLambdaPosteriorMeans(const char* filename, int samples)const{
  ofstream outfile(filename);
  for(vector<double>::const_iterator i = SumLogLambda.begin(); i < SumLogLambda.end(); ++i)
      outfile << exp(*i / samples) << endl;
    outfile.close();
}
