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
#include "gsl/gsl_math.h"
#include "gsl/gsl_specfunc.h"
#include "Comms.h"
#ifdef PARALLEL
#include <mpe.h>//for MPI event logging
#endif
using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

PopHapMix::PopHapMix( AdmixOptions* op, Genome* loci): Population(op, loci)
{
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
    
    if(Comms::isMaster())SumLogLambda.push_back(0.0);
    //set priors
    const vector<double>& priorparams = options->getHapMixLambdaPrior();
    
    LambdaArgs.h = priorparams[0]; // h
    LambdaArgs.beta_shape = priorparams[1];
    LambdaArgs.beta_rate = priorparams[2];
    LambdaArgs.beta = priorparams[1]/priorparams[2];//initialise beta at prior mean
    
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
      
      //set constant args for random walk sampler for h
      h_args.h_shape = 1000.0;
      h_args.h_rate = 1.0;
      h_args.Dlogbeta = Loci->GetLengthOfGenome()*log(LambdaArgs.beta);
      h_args.sum_lngamma_hd = 0.0;
      
      //Hamiltonian sampler
      for(unsigned j = 0; j < numIntervals; ++j){
	HapMixLambdaSampler[j].SetDimensions(1, initial_stepsize, min_stepsize, max_stepsize, num_leapfrog_steps, 
					     target_acceptrate, LambdaEnergy, LambdaGradient);
      }
      //Random Walk sampler for log h
      hTuner.SetParameters( 0.01, 0.0001, 1000.0, 0.44);
      
    }//end sampler initialisation
    //initialise lambda vector
    int locus = 0;
    const double initial_lambda = options->getInitialHapMixLambda();
    for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
      ++locus;//skip first locus on each chromosome
      for(unsigned i = 1; i < Loci->GetSizeOfChromosome(c); ++i){
	h_args.sum_lngamma_hd += lngamma(LambdaArgs.h*Loci->GetDistance(locus));
	if(initial_lambda >0.0)
	  lambda.push_back(initial_lambda);
	else
	  lambda.push_back(Rand::gengam(LambdaArgs.h*Loci->GetDistance(locus), LambdaArgs.beta));
	if(Comms::isMaster()){
	  SumLogLambda.push_back(0.0);
	}
	++locus;
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
  int index = 0;
  
  if(Comms::isMaster()){
#ifdef PARALLEL
    MPE_Log_event(9, 0, "sampleLambda");
#endif
    LambdaPriorArgs.sumlambda = 0.0;
    h_args.sum_dloglambda = 0.0;
    try{
      vector<double>::iterator lambda_iter = lambda.begin();
      vector<double>::iterator sumloglambda_iter = SumLogLambda.begin();
      for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
	++locus;//skip first locus on each chromosome
	for(unsigned i = 1; i < Loci->GetSizeOfChromosome(c); ++i){
	  double loglambda = log(*lambda_iter);//sampler is on log scale
	  
	  LambdaArgs.Distance = Loci->GetDistance(locus);//distance between this locus and last
	  LambdaArgs.SumAncestry = SumAncestry + locus*2;//sums of ancestry for this locus
	  
	  //sample new value
	  HapMixLambdaSampler[index++].Sample(&loglambda, &LambdaArgs);
	  *lambda_iter = exp(loglambda);

	  //accumulate sums of log of lambda
	  if(accumulateLogs)
	    *(sumloglambda_iter++) += loglambda;
	  //accumulate sums, used to sample beta
	  LambdaPriorArgs.sumlambda += *lambda_iter;
	  h_args.sum_dloglambda += Loci->GetDistance(locus)*loglambda;

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
	//NB: assumes always diploid in hapmixmodel
	Loci->getChromosome(c)->SetStateArrivalProbs(options->isRandomMatingModel());
      }
    }
  if(Comms::isMaster()){
      SampleHapmixLambdaPriorParameters();
  }
}

void PopHapMix::SampleHapmixLambdaPriorParameters(){
   try{
//sample h component of lambda shape parameter, using random walk
       const double h = LambdaArgs.h;
       const double logh = log(h);
       const double logproposal = Rand::gennor(logh, hTuner.getStepSize());
       const double proposal = exp(logproposal);

       double logLikelihoodRatio = (proposal - h)*(h_args.Dlogbeta + h_args.sum_dloglambda) + h_args.sum_lngamma_hd;
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
       double logPriorRatio = (h_args.h_shape)*(logproposal - logh) - h_args.h_rate*(proposal - h);
//accept/reject step
       const double LogAccProbRatio = logLikelihoodRatio + logPriorRatio;
       bool accept = false;
       if( (LogAccProbRatio) < 0 ) {
	   if( log(Rand::myrand()) < LogAccProbRatio ) accept = true;
       } else accept = true;  
       
       if(accept){
	   LambdaArgs.h = proposal;
	   h_args.sum_lngamma_hd = sum;//to save recalculating every time
       }
      //update sampler object every w updates
      //if( !( NumberOfUpdates % w ) ){
       hTuner.UpdateStepSize( exp(LogAccProbRatio) );  
	//}


//sample rate parameter with conjugate update
       const double pshape = LambdaArgs.beta_shape + LambdaArgs.h*Loci->GetLengthOfGenome();
       const double prate = LambdaArgs.beta_rate + LambdaPriorArgs.sumlambda;
       LambdaArgs.beta = Rand::gengam(pshape, prate);
       h_args.Dlogbeta = Loci->GetLengthOfGenome()*LambdaArgs.beta;
   }
   catch(string s){
     string err = "Error encountered while sampling lambda prior params:\n" + s;
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
    int sumequal = args->SumAncestry[1], sumnotequal = args->SumAncestry[0];
    double probequal = theta + f*(1.0 - theta);
    E -= sumnotequal * log(1.0-probequal);//constant term in log(1-theta) omitted
    E -= sumequal * log(probequal);

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
    int sumequal = args->SumAncestry[1], sumnotequal = args->SumAncestry[0];
    double probequal = theta + f*(1.0 - theta);

    g[0] +=  sumnotequal / (1.0 - probequal) - sumequal / probequal;//dE / dprobequal
    g[0] *= -lambda*f*(1.0-theta); // dprobequal / df * df / dlambda * dlambda / dx 
    
    //derivative of minus log prior wrt log lambda
    g[0] -= args->h*d  - args->beta * lambda;

  } catch(string s) {
    throw string("Error in LambdaGradient: " + s);
  }
}

void PopHapMix::InitializeOutputFile(const Vector_s& ) {
  if(Comms::isMaster()){
    // Header line of paramfile
    outputstream << "LambdaIntensities.Mean\tLambda.Variance\th\tbeta";
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
  //sumintensities
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
  Log << "Average Expected acceptance rate in sumintensities samplers: "
      << av / (double)(lambda.size());
  
  av = 0;//average stepsize
  for(unsigned j = 0; j < lambda.size(); ++j){
    av += HapMixLambdaSampler[j].getStepsize();
    //if(isnan(HapMixLambdaSampler[j].getStepsize()))cout << j << " " << HapMixLambdaSampler[j].getStepsize() << " " << HapMixLambdaSampler[j].getAcceptanceRate() << endl;  
  }
  Log << "\nwith average final step size of "
      << av / (double)(lambda.size())
      << "\n";
  Log << "Expected acceptance rate in h sampler: " << hTuner.getExpectedAcceptanceRate()
      << "\nwith final step size of " << hTuner.getStepSize() << "\n";
  
}
