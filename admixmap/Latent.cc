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
  RhoSampler = 0;

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
      //set prior of rho prior mean and variance
      RhoPriorArgs.priormeans = &( options->getHapMixRhoPriorMeans()[0]);
      //RhoPriorArgs.priorvars = &( options->getHapMixRhoPriorVars()[0]);

      rhoalpha = RhoPriorArgs.priormeans[0];
      rhobeta0 = RhoPriorArgs.priormeans[1];
      rhobeta1 = RhoPriorArgs.priormeans[2];
      rhobeta = rhobeta0 / rhobeta1;

      unsigned numIntervals = Loci->GetNumberOfCompositeLoci()-Loci->GetNumberOfChromosomes();
      if(Comms::isMaster()){
	RhoArgs.NumPops = K;
	RhoArgs.rhoalpha = rhoalpha;
	RhoArgs.rhobeta = rhobeta;
	RhoArgs.rhobeta0 = rhobeta0;
	RhoArgs.rhobeta1 = rhobeta1;

	RhoArgs.NumIntervals = numIntervals;
	//RhoArgs.sumrho  = numIntervals * rho[0];
	//RhoArgs.sumlogrho = numIntervals * log(rho[0]); // initial value
	
	//set up Hamiltonian sampler for rho
 	RhoSampler = new HamiltonianMonteCarlo[numIntervals];
	const vector<float>& rhosamplerparams = options->getrhoSamplerParams();
	size_t size = rhosamplerparams.size();
	float initial_stepsize = size? rhosamplerparams[0] : 0.05;
	float min_stepsize = size? rhosamplerparams[1] : 0.0001;
	float max_stepsize = size? rhosamplerparams[2] : 1.0;
	float target_acceptrate = size? rhosamplerparams[3] : 0.3;
	int num_leapfrog_steps = size? (int)rhosamplerparams[4] : 20;
	NumberOfUpdates = 0;
	w = 1;
	
//random-walk sampler
	TuneRhoSampler.SetParameters( initial_stepsize, min_stepsize, max_stepsize, target_acceptrate);
	step = initial_stepsize;

//Hamiltonian sampler
  	for(unsigned j = 0; j < numIntervals; ++j){
  	  RhoSampler[j].SetDimensions(1, initial_stepsize, min_stepsize, max_stepsize, num_leapfrog_steps, 
  				      0.9, RhoEnergy, RhoGradient);
 	}

 	RhoPriorArgs.NumIntervals = numIntervals;
 	RhoPriorArgs.rho = &rho; // pointer to vector<double>
	RhoAlphaSampler.SetParameters(0.1, 0.01, 100.0, 0.26, 1); 
	
      }//end sampler initialisation
      //initialise rho vector
      //double initial_rho = rhoalpha / rhobeta;
      double initial_rho = 10000.0;
      rho[0] = initial_rho;//Rand::gengam(rhoalpha, rhobeta);
      RhoPriorArgs.sumlogrho = log(rho[0]);
      for(unsigned j = 0; j < numIntervals-1; ++j){
	rho.push_back(initial_rho/*Rand::gengam(rhoalpha, rhobeta)*/);
	if(Comms::isMaster()){
	    SumLogRho.push_back(0.0);
	    RhoPriorArgs.sumlogrho += log(rho[j+1]);
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
    delete[] RhoSampler;
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
/**
   Updates locus-specific sumintensities in a hapmixmodel by moving all values by the same random step and using
   a Metropolis-Hastings accept/reject step.
*/
// void Latent::UpdateSumIntensitiesByRandomWalk(const IndividualCollection* const IC,bool sumlogrho){
//   NumberOfUpdates++;
//   double LogLikelihood = 0.0;
//   double LogLikelihoodAtProposal = 0.0;
//   double LogLikelihoodRatio = 0.0;
//   double LogPriorRatio = 0.0;
//   double LogAccProbRatio = 0.0;
//   bool accept = false;
//   double newsumlogrho = 0.0;
//   double newsumrho = 0.0;

//   double randomStep = Rand::gennor(0.0, step);
//   for(unsigned i = 0; i < rho.size(); ++i){
//       double logrho = log(rho[i]);
//       double logrhoproposal = logrho + randomStep;
//       rhoproposal[i] = exp ( logrhoproposal );
//       newsumlogrho += logrhoproposal;
//       newsumlogrho += rhoproposal[i];;
//       //compute log ratio of prior densities on log scale
//       LogPriorRatio += rhoalpha* (logrhoproposal-logrho) - 
// 	  (rhobeta0+rhoalpha)*(log(rhobeta1 + rhoproposal[i]) - log(rhobeta1 + rho[i])) ; 
//   }

//   //get log likelihood at current parameter values, annealed if this is an annealing run
//   for(int i = 0; i < IC->getSize(); ++i) {
//     Individual* ind = IC->getIndividual(i);
//     ind->HMMIsBad(true);//to force HMM update
//     LogLikelihood += ind->getLogLikelihood(options, false, true); // don't force update, store result if updated
//     ind->HMMIsBad(true); // HMM probs overwritten by next indiv, but stored loglikelihood still ok
//   }
//   // set ancestry correlations using proposed value of sum-intensities
//   // value for X chromosome set to half the autosomal value 
//   Loci->SetLocusCorrelation(rhoproposal);
  
//   //get log HMM likelihood at proposal rho and current admixture proportions
//   for(int i = 0; i < IC->getSize(); ++i) {
//     Individual* ind = IC->getIndividual(i);
//     LogLikelihoodAtProposal += ind->getLogLikelihood(options, true, false); // force update, do not store result 
//     ind->HMMIsBad(true); // set HMM probs as bad but stored log-likelihood is still ok
//     // line above should not be needed for a forced update with result not stored
//   }
//   LogLikelihoodRatio = LogLikelihoodAtProposal - LogLikelihood;

//   // getGammaLogDensity(rhoalpha, rhobeta, rhoprop) - getGammaLogDensity(rhoalpha, rhobeta, rho[0]);
//   LogAccProbRatio = LogLikelihoodRatio + LogPriorRatio; 
//   //cout << "logLratio = " << LogLikelihoodRatio << " logPriorRatio = " << LogPriorRatio;
  
//   // generic Metropolis step
//   if( LogAccProbRatio < 0 ) {
//     if( log(Rand::myrand()) < LogAccProbRatio ) accept = true;
//   } else accept = true;  
  
//   if(accept) {
//     //cout << " accept";
//     copy(rhoproposal.begin(), rhoproposal.end(), rho.begin());
//     for(int i = 0; i < IC->getSize(); ++i){
//       Individual* ind = IC->getIndividual(i);
//       ind->storeLogLikelihood(false); // store log-likelihoods calculated at rhoprop, but do not set HMM probs as OK 
//     }
//     RhoPriorArgs.sumlogrho = newsumlogrho;
//     RhoPriorArgs.sumrho = newsumrho;
//   } else { 
//     // restore ancestry correlations in Chromosomes using original value of sum-intensities
//     Loci->SetLocusCorrelation(rho);
//   } // stored loglikelihoods are still ok
//   //cout << endl;
//     //update sampler object every w updates
//   if( !( NumberOfUpdates % w ) ){
//     //update sampler object every w updates
//     if( !( NumberOfUpdates % w ) ){
// 	step = TuneRhoSampler.UpdateStepSize( exp(LogAccProbRatio) );  
//     }
//   }
//   //accumulate sums of log of rho
//   if(sumlogrho)
//     for(unsigned i = 0; i < rho.size(); ++i)
//       SumLogRho[i] += log(rho[i]);

// //   if(Comms::isMaster()){
// //       SampleHapmixRhoPriorParameters();
// //   }
// }

///conjugate update of locus-specific sumintensities, conditional on observed numbers of arrivals
// void Latent::SampleSumIntensities(const vector<unsigned> &SumNumArrivals, unsigned NumIndividuals, 
// 				  bool sumlogrho) {
//   double sum = 0.0;
//   int locus = 0;
//   for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
//     ++locus;//skip first locus on each chromosome
//     Chromosome* C = Loci->getChromosome(c); 
//     for(unsigned i = 1; i < C->GetSize(); ++i){
// 	double EffectiveL = Loci->GetDistance(locus) * 2 * NumIndividuals;//length of interval * # gametes
// 	rho[locus] = Rand::gengam( rhoalpha + (double)(SumNumArrivals[locus]), rhobeta + EffectiveL );
// 	sum += rho[locus];
// 	++locus;
//     }
//     //set locus correlation
//     C->SetLocusCorrelation(rho, false, false);
//     C->SetStateArrivalProbs(globaltheta, options->isRandomMatingModel(), true);
//   }
//   //sample rate parameter of gamma prior on rho
//   rhobeta = Rand::gengam( rhoalpha * (double)(rho.size()-1) + rhobeta0, sum + rhobeta1 );
  
//   //accumulate sums of log of rho
//   if(sumlogrho)
//     transform(rho.begin(), rho.end(), SumLogRho.begin(), SumLogRho.begin(), std::plus<double>());
// }

///samples locus-specific sumintensities in a hapmixmodel, using Hamiltonian Monte Carlo
void Latent::SampleSumIntensities(const int* SumAncestry, bool sumlogrho){
  //SumAncestry is a (NumberOfCompositeLoci) * 2 array of counts of gametes with ancestry states that are
  // unequal (element 0) and equal 
  //sumlogrho indicates whether to accumulate sums of log rho

  //RhoArgs.theta = globaltheta;
  //double sum = 0.0;
  int locus = 0;
  int index = 0;
  
  if(Comms::isMaster()){
#ifdef PARALLEL
    MPE_Log_event(9, 0, "sampleRho");
#endif
    RhoPriorArgs.sumrho = 0.0;
    RhoPriorArgs.sumlogrho = 0.0;
    try{
      vector<double>::iterator rho_iter = rho.begin();
      vector<double>::iterator sumlogrho_iter = SumLogRho.begin();
      for(unsigned c = 0; c < Loci->GetNumberOfChromosomes(); ++c){
	++locus;//skip first locus on each chromosome
	for(unsigned i = 1; i < Loci->GetSizeOfChromosome(c); ++i){
	  //for(vector<double>::iterator rho_iter = rho.begin(); rho_iter != rho.end(); ++rho_iter){
	  double logrho = log(*rho_iter);//sampler is on log scale
	  
	  RhoArgs.Distance = Loci->GetDistance(locus);//distance between this locus and last
	  RhoArgs.SumAncestry = SumAncestry + locus*2;//sums of ancestry for this locus
	  
	  //sample new value
	  RhoSampler[index++].Sample(&logrho, &RhoArgs);
	  *rho_iter = exp(logrho);

	  //accumulate sums of log of rho
	  if(sumlogrho)
	    *(sumlogrho_iter++) += logrho;

	  RhoPriorArgs.sumrho += *rho_iter;
	  RhoPriorArgs.sumlogrho += logrho;

	  //accumulate sums, used to sample rhobeta
	  //sum += rho[locus];
	  ++locus;
	  ++rho_iter;

	}
      }
    }
    catch(string s){
      stringstream err;err << "Error encountered while sampling sumintensities " << locus << ":\n" + s;
      throw(err.str());
    }
#ifdef PARALLEL
    MPE_Log_event(10, 0, "sampledRho");
#endif
  }
//broadcast rho to all processes, in Comm (excludes freqsampler)
//no need to broadcast globaltheta if it is kept fixed
#ifdef PARALLEL
  MPE_Log_event(17, 0, "Bcastrho");
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
//   if(Comms::isMaster()){
//     SampleHapmixRhoPriorParameters();
//   }
}

void Latent::SampleHapmixRhoPriorParameters(){
//   //random walk update of rhoalpha
//   const double shape = 4000.0;
//   const double rate = 1.0;

//   double log_proposal = Rand::gennor(log(rhoalpha), RhoAlphaSampler.getStepSize());
//   double proposal = exp(log_proposal);
//   double LogLikelihoodRatio = (rho.size()*log(rhobeta) + RhoPriorArgs.sumlogrho)*(proposal - rhoalpha) + lngamma(rhoalpha) - lngamma(proposal);
//   double LogPriorRatio = shape*(log_proposal - log(rhoalpha)) - rate*(proposal - rhoalpha);

//   if(RhoAlphaSampler.Accept_Reject(LogLikelihoodRatio + LogPriorRatio)){
//     rhoalpha = proposal;
//     RhoArgs.rhoalpha = rhoalpha;
//   }

  //conjugate update of rhobeta
  rhobeta = Rand::gengam(rhobeta0 + rhoalpha*rho.size(), rhobeta1 + RhoPriorArgs.sumrho);
  
  RhoArgs.rhobeta = rhobeta;
}

double Latent::logLikelihoodRhoPriorParams(double shape, double rate, double sumrho, double sumlogrho, double L){
  return(L *(shape*log(rate) + lngamma(shape)) + shape*sumlogrho - rate*sumrho );
}

///energy function for sampling locus-specific sumintensities 
///conditional on locus ancestry and gamma-gamma prior params 
double Latent::RhoEnergy(const double* const x, const void* const vargs){
  //x is log rho
  const RhoArguments* args = (const RhoArguments*)vargs;
  unsigned K = args->NumPops;
  const double d = args->Distance;
  double theta = 1.0 / (double)K;
  double E = 0.0;
  try {
    double rho = myexp(*x);
    double f = myexp(-d*rho);
    int sumequal = args->SumAncestry[1], sumnotequal = args->SumAncestry[0];
    double probequal = theta + f*(1.0 - theta);
    E -= sumnotequal * log(1.0-probequal);//constant term in log(1-theta) omitted
    E -= sumequal * log(probequal);

    // log unnormalized gamma prior on log rho
    E -= args->rhoalpha * *x - args->rhobeta * rho;
    //log unnormalized gamma-gamma prior
    //E -= args->rhoalpha* *x - (args->rhoalpha + args->rhobeta)*log(args->rhobeta1 + rho);
  } catch(string s){
    throw string("Error in RhoEnergy: " + s);
  }
  return E; 
}

///gradient function for sampling locus-specific sumintensities
///conditional on locus ancestry and gamma-gamma prior params
void Latent::RhoGradient( const double* const x, const void* const vargs, double* g ){
  const RhoArguments* args = (const RhoArguments*)vargs;
  unsigned K = args->NumPops;
  const double d = args->Distance;
  double theta = 1.0 / (double)K;
  g[0] = 0.0;
  try {
    double rho = myexp(*x);
    double f = myexp(-d*rho);
    //first compute dE / df
    int sumequal = args->SumAncestry[1], sumnotequal = args->SumAncestry[0];
    double probequal = theta + f*(1.0 - theta);

    g[0] +=  sumnotequal / (1.0 - probequal) - sumequal / probequal;//dprobequal / df
    g[0] *= -rho*d*f*(1.0-theta); // dprobequal / df * df / dx 
    
    //derivative of minus log prior wrt log rho
    g[0] -= args->rhoalpha  - args->rhobeta * rho;

    //g[0] -= args->rhoalpha - rho*(args->rhoalpha + args->rhobeta)/(args->rhobeta1 + rho);

  } catch(string s) {
    throw string("Error in RhoGradient: " + s);
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
      outputstream << "SumIntensities.Mean\tSumIntensities.Variance\tRhoalpha\tRhobeta0\tRhobeta1";
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

    for(vector<double>::const_iterator rho_iter = rho.begin(); rho_iter != rho.end(); ++rho_iter){
	sum += *rho_iter;
	var += (*rho_iter) * (*rho_iter);
    }
    double size = (double)(Loci->GetNumberOfCompositeLoci()-Loci->GetNumberOfChromosomes());
    var = var - (sum*sum) / size;
    (*out) << setiosflags(ios::fixed) << setprecision(6) << sum / size << "\t" << var /size << "\t"
	   << rhoalpha << "\t" << rhobeta << "\t";
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
//random-walk sampler
//     Log << "Expected acceptance rate in sumintensities sampler: "
// 	<< 	TuneRhoSampler.getExpectedAcceptanceRate();
//     Log << "\nwith final step size of "
// 	<< 	TuneRhoSampler.getStepSize()
// 	<< "\n";


//Hamiltonian sampler
    double av = 0;//average acceptance rate
    for(unsigned j = 0; j < rho.size(); ++j){
	av += RhoSampler[j].getAcceptanceRate();
      //cout << j << " " << RhoSampler[j].getAcceptanceRate() << endl;
    }
    Log << "Average Expected acceptance rate in sumintensities samplers: "
	<< av / (double)(Loci->GetNumberOfCompositeLoci()-Loci->GetNumberOfChromosomes());

    av = 0;//average stepsize
    for(unsigned j = 0; j < rho.size(); ++j){
      av += RhoSampler[j].getStepsize();
      if(isnan(RhoSampler[j].getStepsize()))cout << j << " " << RhoSampler[j].getStepsize() << " " << RhoSampler[j].getAcceptanceRate() << endl;  
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
