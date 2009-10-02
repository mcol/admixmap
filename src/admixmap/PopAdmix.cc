/**
 *   ADMIXMAP
 *   PopAdmix.cc
 *   Class to hold and update population admixture and sumintensities parameters and their priors
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *
 * This program is free software distributed WITHOUT ANY WARRANTY.
 * You can redistribute it and/or modify it under the terms of the GNU General Public License,
 * version 2 or later, as published by the Free Software Foundation.
 * See the file COPYING for details.
 */
#include "PopAdmix.h"
#include "Chromosome.h"
#include "AdmixIndividualCollection.h"
#include "bclib/misc.h"
#include "bclib/dist.h"//for log gamma density
#include <algorithm>
#include <numeric>
#include "gsl/gsl_math.h"
#include "gsl/gsl_specfunc.h"

#define DEBUG_OMP		0
#define DEBUG_RHO_PROPOSAL	0

#if DEBUG_OMP
    #include <omp.h>
    #include "CodeTimer.h"
#endif

using namespace std;
using namespace genepi;

#define PR(x) cerr << #x << " = " << x << endl;

PopAdmix::PopAdmix( const AdmixOptions& op, Genome& loci) :
	options ( op			   ),
	Loci	( loci			   ),
	K	( op.getPopulations()	   ),
	poptheta( op.getPopulations(), 0.0 )
    {
    rho.push_back(0.0);
    }

void PopAdmix::Initialise(int Numindividuals, const Vector_s& PopulationLabels, bclib::LogWriter &Log){
  Log.setDisplayMode(bclib::On);

  // ** Initialise population admixture distribution Dirichlet parameters alpha **
  alpha = options.getInitAlpha();
  SumAlpha.resize( K );
  if(!options.getIndAdmixHierIndicator())  copy(alpha[0].begin(), alpha[0].end(), SumAlpha.begin());

  if(K > 1){
    // ** set up sampler for alpha **
    // should be able to pass initial step size to the sampler
    unsigned obs = Numindividuals;
    if( options.isRandomMatingModel() ){
      obs *= 2;//for 2 gametes per individual
    }
    //set values for sampler
    const genepi::cvector<float> & samplerparams = options.getPopAdmixSamplerParams();
    const size_t size = samplerparams.size();
    float initial_stepsize = size? samplerparams[0] : 0.03;
    unsigned num_leapfrogs = size? (unsigned)samplerparams[1] : 40;
    PopAdmixSampler.SetSize( obs, K, initial_stepsize, num_leapfrogs );

    SumLogRho.push_back(0.0);
    // ** get prior on sum-of-intensities parameter rho or on rate parameter of its population distribution

    rhoalpha = options.getRhoalpha();
    if( (options.getIndAdmixHierIndicator() && !options.isGlobalRho() )){
      // get prior on rate parameter beta and initialize it at prior mean
      rhobeta0 = options.getRhobetaShape();
      rhobeta1 = options.getRhobetaRate();
      rhobeta = rhobeta0 / rhobeta1;
      double initial_rho = rhoalpha * rhobeta1 / (rhobeta0 - 1.0);
      rho[0] = initial_rho;

    }
    else{//global rho or non-hierarchical model on rho
      rhobeta = options.getRhobeta();
      if( options.isGlobalRho()){
	// set up sampler for global variable
	rho[0] = rhoalpha / rhobeta ;//initialise global sumintensities parameter at prior mean for globalrho
	// ** set up TuneRW object for global rho updates **
	NumberOfUpdates = 0;
	w = 1;
	step0 = 1.0; // sd of proposal distribution for log rho
	//need to choose sensible value for this initial RW sd
	step = step0;
	const genepi::cvector<float>& rhosamplerparams = options.getrhoSamplerParams();
	const size_t size = rhosamplerparams.size();
	float initial_stepsize = size? rhosamplerparams[0] : step0;
	float min_stepsize = size? rhosamplerparams[1] : 0.01;
	float max_stepsize = size? rhosamplerparams[2] : 10;
	float target_acceptrate = size? rhosamplerparams[3] : 0.44;
	TuneRhoSampler.SetParameters( initial_stepsize, min_stepsize, max_stepsize, target_acceptrate);
      }
    }


#define DEBUG_RHO 0
#if DEBUG_RHO
  //init = numeric_limits<double>::infinity();
  cerr << "\n"
	"*********************************************\n"
	"********** FORCING POP-ADMIX-RHO TO 0 **********\n"
	"*********************************************\n"
	"\n";
  for ( genepi::RhoType::iterator it = rho.begin() ; it != rho.end() ; ++it )
    *it = 0.0;
#endif

    // ** Open paramfile **
    if ( options.getIndAdmixHierIndicator()){
      Log.setDisplayMode(bclib::Quiet);
      if( strlen( options.getParameterFilename() ) ){
	outputstream.open( options.getParameterFilename());
	if( !outputstream )
	  {
	    Log.setDisplayMode(bclib::On);
	    Log << "ERROR: Couldn't open paramfile\n";
	    exit( 1 );
	  }
	else{
	  Log << "Writing population-level parameters to " << options.getParameterFilename() << "\n";
	  InitializeOutputFile(PopulationLabels);
	}
      }
      else{
	Log << "No paramfile given\n";
      }
    }
  }//end if Populations > 1
}

void PopAdmix::resetStepSizeApproximator(int k) {
    TuneRhoSampler.resetStepSizeApproximator(k);
}


PopAdmix::~PopAdmix()
{
}

/** Samples for population admixture distribution Dirichlet parameters, alpha **
 For a model in which the distribution of individual admixture in the population is a mixture
 of components, we have one Dirichlet parameter vector for each component,
 updated only from those individuals who belong to the component
*/
void PopAdmix::UpdatePopAdmixParams(int iteration, const AdmixIndividualCollection* const individuals, bclib::LogWriter &Log)
 {
   if( options.getPopulations() > 1 && individuals->getSize() > 1 &&
       options.getIndAdmixHierIndicator() ){
     const double* sumlogtheta = individuals->getSumLogTheta();

     //sample alpha conditional on individual admixture proportions
     //cout << "alpha " << alpha[0][0] << " " << alpha[0][1] <<  " sumlogtheta "
     //	    << individuals->getSumLogTheta()[0] << " " <<  individuals->getSumLogTheta()[1] << endl;
     try{
       if(options.PopAdmixturePropsAreEqual())
	 //sample only dispersion
	 PopAdmixSampler.SampleEta(sumlogtheta, alpha[0].getVector_unsafe());
       else//sample proportions and dispersion
	 PopAdmixSampler.Sample( sumlogtheta, alpha[0].getVector_unsafe());
     }
     catch(string s){
       throw string("Error encountered while sampling population admixture parameters:\n" +s);
     }

     copy(alpha[0].begin(), alpha[0].end(), alpha[1].begin()); // alpha[1] = alpha[0]

  }
   // ** accumulate sum of Dirichlet parameter vector over iterations  **
   transform(alpha[0].begin(), alpha[0].end(), SumAlpha.begin(), SumAlpha.begin(), std::plus<double>());//SumAlpha += alpha[0];

   if( iteration == options.getBurnIn() && options.getPopulations() > 1) {
     if(options.getNumberOfOutcomes() > 0){
       Log << bclib::Off << "Individual admixture centred in regression model around: ";
       for(int i = 0; i < options.getPopulations(); ++i)Log << poptheta[i] << "\t";
       Log << "\n";
     }
     fill(SumAlpha.begin(), SumAlpha.end(), 0.0);
   }


   if( iteration < options.getBurnIn() && options.getPopulations() > 1) {
     // accumulate ergodic average of population admixture, which is used to centre
     // the values of individual admixture in the regression model
     double sum = accumulate(SumAlpha.begin(), SumAlpha.end(), 0.0);
     if(options.getNumberOfOutcomes() > 0)for( int j = 0; j < options.getPopulations(); j++ )poptheta[j] = SumAlpha[j] / sum;
   }
}



//-----------------------------------------------------------------------------
///updates global sumintensities in a globalrho model, using random-walk Metropolis-Hastings
//-----------------------------------------------------------------------------

void PopAdmix::UpdateGlobalSumIntensities(const AdmixIndividualCollection* const IC, bool sumlogrho) {
  using bclib::Rand;
  if( options.isGlobalRho() ) {
    double LogLikelihood = 0.0;
    double LogLikelihoodAtProposal = 0.0;
    #if DEBUG_RHO
	double logrho0 = log(0.5);
	double logrhoprop = log(0.5);
	double rhoprop = rho[0];
    #else
	double logrho0 = log( rho[0] );
	double logrhoprop = Rand::gennor( logrho0, step );
	double rhoprop = exp( logrhoprop ); // propose log rho from normal distribution with SD step
    #endif

    NumberOfUpdates++;

    //get log likelihood at current parameter values, annealed if this is an annealing run
    #if defined(_OPENMP) && PARALLELIZE_PEDIGREE_LOOP
      #pragma omp parallel for default(shared) if(options.getUsePedForInd())
    #endif
    for(int i = 0; i < IC->getSize(); i++) {
      PedBase & ind = IC->getElement(i);

      #if DEBUG_OMP
	  genepi::CodeTimer ct;
      #endif

      ind.HMMIsBad(true);//to force HMM update
      LogLikelihood += ind.getLogLikelihood(options, false, true); // don't force update, store result if updated
      ind.HMMIsBad(true); // HMM probs overwritten by next indiv, but stored loglikelihood still ok

      #if DEBUG_OMP
	  #pragma omp critical
	    {
	    if ( i == 0 )
		cout << "DEBUG-OMP-2: " << omp_get_num_threads() << " threads\n";
	    std::cout << "DEBUG-OMP-2: ped#" << i << " (" << dynamic_cast<Pedigree&>(ind).getId() <<
		    ") is thread#" << omp_get_thread_num() << " started at: " << ct.local_started()
		    << "  time: " << ct.local_elapsed() << '\n';
	    }
      #endif
    }

    // set ancestry correlations using proposed value of sum-intensities
    // value for X chromosome set to half the autosomal value
    Loci.SetLocusCorrelation( rhoprop );
    for ( int idx = IC->getSize() ; idx-- != 0 ; )
	IC->getElement(idx).setRho( rhoprop );

    //get log HMM likelihood at proposal rho and current admixture proportions
    #if defined(_OPENMP) && PARALLELIZE_PEDIGREE_LOOP
      #pragma omp parallel for default(shared) if(options.getUsePedForInd())
    #endif
    for(int i = 0; i < IC->getSize(); i++) {
      PedBase & ind = IC->getElement(i);

      #if DEBUG_OMP
	  genepi::CodeTimer ct;
      #endif

      LogLikelihoodAtProposal += ind.getLogLikelihood(options, true, false); // force update, do not store result
      ind.HMMIsBad(true); // set HMM probs as bad but stored log-likelihood is still ok
      // line above should not be needed for a forced update with result not stored

      #if DEBUG_OMP
	  #pragma omp critical
	    {
	    if ( i == 0 )
		cout << "DEBUG-OMP-3: " << omp_get_num_threads() << " threads\n";
	    std::cout << "DEBUG-OMP-3: ped#" << i << " (" << dynamic_cast<Pedigree&>(ind).getId() <<
		    ") is thread#" << omp_get_thread_num() << " started at: " << ct.local_started()
		    << "  time: " << ct.local_elapsed() << '\n';
	    }
      #endif
    }

    const double LogLikelihoodRatio = LogLikelihoodAtProposal - LogLikelihood;

    //compute log ratio of prior densities in log rho basis
    double LogPriorRatio = rhoalpha * (logrhoprop - logrho0) - rhobeta * (rhoprop - rho[0]);
    // getGammaLogDensity(rhoalpha, rhobeta, rhoprop) - getGammaLogDensity(rhoalpha, rhobeta, rho[0]);
    const double LogAccProbRatio = LogLikelihoodRatio + LogPriorRatio;

    // generic Metropolis step
    const bool accept = (LogAccProbRatio >= 0) || (log(Rand::myrand()) < LogAccProbRatio);

#if DEBUG_RHO_PROPOSAL
  fprintf( stderr, "Rho-prop: LogLikelihoodAtProposal=%.9lf; LogLikelihood=%.9lf; LogLikelihoodRatio=%.9lf; LogPriorRatio=%.9lf; LogAccProbRatio=%.9lf; accept=%s\n",
	    LogLikelihoodAtProposal, LogLikelihood, LogLikelihoodRatio, LogPriorRatio, LogAccProbRatio, accept ? "yes" : "no" );
#endif

    if(accept){
      rho[0] = rhoprop;
      logrho0 = logrhoprop;
    }
    //update sampler object every w updates
    if( !( NumberOfUpdates % w ) ){
      step = TuneRhoSampler.UpdateStepSize( exp(LogAccProbRatio) );
    }
    if(sumlogrho )SumLogRho[0] += logrho0;// accumulate sum of log of sumintensities after burnin.

    if(accept) {
      for(int i = 0; i < IC->getSize(); i++){
	IC->getElement(i).storeLogLikelihood(false); // store log-likelihoods calculated at rhoprop, but do not set HMM probs as OK
      }
    } else {
      // restore ancestry correlations in Chromosomes using original value of sum-intensities
      Loci.SetLocusCorrelation(rho);
      for ( int idx = IC->getSize() ; idx-- != 0 ; )
	  IC->getElement(idx).setRho( rho[0] );
    } // stored loglikelihoods are still ok
  }//end if global rho model

  else { //individual- or gamete-specific rho model
    if(IC->getSize()>1 && options.getIndAdmixHierIndicator() ) { // >1 individual and hierarchical model
      double sumrho = IC->GetSumrho();

      // update scale parameter of gamma distribution of sumintensities in population
      if( options.isRandomMatingModel() )
	rhobeta = Rand::gengam( 2*rhoalpha * IC->getSize() + rhobeta0, sumrho + rhobeta1 );
      else
	rhobeta = Rand::gengam( rhoalpha* IC->getSize() + rhobeta0, sumrho + rhobeta1 );

    } // otherwise do not update rhobeta

    // accumulate sum of log of mean of sumintensities after burnin.
    if(sumlogrho)
      SumLogRho[0] += log(rhoalpha) - log(rhobeta);
  }
}



//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void PopAdmix::InitializeOutputFile(const Vector_s& PopulationLabels) {
  // Header line of paramfile
  //Pop. Admixture
  outputstream.delimit(false);
  for( int i = 0; i < options.getPopulations(); i++ ) {
    outputstream << "Dirichlet." << PopulationLabels[i] << "\t";
  }
  //SumIntensities
  if( options.isGlobalRho() ) outputstream << "sumIntensities\t";
  else outputstream << "sumIntensities.mean\t";

  outputstream.delimit(true);
  outputstream << bclib::newline;
}

void PopAdmix::OutputErgodicAvg( int samples, std::ofstream *avgstream)
{
  if(options.getPopulations()>1){
    for( int j = 0; j < options.getPopulations(); j++ ){
      avgstream->width(9);
      *avgstream << setprecision(6) << SumAlpha[j] / samples << "\t";
    }
    avgstream->width(9);
    *avgstream << setprecision(6) << exp(SumLogRho[0] / samples) << "\t";
  }
}

//output to given output stream
void PopAdmix::OutputParams(bclib::Delimitedostream& out){
  //pop admixture params
  for( int j = 0; j < options.getPopulations(); j++ ){
    out << setprecision(6) << alpha[0][ j ];
  }
  //sumintensities
  if( options.isGlobalRho() )
    out << setprecision(6) << rho[0] ;
  else
    out << setprecision(6) << rhoalpha / rhobeta ;
}

void PopAdmix::OutputParams(){
  //Output to paramfile
  OutputParams(outputstream);
  outputstream << bclib::newline;
}

const cvector<double> & PopAdmix::getalpha0() const
    {
    return alpha[0];
    }

const cvector<cvector<double> > &PopAdmix::getalpha()const{
  return alpha;
}

double PopAdmix::getrhoalpha()const{
  return rhoalpha;
}
double PopAdmix::getrhobeta()const{
  return rhobeta;
}
double PopAdmix::getglobalrho()const{
  return rho[0];
}
const genepi::RhoType & PopAdmix::getrho() const {
  return rho;
}
const genepi::RhoType & PopAdmix::getSumLogRho() const {
  return SumLogRho;
}


void PopAdmix::printAcceptanceRates(bclib::LogWriter &Log) {
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

    if( options.isGlobalRho() ){
      Log << "Expected acceptance rate in global sumintensities sampler: "
	  << TuneRhoSampler.getExpectedAcceptanceRate()
	  << "\nwith final step size of "
	  << step	<< "\n";
    }

}
