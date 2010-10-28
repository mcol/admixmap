/*
 *   ADMIXMAP
 *   PopAdmix.cc
 *   Class to hold and update population admixture and sumintensities parameters and their priors
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *   Portions Copyright (c) 2009-2010 Marco Colombo
 *
 * This program is free software distributed WITHOUT ANY WARRANTY.
 * You can redistribute it and/or modify it under the terms of the GNU General Public License,
 * version 2 or later, as published by the Free Software Foundation.
 * See the file COPYING for details.
 */

//=============================================================================
/// \file PopAdmix.cc
/// Implementation of the PopAdmix class.
//=============================================================================

#include "PopAdmix.h"
#include "Chromosome.h"
#include "AdmixIndividualCollection.h"
#include <algorithm>


using namespace std;
using namespace genepi;


// Print the details of the update of the sum intensities parameter rho
#define DEBUG_UPDATE_GLOBAL_SUM_INTENSITIES 0

// Print the details of the update of the odds ratios vector psi
#define DEBUG_UPDATE_ODDS_RATIOS 0


PopAdmix::PopAdmix( const AdmixOptions& op, Genome& loci) :
	options ( op			   ),
	Loci	( loci			   ),
	K	( op.getPopulations()	   ),
	poptheta( op.getPopulations(), 0.0 )
    {
    rho.push_back(0.0);
    }


void PopAdmix::Initialise(int NumObservations, const Vector_s& PopulationLabels,
                          bclib::LogWriter& Log) {

  Log.setDisplayMode(bclib::On);

  // ** Initialise population admixture distribution Dirichlet parameters alpha **
  alpha = options.getInitAlpha();
  SumAlpha.resize( K );
  if(!options.getIndAdmixHierIndicator())  copy(alpha[0].begin(), alpha[0].end(), SumAlpha.begin());

  if(K > 1){
    // ** set up sampler for alpha **
    //set values for sampler
    const genepi::cvector<float> & samplerparams = options.getPopAdmixSamplerParams();
    const size_t size = samplerparams.size();
    float initial_stepsize = size? samplerparams[0] : 0.02;
    unsigned num_leapfrogs = size? (unsigned)samplerparams[1] : 40;
    PopAdmixSampler.SetSize(NumObservations, K,
                            initial_stepsize, num_leapfrogs);

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
	const genepi::cvector<float> & rhosamplerparams = options.getrhoSamplerParams();
	const size_t size = rhosamplerparams.size();
	float initial_stepsize = size? rhosamplerparams[0] : step0;
	float min_stepsize = size? rhosamplerparams[1] : 0.01;
	float max_stepsize = size? rhosamplerparams[2] : 10;
	float target_acceptrate = size? rhosamplerparams[3] : 0.44;
	TuneRhoSampler.SetParameters( initial_stepsize, min_stepsize, max_stepsize, target_acceptrate);
      }
    }

    if (Loci.isX_data()) {

      NumberOfPsiUpdates = 0;
      psi.resize(K, 1.0);
      SumLogPsi.resize(K);
      w = 1;

      // prior parameters on log psi
      psimean0 = options.getLogPsiPriorMean();
      psiprec0 = options.getLogPsiPriorPrec();

      // global psi (one odds ratio vector)
      if (options.isGlobalPsi()) {
        psistep.resize(K, 1.0);

        // initialise the tuner for psi: we create a stepsize tuner for each
        // element of psi, but we never use the first element as it corresponds
        // to psi[0] which is 1.0 by definition
        TunePsiSampler.resize(K);
        for (int i = 1; i < K; ++i)
          TunePsiSampler[i].SetParameters(1.0, 0.01, 10, 0.44);
      }

      // hierarchical model on psi (one odds ratio vector per individual)
      else {
        psimu.resize(K, psimean0);
        psitau.resize(K, psiprec0);

        // gamma prior parameters on the precision
        psialpha0 = options.getLogPsigammaShape();
        psibeta0  = options.getLogPsigammaRate();
      }
    }

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
void PopAdmix::UpdatePopAdmixParams(int iteration,
                                    const AdmixIndividualCollection* const individuals,
                                    bclib::LogWriter &Log) {
  const int Populations = options.getPopulations();

  if ( Populations > 1 && individuals->getSize() > 1 &&
       options.getIndAdmixHierIndicator() ) {
     const double* sumlogtheta = individuals->getSumLogTheta();

#if 0
     cout << "alpha " << alpha[0][0] << " " << alpha[0][1]
          << " sumlogtheta " << sumlogtheta[0] << " " << sumlogtheta[1] << endl;
#endif

     // sample alpha conditional on individual admixture proportions
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

   if ( iteration == options.getBurnIn() && Populations > 1 ) {
     if(options.getNumberOfOutcomes() > 0){
       Log << bclib::Off << "Individual admixture centred in regression model around: ";
       for (int i = 0; i < Populations; ++i)
         Log << poptheta[i] << "\t";
       Log << "\n";
     }
     fill(SumAlpha.begin(), SumAlpha.end(), 0.0);
   }

   if ( iteration < options.getBurnIn() && Populations > 1 ) {
     // accumulate ergodic average of population admixture, which is used to centre
     // the values of individual admixture in the regression model
     double sum = accumulate(SumAlpha.begin(), SumAlpha.end(), 0.0);
     if ( options.getNumberOfOutcomes() > 0 )
       for (int j = 0; j < Populations; ++j)
         poptheta[j] = SumAlpha[j] / sum;
   }
}


//-----------------------------------------------------------------------------
/// Updates global sumintensities in a globalrho model, using random-walk Metropolis-Hastings.
//-----------------------------------------------------------------------------

void PopAdmix::UpdateGlobalSumIntensities( const AdmixIndividualCollection & IC, bool sumlogrho ) {

  using bclib::Rand;

  const int IC_size = IC.getSize();

  if( options.isGlobalRho() ) {

    double LogLikelihood = 0.0;
    double LogLikelihoodAtProposal = 0.0;
    double logrho0 = log( rho[0] );
    double logrhoprop = Rand::gennor( logrho0, step );
    double rhoprop = exp( logrhoprop ); // propose log rho from normal distribution with SD step

    NumberOfUpdates++;


    // Get log likelihood at current parameter values, annealed if this is an annealing run
    #if defined(_OPENMP) && PARALLELIZE_PEDIGREE_LOOP
      #pragma omp parallel for reduction(+:LogLikelihood) default(shared) PED_LOOP_OMP_SCHED if(options.getUsePedForInd())
    #endif
    for( int i = 0 ; i < IC_size ; i++ ) {

      PedBase & ind = IC.getElement(i);

      ind.HMMIsBad(true);//to force HMM update
      LogLikelihood += ind.getLogLikelihood(options, false, true); // don't force update, store result if updated
      ind.HMMIsBad(true); // HMM probs overwritten by next indiv, but stored loglikelihood still ok

    }


    // set ancestry correlations using proposed value of sum-intensities
    // value for X chromosome set to half the autosomal value
    Loci.SetLocusCorrelation( rhoprop );	  // For individuals
    for ( int idx = IC_size ; idx-- != 0 ; ) // For pedigrees
	{
	IC.getElement( idx ).startRhoProposal();
	IC.getElement( idx ).setRho( rhoprop );
	}


    // Get log HMM likelihood at proposal rho and current admixture proportions
    #if defined(_OPENMP) && PARALLELIZE_PEDIGREE_LOOP
      #pragma omp parallel for reduction(+:LogLikelihoodAtProposal) default(shared) PED_LOOP_OMP_SCHED if(options.getUsePedForInd())
    #endif
    for (int i = 0; i < IC_size; i++) {
      PedBase & ind = IC.getElement(i);

      LogLikelihoodAtProposal += ind.getLogLikelihood( options, true, false ); // force update, do not store result
      ind.HMMIsBad(true); // set HMM probs as bad but stored log-likelihood is still ok
      // line above should not be needed for a forced update with result not stored
    }

    const double LogLikelihoodRatio = LogLikelihoodAtProposal - LogLikelihood;

    //compute log ratio of prior densities in log rho basis
    double LogPriorRatio = rhoalpha * (logrhoprop - logrho0) - rhobeta * (rhoprop - rho[0]);
    // getGammaLogDensity(rhoalpha, rhobeta, rhoprop) - getGammaLogDensity(rhoalpha, rhobeta, rho[0]);
    const double LogAccProbRatio = LogLikelihoodRatio + LogPriorRatio;

    // generic Metropolis step
    const bool accept = (LogAccProbRatio >= 0) || (log(Rand::myrand()) < LogAccProbRatio);

#if DEBUG_UPDATE_GLOBAL_SUM_INTENSITIES
    fprintf(stderr, "Current rho: %.9lf  Proposed: %.9lf  Step: %.9lf\n",
            rho[0], rhoprop, step);
    fprintf(stderr, "Curr-LL: %.6lf  Prop-LL: %.6lf"
                    "  Ratio: %.6lf  PriorRat: %.6lf\n",
            LogLikelihood, LogLikelihoodAtProposal,
            LogLikelihoodRatio, LogPriorRatio);
    fprintf(stderr, "Accept: %s\n", accept ? "yes" : "no");
#endif

    if ( accept ) {
      rho[0] = rhoprop;
      logrho0 = logrhoprop;
    }

    // Update sampler object every w updates
    if( !( NumberOfUpdates % w ) ){
      step = TuneRhoSampler.UpdateStepSize( exp(LogAccProbRatio) );
    }

    // Accumulate sum of log of sumintensities after burnin.
    if ( sumlogrho ) SumLogRho[0] += logrho0;

    if ( accept ) {
      for (int i = 0; i < IC_size; i++) {
	// Old-style:
	IC.getElement(i).storeLogLikelihood(false); // store log-likelihoods calculated at rhoprop, but do not set HMM probs as OK
	// New-style:
	//IC.getElement(i).acceptRhoProposal();
      }
    } else {
      // restore ancestry correlations in Chromosomes using original value of sum-intensities
      Loci.SetLocusCorrelation( rho );		    // For individuals
      for ( int idx = IC_size ; idx-- != 0 ; ) // For pedigrees
	  IC.getElement( idx ).rejectRhoProposal(); // was: IC.getElement(idx).setRho( rho[0] );
    } // stored loglikelihoods are still ok
  }//end if global rho model

  else { //individual- or gamete-specific rho model
    if (IC_size > 1 && options.getIndAdmixHierIndicator() ) { // >1 individual and hierarchical model
      double sumrho = IC.GetSumrho();

      // update scale parameter of gamma distribution of sumintensities in population
      if( options.isRandomMatingModel() )
	rhobeta = Rand::gengam( 2*rhoalpha * IC_size + rhobeta0, sumrho + rhobeta1 );
      else
	rhobeta = Rand::gengam( rhoalpha * IC_size + rhobeta0, sumrho + rhobeta1 );

    } // otherwise do not update rhobeta

    // accumulate sum of log of mean of sumintensities after burnin.
    if(sumlogrho)
      SumLogRho[0] += log(rhoalpha) - log(rhobeta);
  }
}


/// Updates the population-level vector psi of odds ratio female/male between
/// ancestral populations.
/// Depending on the value of the "globalpsi" option, we either use a
/// random-walk Metropolis-Hastings algorithm on the population-level parameter
/// (globalpsi=1), or we use a hierarchical model in which we accumulate the
/// individual parameters using a conjugate normal-gamma prior (globalpsi=0).
void PopAdmix::UpdateOddsRatios(const AdmixIndividualCollection& IC,
                                bool sumlogpsi) {

  using bclib::Rand;

  const int IC_size = IC.getSize();

  if (sumlogpsi)
    NumberOfPsiUpdates++;

  // Individual-level odds ratios
  if (!options.isGlobalPsi()) {

    // sample the individual-level parameters
    for (int i = 0; i < IC_size; ++i)
      IC.getElement(i).SamplePsi(options, psimu, psitau, sumlogpsi);

    // Update the population-level parameters by using the conjugate prior:
    // we skip element 0 as it is 1.0 by definition
    for (int el = 1; el < K; ++el) {

      // accumulate the sufficient statistics
      double samplemean = 0.0, sumlogs2 = 0.0;
      for (int i = 0; i < IC_size; ++i) {
        double logpsi = log(IC.getElement(i).getPsi()[el]);
        samplemean += logpsi;
        sumlogs2 += logpsi * logpsi;
      }
      samplemean /= IC_size;

      // posterior hyperparameters
      double prec = psiprec0 + IC_size;
      double mean = (psiprec0 * psimean0 + IC_size * samplemean) / prec;

      double alpha = psialpha0 + 0.5 * IC_size;
      double sumsqdevn = sumlogs2 - IC_size * samplemean * samplemean;
      double otherterm = IC_size * psiprec0 / prec
                       * (samplemean - psimean0) * (samplemean - psimean0);
      double beta = psibeta0 + 0.5 * (sumsqdevn + otherterm);

      // sample the population mean from a gaussian
      psimu[el] = Rand::gennor(mean, 1/sqrt(prec));

      // sample the population precision from a gamma
      psitau[el] = Rand::gengam(alpha, beta);

      // accumulate sum of log of psi after burnin
      if (sumlogpsi)
        SumLogPsi[el] += psimu[el];

      // store the current value of psi
      psi[el] = exp(psimu[el]);
    }

    return;
  }

  // Population-level odds ratio
  // In the update of the odds ratio vector psi, we skip element 0 as it
  // is 1.0 by definition, and we do a random-walk on the other elements
  for (int el = 1; el < K; ++el) {

    double LogLikelihood = 0.0;
    double LogLikelihoodAtProposal = 0.0;

    // propose log psi from normal distribution with SD step
    const double logpsi = log(psi[el]);
    const double logpsiprop = Rand::gennor(logpsi, psistep[el]);
    const double psiprop = exp(logpsiprop);

    // get log likelihood at current parameter values
    for (int i = 0; i < IC_size; ++i) {

      PedBase& ind = IC.getElement(i);

      // force update, don't store the result
      LogLikelihood += ind.getLogLikelihoodXChr(options, true, false);

      // HMM probs overwritten by next indiv
      ind.HMMIsBad(true);
    }

    // store the current psi and set the proposed value
    const double storepsi = psi[el];
    psi[el] = psiprop;

    // get log likelihood at proposed values
    for (int i = 0; i < IC_size; ++i) {

      PedBase& ind = IC.getElement(i);

      ind.startPsiProposal();

      // set the proposed values for psi
      ind.setPsi(psi);

      // force update, don't store the result
      LogLikelihoodAtProposal += ind.getLogLikelihoodXChr(options, true, false);

      // HMM probs overwritten by next indiv
      ind.HMMIsBad(true);
    }

    const double LogLikelihoodRatio = LogLikelihoodAtProposal - LogLikelihood;

    // gaussian prior: -0.5 * tau * (logpsi - mu)^2
    const double LogPriorRatio = -0.5 * psiprec0 * (logpsiprop - logpsi)
                                      * (logpsiprop + logpsi - 2 * psimean0);
    const double LogAccProbRatio = LogLikelihoodRatio + LogPriorRatio;

    // generic Metropolis step
    const bool accept = (LogAccProbRatio >= 0) ||
                        (log(Rand::myrand()) < LogAccProbRatio);

#if DEBUG_UPDATE_ODDS_RATIOS
    fprintf(stderr, "Current psi[%d]: %.9lf  Proposed: %.9lf\n",
            el, storepsi, psiprop);
    fprintf(stderr, "Curr-LL: %.6lf  Prop-LL: %.6lf"
                    "  Ratio: %.6lf  PriorRat: %.6lf\n",
            LogLikelihood, LogLikelihoodAtProposal,
            LogLikelihoodRatio, LogPriorRatio);
    fprintf(stderr, "Accept-Psi: %s\n", accept ? "yes" : "no");
#endif

    if (!accept) {
      psi[el] = storepsi;
    }

    for (int i = 0; i < IC_size; ++i) {
      IC.getElement(i).setPsi(psi);
      if (accept)
        IC.getElement(i).acceptPsiProposal();
      else
        IC.getElement(i).rejectPsiProposal();
    }

    // update sampler object every w updates
    if( !(NumberOfPsiUpdates % w) )
      psistep[el] = TunePsiSampler[el].UpdateStepSize(exp(LogAccProbRatio));

    // Accumulate sum of log psi after burnin
    if (sumlogpsi)
      SumLogPsi[el] += logpsi;
  }
}

/// Computes the posterior mean for the odds ratio vector and stores it in
/// the Individual.
/// This must be called before computing any other posterior means, otherwise
/// we end up using the psi set at the last iteration.
void PopAdmix::StoreOddsRatiosPosteriorMean(const AdmixIndividualCollection& IC) {

  // compute the posterior mean
  for (int i = 0; i < K; ++i)
    psi[i] = exp(SumLogPsi[i] / NumberOfPsiUpdates);

  // store it in all individuals
  for (int i = 0; i < IC.getSize(); ++i)
    IC.getElement(i).setPsi(psi);
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

  // Odds ratios for the X chromosome
  if (Loci.isX_data()) {
    for (int i = 0; i < options.getPopulations(); i++) {
      outputstream << "Psi." << PopulationLabels[i] << "\t";
    }
  }

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

  int prec = out.precision(6);  

  //pop admixture params
  for( int j = 0; j < options.getPopulations(); j++ ){
    out << alpha[0][j];
  }
  //sumintensities
  if( options.isGlobalRho() )
    out << rho[0] ;
  else
    out << rhoalpha / rhobeta ;

  // odds ratios
  if (Loci.isX_data()) {
    out.precision(3);
    for (size_t i = 0; i < psi.size(); ++i) {
      out << psi[i];
      if (!options.isGlobalPsi() && i > 0)
        out << "(" << psitau[i] << ")";
    }
  }

  // restore the original precision
  out.precision(prec);
}

void PopAdmix::OutputParams(){
  //Output to paramfile
  OutputParams(outputstream);
  outputstream << bclib::newline;
}



void PopAdmix::printAcceptanceRates(bclib::LogWriter& Log,
                                    const Vector_s& PopulationLabels) {

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

    if (Loci.isX_data()) {
      if (options.isGlobalPsi()) {
        Log << "Expected acceptance rate in the odds ratios sampler:\n";
        for (int i = 1; i < K; ++i)
          Log << PopulationLabels[i] << ": "
              << TunePsiSampler[i].getExpectedAcceptanceRate()
              << " with final step size of " << psistep[i] << "\n";
      }
      Log << "Odds ratios female/male between ancestral populations:\n";
      for (int i = 1; i < K; ++i)
        Log << PopulationLabels[i] << " vs " << PopulationLabels[0] << ": "
            << exp(SumLogPsi[i] / NumberOfPsiUpdates) << "\n";
    }
}
