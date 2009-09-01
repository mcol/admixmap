//=============================================================================
//
// Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
// Portions Copyright (C) 2009 David D. Favro
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 as published by the Free
// Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
//
/// \file PedigreeAdmix.cc
/// Implementation of the admixmap-specific methods of the Pedigree class.
///
/// Additional methods to run the admixmap model.  These should be moved
/// into a derived class, not here in Pedigree.
//
//=============================================================================


//=============================================================================
/// Issues yet to resolve:
///	* Rho indexing.  Is it indexed on for pedigrees?  Which to pass into TPC constructor?
///	*
//=============================================================================

#include "Pedigree.h"


#include <cmath> // exp(), log()

#include <bclib/misc.h>
#include <bclib/rand.h>

#include "TransProbCache.h"
#include "HiddenMarkovModel.new.h"


#if 0
    /// The way it was for individuals.
    #define IS_ADMIXED()  options.isAdmixed( this_gamete )
#else
    /// For pedigrees, we assume that all founders are admixed.
    #define IS_ADMIXED()  true
#endif



//------------------------------------------------------------------------
// NOTE X1:
// These were defined as forward references in the header file:
#include "AdmixOptions.h"
#include "CopyNumberAssocTest.h"
#include "bclib/DataMatrix.h"
#include "Chromosome.h"
#include "AdmixOptions.h"
#include "Individual.h" // For Individual::Loci
//------------------------------------------------------------------------



using bclib::pvector;


#define NOT_YET_PORTED			0

#define NON_GLOBAL_RHO_WORKS		0

#define NOT_NEEDED_FOR_FIRST_VERSION	0

#if ! NOT_NEEDED_FOR_FIRST_VERSION
    #define IAmUnderTest    0
#endif

#define THIS_IS_HOW_IT_WAS  0

#define USE_SINGLE_POPULATION_SPECIAL_CASE  0



namespace genepi { // ----



#if 0
    static bool greater_than_zero( double x )
	{
	return x > 0.0;
	}
#endif



//-----------------------------------------------------------------------------
// getGenome() [static]
// Context accesss.
//-----------------------------------------------------------------------------
inline static Genome & getGenome()
    {
    return Individual::getGenome();
    }



//-----------------------------------------------------------------------------
// getK() [perhaps could be static but currently isn't]
// Nasty hack for context access.
//-----------------------------------------------------------------------------

static PopIdx K = 0;

void Pedigree::setK( PopIdx _K )
    {
    if ( K != 0 )
	throw std::runtime_error( genepi::estr("Set K to ") + _K
				    + " after setting to " + K );

    K = _K;
    }

inline PopIdx Pedigree::getK() const
    {
    gp_assert_ne( K, 0 ); // Assert: has been set already
    return K;
    }



//-----------------------------------------------------------------------------
// getNumChromosomes() [perhaps could be static but currently isn't]
// Context access.
//-----------------------------------------------------------------------------

#if SEPARATE_HMM_FOR_EACH_CHROM
    ChromIdxType Pedigree::getNumChromosomes() const
	{
	#error Not yet implemented.
	}
#endif



//-----------------------------------------------------------------------------
/// This stuff should be in the constructor, when we have a derived class.
//-----------------------------------------------------------------------------

void Pedigree::InitialiseAdmixedStuff( const AdmixOptions & options )
    {
    logLikelihood.value = 0.0;
    logLikelihood.ready = false;
    logLikelihood.HMMisOK = false;
    const double step0 = step;
    ThetaTuner.SetParameters( step0, 0.0001, 10.0, 0.44 );


    #if ! NON_GLOBAL_RHO_WORKS
	if ( ! options.isGlobalRho() )
	    throw std::runtime_error( "Only global-rho currently is supported." );
    #endif


    // This was AdmixedIndividual::InitialiseSumIntensities():
    double init;
    if ( !options.isGlobalRho() && options.getIndAdmixHierIndicator())
	{//model with individual- or gamete-specific sumintensities

	//set prior mean as initial value for rho
	if ( options.getRhobetaShape() > 1 )
	    init = options.getRhoalpha() * options.getRhobetaRate() / (options.getRhobetaShape() - 1 );
	else
	    init = options.getRhoalpha() * options.getRhobetaRate() / options.getRhobetaShape() ;//conditional prior mean
	}
    else // no hierarchical model or globalrho prior
	{
	init = options.getRhoalpha() / options.getRhobeta();
	}


    _rho  .assign( NumGametes, 0.0 ); // set to 0 for unadmixed gametes
    rhohat.assign( NumGametes, 0.0 ); // set to 0 for unadmixed gametes

    for( unsigned g = 0; g < NumGametes; ++g)
	{
	if ( options.isAdmixed(g) )
	    {
	    _rho  [g] = init;
	    rhohat[g] = init;
	    }
	}

    sumlogrho.assign( _rho.size(), 0.0 );
    }



//-----------------------------------------------------------------------------
// HMMIsBad()
// Copied directly from Individual, could be placed in some base class.
//-----------------------------------------------------------------------------

void Pedigree::HMMIsBad( bool loglikisbad )
    {
    logLikelihood.HMMisOK = false;
    if ( loglikisbad )
	logLikelihood.ready = false;
    }



//-----------------------------------------------------------------------------
// getHMM()
//-----------------------------------------------------------------------------

#if ! SEPARATE_HMM_FOR_EACH_CHROM

    TransProbCache & Pedigree::getTPC() const
	{
	#if NON_GLOBAL_RHO_WORKS
	    #error Which rho to use here?
	#else
	    // Assume that _rho[0] is the global rho to use:
	    if ( tpCache == 0 )
		tpCache = new TransProbCache( *this, _rho[0], Theta );
	#endif
	return *tpCache;
	}

    /// Allocate the HiddenMarkovModel for the Pedigree object (only if have one HMM
    /// for each pedigree rather than one per chromosome).

    HiddenMarkovModel & Pedigree::getHMM() const
	{
	if ( hmm == 0 )
	    hmm = new HiddenMarkovModel( *this, getTPC(), &Theta );
	return *hmm;
	}

    void Pedigree::freeHMM() const
	{
	delete tpCache;
	tpCache = 0;
	delete hmm;
	hmm = 0;
	}

#endif



//-----------------------------------------------------------------------------
//
// SampleTheta()
//
/// Samples admixture proportions.  Called with @a RW true for a random-walk
/// proposal, false for a conjugate proposal
//
//-----------------------------------------------------------------------------

void Pedigree::SampleTheta(
		 int				    iteration	      ,
		 double *			    SumLogTheta       ,
		 const bclib::DataMatrix *	    Outcome	      ,
		 const DataType *		    OutcomeType       ,
		 const std::vector<double> &	    lambda	      ,
		 int				    NumCovariates     ,
		 bclib::DataMatrix *		    Covariates	      ,
		 const std::vector<const double*> & beta	      ,
		 const PopThetaType &		    poptheta	      ,
		 const AdmixOptions &		    options	      ,
		 const AlphaType &		    alpha	      ,
		 double				    /*DInvLink*/      ,
		 double				    /*dispersion*/    ,
		 CopyNumberAssocTest &		    /*ancestryAssocTest*/,
		 bool				    RW		      ,
		 bool				    anneal	      )
    {
    const PopIdx K = getK();

    double logpratio = 0.0;
    try
	{
	if ( RW )
	    {
	    NumberOfUpdates++;
	    logpratio += ProposeThetaWithRandomWalk( options, alpha );
	    }
	    #if NOT_NEEDED_FOR_FIRST_VERSION
		else
		  ProposeTheta(options, alpha, SumLocusAncestry, SumLocusAncestry_X);
	    #endif
	}
    catch ( string & s )
	{
	throw std::runtime_error( "Error encountered while generating proposal individual admixture proportions:\n" + s );
	}

    // Calculate Metropolis acceptance probability ratio for proposal theta
    if ( (! options.getTestForAdmixtureAssociation()) && (getMyNumber() < Outcome->nCols()) )
	{
	const int NumOutcomes = Outcome->nCols();
	for ( int k = 0; k < NumOutcomes; k++ )
	    {
	    const RegressionType RegType = (OutcomeType[k] == Binary) ? Logistic : Linear;
	    logpratio += LogAcceptanceRatioForRegressionModel( RegType, options.isRandomMatingModel(), K,
							NumCovariates,
							 Covariates, beta[k], Outcome->get( getIndex(), k ),
							 poptheta, lambda[k]);
	    }
	}

    // Accept or reject proposed value - if conjugate update and no regression
    // model, proposal will be accepted because logpratio = 0
    Accept_Reject_Theta(logpratio, K, options.isRandomMatingModel(), RW );

    // Update the value of admixture proportions used in the regression model, but
    // only if pedigree is a single individual:
    if ( (options.getNumberOfOutcomes() > 0) && (getNMembers() == 1) )
	UpdateAdmixtureForRegression(K, NumCovariates, poptheta, options.isRandomMatingModel(), Covariates);

    if ( (! anneal) && (iteration > options.getBurnIn()) )
	{
	// Accumulate sums in softmax basis for calculation of posterior means

	double a[ K ]; // Or could use: bclib::pvector<double> a;

	for ( ThetaIdx tIdx = 0 ; tIdx < getNTheta() ; ++tIdx )
	    {
	    Theta[tIdx].inv_softmax_gt0( a ); // Equivalent: Theta[tIdx].inv_softmax( a, gt_0 );
	    SumSoftmaxTheta[tIdx] += a;
	    }
	}


    if ( ! IAmUnderTest )
	{
	for ( ThetaIdx tIdx = 0 ; tIdx < getNTheta() ; ++tIdx )
	    {
	    //SumLogTheta.addAfterFunction( Theta[tIdx], log );
	    const pvector<double> & th = Theta[ tIdx ];
	    for ( size_t k = 0 ; k < th.size() ; ++k )
		SumLogTheta[k] += log( th[k] );
	    }

	#if NOT_NEEDED_FOR_FIRST_VERSION
	    //increment B using new Admixture Props
	    //Xcov is a vector of admixture props as covariates as in UpdateScoreForAncestry
	    if(iteration >= options.getBurnIn() && options.getTestForLinkageWithAncestry()){
	      double* admixtureCovars = new double[K-1];
	      for(int t = 0; t < K-1; ++t)admixtureCovars[t] = Covariates->get(getIndex(), Covariates->nCols()-K+1+t);
	       ancestryAssocTest.UpdateB(DInvLink, dispersion, admixtureCovars);
	      delete[] admixtureCovars;
	    }
	#endif
	}

    }



//-----------------------------------------------------------------------------
// LogAcceptanceRatioForRegressionModel()
//
/// Returns log of ratio of likelihoods of new and old values of population
/// admixture in regression models.  individual admixture theta is
/// standardized about the mean poptheta calculated during burn-in.
//-----------------------------------------------------------------------------

double Pedigree::LogAcceptanceRatioForRegressionModel( RegressionType /*RegType*/, bool /*RandomMatingModel*/,
							PopIdx /*K*/, int /*NumCovariates*/,
							const bclib::DataMatrix * /*Covariates*/, const double * /*beta*/,
							const double /*Outcome*/, const PopThetaType & /*poptheta*/,
							double /*lambda*/ ) const
    {

    #if NOT_NEEDED_FOR_FIRST_VERSION

	double logprobratio = 0.0;
	double XBeta	= 0.0;
	double currentXBeta = 0.0;

	vector<double> avgtheta( K );
	avgtheta[0] = 0.0;

	vector<double> currentavgtheta( K );
	currentavgtheta[0] = 0.0;

	// Resume here...
	if( RandomMatingModel && NumGametes==2)
	  for(int k = 1;k < NumHiddenStates; ++k){
	    avgtheta[k] = (ThetaProposal[k] + ThetaProposal[k + NumHiddenStates ])/ 2.0 - poptheta[k];
	    currentavgtheta[k] = (Theta[k] + Theta[k + NumHiddenStates ])/ 2.0 - poptheta[k];
	  }
	else
	  for(int k = 1;k < NumHiddenStates; ++k){
	    avgtheta[k] = ThetaProposal[k]  - poptheta[k];
	    currentavgtheta[k] = Theta[k]  - poptheta[k];
	  }

	for( int jj = 0; jj < NumCovariates - NumHiddenStates + 1; jj++ ){
	  XBeta += Covariates->get( getIndex(), jj ) * beta[jj];
	  currentXBeta += Covariates->get( getIndex(), jj ) * beta[jj];
	}
	for( int k = 1; k < NumHiddenStates; k++ ){
	  XBeta += avgtheta[ k ] * beta[NumCovariates - NumHiddenStates + k ];
	  currentXBeta += currentavgtheta[ k ] * beta[NumCovariates - NumHiddenStates + k];
	}
	if(RegType == Linear) {
	  logprobratio = 0.5 * lambda * (XBeta - currentXBeta) * (Outcome + Outcome - currentXBeta - XBeta);
	  // ( currentXBeta - Outcome )^2 - ( XBeta - Outcome )^2 factorized
	}
	else if(RegType == Logistic) {
	  if( Outcome == 1 ) {
	    logprobratio =  log( ( 1.0 + exp( -currentXBeta ) ) / ( 1.0 + exp( -XBeta ) ) );
	  } else {
	    logprobratio =  log( ( 1.0 + exp( currentXBeta ) ) / ( 1.0 + exp( XBeta ) ) );
	  }
	}
	return( logprobratio );
    #else
	throw std::runtime_error( "Not implemented" );
    #endif
    }



//-----------------------------------------------------------------------------
// UpdateAdmixtureForRegression()
// Update individual admixture values (mean of both gametes) used in the regression model
//-----------------------------------------------------------------------------

void Pedigree::UpdateAdmixtureForRegression( int NumHiddenStates, int NumCovariates,
                                               const PopThetaType & poptheta, bool /*RandomMatingModel*/,
					       bclib::DataMatrix * Covariates )
    {
    cvector<double> avgtheta( NumHiddenStates );


    for ( PopIdx k = getK() ; k-- != 0 ; )
	{
	double sum = 0.0;
	for ( ThetaIdx tIdx = getNTheta() ; tIdx-- != 0 ; )
	    sum += Theta[ tIdx ][ k ];
	avgtheta[ k ] = sum / getNTheta();
	}


    const PopIdx K = getK();
    const size_t NumberOfInputCovariates = NumCovariates - K;
    // Seems strange, why does this loop start at 1?
    for( PopIdx k = 1 , covIdx = NumberOfInputCovariates ; k < K ; k++, covIdx++ )
	Covariates->set( getIndex(), covIdx, avgtheta[ k ] - poptheta[ k ] );
    }



//-----------------------------------------------------------------------------
// ProposeThetaWithRandomWalk()
//-----------------------------------------------------------------------------

double Pedigree::ProposeThetaWithRandomWalk( const AdmixOptions & options, const AlphaType & alpha )
    {
    double LogLikelihoodRatio = 0.0;
    double LogPriorRatio = 0.0;


    //-------------------------------------------------------------------------
    // Generate proposals
    //-------------------------------------------------------------------------

    const PopIdx K = getK();

    pvector<double> a( K ); // Maybe should be at class scope?

    for ( ThetaIdx tIdx = 0 ; tIdx < getNTheta() ; ++tIdx )
	{

	if ( IS_ADMIXED() )
	    {
	    const pvector<double> & th	    = Theta	   [ tIdx ];
	    pvector<double> &	    th_prop = ThetaProposal[ tIdx ];

	    th.inv_softmax_gt0( a );

	    // Random walk step - on all elements of array a
	    for ( PopIdx k = 0 ; k < K ; ++k )
		if ( a[k] != 0.0 )
		    a[k] = bclib::Rand::gennor( a[k], step );

	    // Reverse transformation from numbers on real line to proportions
	    a.softmax_gt0( th_prop );

	    // Compute contribution of this founder to log prior ratio
	    // Prior densities must be evaluated in softmax basis
	    for ( PopIdx k = 0; k < K; ++k )
		LogPriorRatio += alpha[tIdx][k] * (log(th_prop[k]) - log(th[k]));
	    }

	else // IS_ADMIXED()
	    Theta[ tIdx ] = ThetaProposal[ tIdx ];

	} // End loop over Theta's elements (i.e. founders)


    // Get log likelihood at current parameter values - do not force update, store result of update
    LogLikelihoodRatio -= getLogLikelihood(options, false, true);

    // Get log likelihood at proposal theta and current rho - force update
    // store result in loglikelihood.tempvalue, and accumulate loglikelihood ratio
    logLikelihood.tempvalue = getLogLikelihood( options, ThetaProposal, _rho, true );
    LogLikelihoodRatio += logLikelihood.tempvalue;

    return LogLikelihoodRatio + LogPriorRatio;// log ratio of full conditionals
    }



//-----------------------------------------------------------------------------
// Accept_Reject_Theta() [private]
//
/// Metropolis update for admixture proportions theta, taking log of acceptance
/// probability ratio as argument uses log ratio because this is easier than
/// ratio to calculate for linear regression model if conjugate update and no
/// regression model, logpratio remains set to 0, so all proposals are accepted
//
//-----------------------------------------------------------------------------

void Pedigree::Accept_Reject_Theta( double logpratio, int /*NumHiddenStates*/, bool /*RandomMatingModel*/, bool RW )
    {
    // Loop over populations: if any element of proposed Dirichlet parameter
    // vector is too small, reject update without test step
    bool should_do_test = true;
    for ( ThetaIdx tIdx = getNTheta() ; should_do_test && (tIdx-- != 0) ; )
	{
	const pvector<double> & th	= Theta	       [ tIdx ];
	const pvector<double> & th_prop = ThetaProposal[ tIdx ];

	for ( PopIdx k = getK() ; should_do_test && (k-- != 0) ; )
	    if ( (th[k] > 0.0) && (th_prop[k] < 0.0001) )
		should_do_test = false;
	}


    bool accept = false;
    double AccProb = 1.0;

    if ( should_do_test )
	{
	// Generic Metropolis step
	if ( logpratio < 0 )
	    {
	    AccProb = exp( logpratio );
	    if ( bclib::Rand::myrand() < AccProb )
		accept = true;
	    }
	else
	    accept = true;
	}


    if ( accept )
	{
	// Set proposed values as new values
	setAdmixtureProps( ThetaProposal );
	if ( RW )
	    {
	    // If random-walk update, store the temp log-likelihood and set loglikelihood.HMMisOK to true
	    storeLogLikelihood( true );
	    }
	else
	    {
	    // conjugate update of parameters invalidates both HMM forward probs and stored loglikelihood
	    logLikelihood.HMMisOK = false;
	    logLikelihood.ready = false;
	    }
	}
    else
	{
	// If RW proposal is rejected, loglikelihood.HMMisOK is already set to
	// false, and stored log-likelihood is still valid
	}


    if ( RW )
	{
	// Update step size in tuner object every w updates
	if ( (NumberOfUpdates % w) == 0 )
	    step = ThetaTuner.UpdateStepSize( AccProb );
	}

    }



//-----------------------------------------------------------------------------
// storeLogLikelihood()
//
/// To call if a Metropolis proposal is accepted
//-----------------------------------------------------------------------------

void Pedigree::storeLogLikelihood( bool setHMMAsOK )
    {
    logLikelihood.value = logLikelihood.tempvalue;
    logLikelihood.ready = true;
    if ( setHMMAsOK )
	logLikelihood.HMMisOK = true;
    }



//-----------------------------------------------------------------------------
/// Set admixture proportions
//-----------------------------------------------------------------------------

void Pedigree::setAdmixtureProps( const ThetaType & rhs )
    {
    gp_assert( Theta.size() == getNTheta() );
    gp_assert( rhs  .size() == getNTheta() );

    for ( ThetaIdx tIdx = getNTheta() ; tIdx-- != 0 ; )	
{
	pvector<double> &	sub_theta = Theta [ tIdx ];
	const pvector<double> & sub_rhs	  = rhs	  [ tIdx ];

	gp_assert( sub_theta.size() == getK() );
	gp_assert( sub_rhs  .size() == getK() );

	for ( PopIdx k = getK() ; k-- != 0 ; )
	    sub_theta[ k ] = sub_rhs[ k ];
	}
    }



//*****************************************************************************
//***************************  Log-Likelihoods  *******************************
//*****************************************************************************

//-----------------------------------------------------------------------------
// public function:
// Calls private function to get log-likelihood at current parameter values, and
// stores it either as loglikelihood.value or as loglikelihood.tempvalue store
// should be false when calculating energy for an annealed run, or when
// evaluating proposal for global sum-intensities
// **** COPIED FROM AdmixedIndividual, NOT Individual ****
//-----------------------------------------------------------------------------

double Pedigree::getLogLikelihood( const Options& options, const bool forceUpdate, const bool store)
    {
    if ( ! logLikelihood.ready || forceUpdate )
	{
	logLikelihood.tempvalue = getLogLikelihood(options, Theta, _rho, true);
	if ( store )
	    {
	    logLikelihood.value = logLikelihood.tempvalue;
	    logLikelihood.ready = false;    //true;
	    logLikelihood.HMMisOK = false;  //true; //because forward probs now correspond to current parameter values
	    }					    //and call to UpdateHMM has set this to false
	return logLikelihood.tempvalue;
	}
    else
	return logLikelihood.value; // nothing was changed
    }



//-----------------------------------------------------------------------------
//
// Private function: gets log-likelihood at parameter values specified as
// arguments, but does not update loglikelihoodstruct
//
// This is a combination of the AdmixedIndividual:: and Individual:: versions,
// because the AdmixedIndividual version called the superclass's version
// internally.
//
//-----------------------------------------------------------------------------

double Pedigree::getLogLikelihood( const Options & options, const ThetaType & theta,
				    const RhoType & rho,  bool updateHMM) const
    {
    double rv;

    #if USE_SINGLE_POPULATION_SPECIAL_CASE
      if ( getK() == 1 )
	rv = getLogLikelihoodOnePop();
      else
    #endif
	{
	// ==== This was the body of the Individual:: version, a superclass-call
	// ==== from AdmixedIndividual::getLogLikelihood()
	// ==== rv = Individual::getLogLikelihood(options, theta, rho, updateHMM);


	#if SEPARATE_HMM_FOR_EACH_CHROM
	    rv = 0.0;
	    for ( ChromIdxType chromIdx = 0 ; chromIdx < getNumChromosomes() ; chromIdx++ )
		{
		if ( updateHMM )
		    {// force update of forward probs
		    UpdateHMMInputs( chromIdx, options, theta, rho );
		    }

		// For individuals: is_diploid = !isHaploid && (!Loci->isXChromosome(j) || SexIsFemale)
		// Should be: getChromosome(j).isDiploid()
		// But in fact, since we have only one HMM per chromosome, it should
		// remember its state-space, so not need to have is-diploid passed
		// into each method-call.
		rv += getHMM( chromIdx ).getLogLikelihood();
		}
	#else
	    if ( updateHMM )
		updateHMMInputs( options, theta, rho );
	    rv = getHMM().getLogLikelihood();
	#endif

	// FIXME: under what circumstances does this happen?  Throw the proper
	// exception from within HMM instead.
	if ( isnan(rv) )
	    throw std::runtime_error( "HMM returns log-likelihood as nan (not a number)\n" );
	}

    return rv;
    }



//-----------------------------------------------------------------------------
// getLogLikelihoodAtPosteriorMeans()
//-----------------------------------------------------------------------------

double Pedigree::getLogLikelihoodAtPosteriorMeans( const Options & options )
    {

    double rv;

    // Should set allele freqs also to posterior means, and recalculate prob
    // genotypes at these freqs before calling getloglikelihood

    #if USE_SINGLE_POPULATION_SPECIAL_CASE
      if ( getK() == 1 )
	rv = getLogLikelihoodOnePop();
      else
    #endif
	{
	ThetaType thetabar;
	RhoType	  rhobar;

	rhobar.reserve( sumlogrho.size() );
	for ( unsigned i = 0 ; i < sumlogrho.size() ; ++i )
	    rhobar.push_back( exp(sumlogrho[i]/(double)(options.getTotalSamples() - options.getBurnIn())) );

	const double scale_factor = options.getTotalSamples() - options.getBurnIn();

	SumSoftmaxTheta /= scale_factor;

	// apply softmax transformation to obtain thetabar
	for ( ThetaIdx tIdx = 0 ; tIdx < getNTheta() ; ++tIdx )
	    SumSoftmaxTheta[tIdx].softmax_gt0( thetabar[tIdx] ); // Equivalent: SumSoftmaxTheta[tIdx].softmax( thetabar[tIdx], gt_0 );

	// Rescale sumsoftmaxtheta back
	SumSoftmaxTheta *= scale_factor;

	#if SEPARATE_HMM_FOR_EACH_CHROM
	    rv = 0.0;
	    for ( ChromIdxType chromIdx = 0 ; chromIdx < getNumChromosomes() ; chromIdx++ )
		{
		UpdateHMMInputs( chromIdx, options, thetabar, rhobar );
		#if THIS_IS_HOW_IT_WAS
		    LogLikelihood += getChromosome(j).getHMM()->getLogLikelihood( !isHaploid && (!Loci->isXChromosome(j) || SexIsFemale) );
		#else
		    rv += getHMM(chromIdx).getLogLikelihood();
		#endif
		}
	#else
	    updateHMMInputs( options, thetabar, rhobar );
	    rv = getHMM().getLogLikelihood();
	#endif

	}

    return rv;
    }



//-----------------------------------------------------------------------------
// getLogLikelihoodOnePop()
//
// Convenient for a single population as no arguments required
// Not as yet converted to pedigrees.  The "GetGenotypeProbs()" is the emission
// probabilities for individuals.
//-----------------------------------------------------------------------------

#if USE_SINGLE_POPULATION_SPECIAL_CASE
    double Pedigree::getLogLikelihoodOnePop() const
	{
	double rv = 0.0;

	// Loop over composite loci:
	for ( unsigned int j = 0; j < getGenome().GetNumberOfCompositeLoci(); j++ )
	    {
	    if ( ! GenotypeIsMissing(j) )
		{
		double Prob; // One pop so 1x1 array
		(*Loci)(j)->GetGenotypeProbs( &Prob, getPossibleHapPairs(j), false );
		rv += log( Prob );
		}
	    }

	return rv;
	}
#endif



//-----------------------------------------------------------------------------
// UpdateHMMInputs()
//-----------------------------------------------------------------------------

#if SEPARATE_HMM_FOR_EACH_CHROM

    void Pedigree::UpdateHMMInputs(
	    ChromIdxType	chromIdx,
	    const Options &	options ,
	    const ThetaType &	theta	,
	    const RhoType &	rho	) const
	{
	// Updates inputs to HMM for chromosome j also sets Diploid flag in
	// Chromosome (last arg of SetStateArrivalProbs)
	const bool diploid = !isHaploid && (j!=X_posn || SexIsFemale);
	const bool isRandomMating = options.isRandomMatingModel();

	const Chromosome & chrom = getChromosome( chromIdx );

	C->HMM->SetGenotypeProbs( GenotypeProbs[j], GenotypesMissing[j] );
	//   if(!SexIsFemale) cout << "chr " << j << " " << X_posn << " ";
	//   if(!diploid) cout << "haploid GenotypeProbs set for chromosome " << j << "\n";

	if ( ! options.isGlobalRho() )
	    {
	    // Set locus correlation, f, if individual- or gamete-specific rho
	    chrom.SetLocusCorrelation( rho, isRandomMating );
	    }

	// Set the state-arrival-probabilities in the HMM:
	// For individuals we do this:
	//	    "in the haploid case in random mating model, pass pointer to maternal admixture props."
	// For pedigrees, we will pass the entire Theta array into the HMM, regardless.
	#if THIS_IS_HOW_IT_WAS
	    if ( diploid || !isRandomMating )
		C->HMM->SetStateArrivalProbs(theta, (int)isRandomMating, diploid);
	    else
		C->HMM->SetStateArrivalProbs(theta + NumHiddenStates, true, false);
	#else
	    C->HMM->SetStateArrivalProbs(theta + NumHiddenStates, true, false);
	#endif

	logLikelihood.HMMisOK = false;//because forward probs in HMM have been changed
	}

#else

    /// Updates inputs to the HMM (which may seem odd considering the name).
    void Pedigree::updateHMMInputs( const Options & options, const ThetaType & /*theta*/, const RhoType & /*rho*/ ) const
	{

	#if NON_GLOBAL_RHO_WORKS
	    if ( ! options.isGlobalRho() )
		getHMM().getTPC().setRho( rho );
	#else
	    if ( ! options.isGlobalRho() )
		throw std::runtime_error( "non-global-rho not yet supported." );
	#endif

	#if SEPARATE_HMM_FOR_EACH_CHROM
	    getHMM().setTheta( theta );
	#else
	    getHMM().thetaChanged();
	#endif

	logLikelihood.HMMisOK = false; // Because forward probs in HMM have been changed
	}

#endif



//-----------------------------------------------------------------------------
// SetUniformAdmixtureProps() [protected]
// Set the initial values of theta.
//-----------------------------------------------------------------------------

void Pedigree::SetUniformAdmixtureProps()
    {
    bclib::pvector<double> uniform;
    uniform.resize( getK(), 1.0 / K );

    for ( size_t fIdx = getNTheta() ; fIdx-- != 0 ; )
	Theta[ fIdx ] = uniform;
    }



//=============================================================================
// Methods from AdmixedIndividual (overridden from PedBase):
//=============================================================================

//-----------------------------------------------------------------------------
/// @param chibIndicator Passed to CompositeLocus object.  If set to true,
///		CompositeLocus will use HapPairProbsMAP instead of HapPairProbs
///		when allelefreqs are not fixed.
//-----------------------------------------------------------------------------

void Pedigree::SetGenotypeProbs( /*ChromIdxType*/ int chrmIdx		,
				 int		      locusWithinChrom  ,
				 unsigned int	      absCompLocIdx	,
				 bool		      chibIndicator	)
    {

    #if 1 // Not finished porting yet
	fprintf( stderr, "Warning: unimplemented Pedigree::SetGenotypeProbs()\n" );
    #else


	#if THE_WAY_IT_WAS_BUT_CLEANED_UP_A_LITTLE

	    // === Code copied from AdmixedIndividual, may or may not be appropriate here ===

	    if( ! GenotypesMissing[chrmIdx][locusWithinChrom] )
		{
		if ( !isHaploid && ((chrmIdx != (int)X_posn) || SexIsFemale) )
		    { //diploid genotype
		    (*Loci)(locus)->GetGenotypeProbs(   GenotypeProbs[chrmIdx] + locusWithinChrom*NumHiddenStates*NumHiddenStates,
							PossibleHapPairs[locus],
							chibindicator);
		    }
		else
		    { //haploid genotype
		    (*Loci)(locus)->GetHaploidGenotypeProbs(
							GenotypeProbs[chrmIdx] + locusWithinChrom*NumHiddenStates,
							PossibleHapPairs[locus],
							chibindicator);
		    }
		}
	    else
		{
		if ( !isHaploid && (chrmIdx!=(int)X_posn || SexIsFemale))  //diploid genotype
		    for ( int k = 0; k < NumHiddenStates*NumHiddenStates; ++k )
			GenotypeProbs[chrmIdx][locusWithinChrom*NumHiddenStates*NumHiddenStates + k] = 1.0;
		else //haploid genotype
		    for ( int k = 0; k < NumHiddenStates; ++k )
			GenotypeProbs[chrmIdx][locusWithinChrom*NumHiddenStates + k] = 1.0;
		}

	#else // ==== Re-implemented for pedigrees: ====

	    if ( isHaploid(  ) )
		{ //diploid genotype
		(*Loci)(locus)->GetGenotypeProbs(	GenotypeProbs[chrmIdx] + locusWithinChrom*NumHiddenStates*NumHiddenStates,
						    PossibleHapPairs[locus],
						    chibindicator);
		}
	    else
		{ //haploid genotype
		(*Loci)(locus)->GetHaploidGenotypeProbs(GenotypeProbs[chrmIdx] + locusWithinChrom*NumHiddenStates,
							PossibleHapPairs[locus],
							chibindicator);
		}

	#endif

    #endif

    }



// draw initial values for admixture proportions theta from Dirichlet prior
void Pedigree::drawInitialAdmixtureProps(const std::vector<std::vector<double> > &alpha)
    {
      const PopIdx K = getK();

      for ( FounderIdx fIdx = 0 ; fIdx < getNFounders() ; ++fIdx )
	{

	double sum = 0.0;

	for ( PopIdx k = 0 ; k < K ; ++k )
	    sum += alpha[fIdx][k];

	for ( PopIdx k = 0 ; k < K ; ++k )
	    {
	    #if 1 // DEBUG
		const double a = alpha[fIdx][k];
		const double b = a / sum;
		pvector<double> & pv = thetahat[fIdx];
		pv.at(k) = b;
	    #else
		thetahat[fIdx][k] = alpha[fIdx][k] / sum; // set thetahat to prior mean
	    #endif
	    dirparams[k] = alpha[fIdx][k];
	    }

	// draw theta from Dirichlet with parameters dirparams
	bclib::Rand::gendirichlet< const std::vector<double> , bclib::pvector<double> >( K, dirparams.getVector_unsafe(), Theta[fIdx] );

	}

    }



#if NEEDED_ONLY_FOR_CONJUGATE_UPDATE
    //************** Updating (Public) **********************************************************
    void Pedigree::ResetSufficientStats() {
	if(NumHiddenStates>1) {
	  // ** reset SumLocusAncestry to zero
	  for(int j = 0; j < NumHiddenStates *2; ++j) {
	    SumLocusAncestry[j] = 0;
	    SumLocusAncestry_X[j] = 0;
	  }

	  //SumNumArrivals is the number of arrivals between each pair of adjacent loci
	  fill(SumNumArrivals.begin(), SumNumArrivals.end(), 0);
	}
      }
#endif



} // ---- end namespace genepi
