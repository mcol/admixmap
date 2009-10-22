//=============================================================================
//
// Copyright (C) 2009  David D. Favro
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 3 as published by the Free
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
/// \file HiddenMarkovModel.new.cc
/// Implementation of the HiddenMarkovModel class.
//=============================================================================

#include "HiddenMarkovModel.new.h"

#include <cmath> // log()


#define CHECK_ALL_LIKELIHOODS	0 ///< Validation check that dot-products of all alpha & beta are "equal"
				  ///< Which only works if we're _not_ renormalizing on-the-fly
#define USE_LIBOIL		0



#if USE_LIBOIL
    #include <liboil/liboil.h>
#endif



namespace genepi { // ----



//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------

HiddenMarkovModel::HiddenMarkovModel( const Pedigree &	_ped	 ,
				      TransProbCache &	_tpCache ,
				      const ThetaType *	_theta	 ) :
	ped	      ( &_ped	  ) ,
	tpCache	      ( &_tpCache ) ,
	theta	      ( _theta	  ) ,
	dirtyForwards ( true	  ) ,
	dirtyBackwards( true	  ) ,
	dirtyCondStateProbs( true )
    {
    #if USE_LIBOIL
	oil_init();
    #endif

    condStateProbs.resize( getNLoci() );
    }



//-----------------------------------------------------------------------------
// setTheta()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::setTheta( const ThetaType * nv )
    {
    theta = nv;
    tpCache->setMu( *nv );
    thetaChanged();
    }



//-----------------------------------------------------------------------------
// thetaChanged()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::thetaChanged() const
    {
    transProbsChanged();
    }



//-----------------------------------------------------------------------------
// transProbsChanged()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::transProbsChanged() const
    {
    dirtyForwards = true;
    dirtyBackwards = true;

    // This is not really necessary, forwards-backwards recomputation invalidates this cache:
    dirtyCondStateProbs = true;
    }



//-----------------------------------------------------------------------------
// computeForwardsBackwards()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::computeForwardsBackwards() const
    {

    #if HMM_PARALLELIZE_FWD_BKWD && defined(_OPENMP)
      #pragma omp parallel default(shared)
      { // begin parallel

      #pragma omp sections
      { // begin sections

      #pragma omp section
      { // begin compute-forwards section
    #endif


    if ( dirtyForwards )
	computeForwards();


    #if HMM_PARALLELIZE_FWD_BKWD && defined(_OPENMP)
      } // end compute-forwards section

      #pragma omp section
      { // begin compute-backwards section
    #endif


    if ( dirtyBackwards )
	computeBackwards();


    #if HMM_PARALLELIZE_FWD_BKWD && defined(_OPENMP)
      } // end compute-backwards section
      } // end sections
      } // end parallel
    #endif


    }



//-----------------------------------------------------------------------------
// computeForwards()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::computeForwards() const
    {

    // Verify that all input for the model has been specified:
    gp_assert( theta != 0 );


    const SLocIdxType T = getNLoci();



    //------------------------------------------------------------------------
    // Reserve memory space
    //------------------------------------------------------------------------

    alpha.resize( T );



    //------------------------------------------------------------------------
    // Compute alpha_0 (forward probabilities for locus #0):
    //------------------------------------------------------------------------

    const HiddenStateSpace & hss_0   = getPed().getStateProbs( 0 );
    ProbsAtLocusType &	     alpha_0 = alpha[ 0 ];

    #if HMM_OTF_RENORM
	double normalize_sum = 0.0;
    #endif

    alpha_0.resize( hss_0.getNNon0() );

    // Probability of any given IV is 1/(2^M): perhaps this should be cached in
    // the pedigree or HSS?
    const double prob_of_each_iv = 1.0 / (1LL << getPed().getNMeiosis());

    const ThetaType th = *theta;
    for ( HiddenStateSpace::Iterator it( hss_0 ) ; it ; ++it )
	{
	double pi = prob_of_each_iv;
	const AncestryVector & av = it.getAV();
	for ( AncestryVector::IdxType idx = av.size() ; idx-- != 0 ; )
	    pi *= th[ AncestryVector::founderOf(idx) ] [ av.at(idx) ];
	const double val = pi * it.getEProb();
	alpha_0[ it.getNon0Index() ] = val;
	#if HMM_OTF_RENORM
	    normalize_sum += val;
	#endif
	}



    //------------------------------------------------------------------------
    // Fill in the alpha array for the rest of the loci:
    //------------------------------------------------------------------------

    #if HMM_OTF_RENORM
	norm_log_sum_alpha = 0.0;
    #endif

    for ( SLocIdxType t = 1 ; t < T ; ++t )
	{
	const SLocIdxType t_m1 = t - 1;

	const HiddenStateSpace & hss_t	    = getPed().getStateProbs( t );
	ProbsAtLocusType &	 alpha_t    = alpha[ t ];
	const HiddenStateSpace & hss_t_m1   = getPed().getStateProbs( t_m1 );
	ProbsAtLocusType &	 alpha_t_m1 = alpha[ t_m1 ];

	#if HMM_OTF_RENORM
	    norm_log_sum_alpha += log( normalize_sum );
	    gp_assert( normalize_sum != 0.0 );
	    normalize_sum = 1.0 / normalize_sum;
	    for ( HiddenStateSpace::Non0IdxType j = hss_t_m1.getNNon0() ; j-- != 0 ; )
		alpha_t_m1[ j ] *= normalize_sum;
	    normalize_sum = 0.0;
	#endif

	alpha_t.resize( hss_t.getNNon0() );

	for ( HiddenStateSpace::Iterator to_it( hss_t ) ; to_it ; ++to_it )
	    {
	    double val = 0.0;
	    for ( HiddenStateSpace::Iterator fr_it( hss_t_m1 ) ; fr_it ; ++fr_it )
		val += alpha_t_m1[ fr_it->getNon0Index() ] * tpCache->getProb( t_m1, fr_it, to_it );
	    val *= to_it->getEProb();
	    #if HMM_OTF_RENORM
		normalize_sum += val;
	    #endif
	    alpha_t[ to_it->getNon0Index() ] = val;
	    }
	}

    dirtyForwards = false;
    dirtyCondStateProbs = true;

    }



//-----------------------------------------------------------------------------
// computeBackwards()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::computeBackwards() const
    {

    // Verify that all input for the model has been specified:
    gp_assert( theta != 0 );


    const SLocIdxType T = getNLoci();



    //------------------------------------------------------------------------
    // Reserve memory space
    //------------------------------------------------------------------------

    beta.resize( T );


    //------------------------------------------------------------------------
    // Compute beta_T-1 (forward probabilities for locus # T-1):
    //------------------------------------------------------------------------

    gp_assert( T >= 1 );
    SLocIdxType t = T - 1;
    const HiddenStateSpace & hss_Tm1  = getPed().getStateProbs( t );
    ProbsAtLocusType &	     beta_Tm1 = beta[ t ];

    beta_Tm1.resize( hss_Tm1.getNNon0() );

    for ( HiddenStateSpace::Non0IdxType i = hss_Tm1.getNNon0() ; i-- != 0 ; )
	beta_Tm1[ i ] = 1.0;

    #if HMM_OTF_RENORM
	double normalize_sum = hss_Tm1.getNNon0(); // 1.0 + 1.0 + ... + 1.0
    #endif



    //------------------------------------------------------------------------
    // Fill in the beta array for the rest of the loci:
    //------------------------------------------------------------------------

    #if HMM_OTF_RENORM
	norm_log_sum_beta = 0.0;
    #endif

    SLocIdxType t_p1 = t; // t plus 1
    while ( t-- != 0 )
	{
	const HiddenStateSpace & hss_t	   = getPed().getStateProbs( t );    // HSS at locus t
	ProbsAtLocusType &	 beta_t	   = beta[ t ];			     // beta at locus t
	const HiddenStateSpace & hss_t_p1  = getPed().getStateProbs( t_p1 ); // HSS at locus t+1
	ProbsAtLocusType &	 beta_t_p1 = beta[ t_p1 ];		     // beta at locus t+1

	// Normalize the probabilities for beta at t+1
	#if HMM_OTF_RENORM
	    norm_log_sum_beta += log( normalize_sum );
	    gp_assert( normalize_sum != 0.0 );
	    normalize_sum = 1.0 / normalize_sum;
	    for ( HiddenStateSpace::Non0IdxType j = hss_t_p1.getNNon0() ; j-- != 0 ; )
		beta_t_p1[ j ] *= normalize_sum;
	    normalize_sum = 0.0;
	#endif

	beta_t.resize( hss_t.getNNon0() );

	for ( HiddenStateSpace::Iterator to_it( hss_t ) ; to_it ; ++to_it )
	    {
	    double val = 0.0;

	    for ( HiddenStateSpace::Iterator fr_it( hss_t_p1 ) ; fr_it ; ++fr_it )
		val += beta_t_p1[ fr_it->getNon0Index() ] * tpCache->getProb( t, to_it, fr_it ) * fr_it->getEProb();

	    #if HMM_OTF_RENORM
	      normalize_sum += val;
	    #endif

	    beta_t[ to_it->getNon0Index() ] = val;
	    }

	t_p1 = t;
	}

    dirtyBackwards = false;
    dirtyCondStateProbs = true;

    }



//-----------------------------------------------------------------------------
// getLogLikelihood()
//-----------------------------------------------------------------------------

double HiddenMarkovModel::getLogLikelihood() const
    {
    if ( dirtyForwards )
	computeForwards();


    SLocIdxType t = getNLoci();
    const ProbsAtLocusType & alpha_Tm1 = alpha[ --t ];

    #if USE_LIBOIL
	double rv;
	oil_sum_f64( &rv, alpha_Tm1.data(), 1, alpha_Tm1.size() );
    #else
	double rv = 0.0;

	// Parallelizing this probably gains little after overhead, but works:
	#if 0 && defined(_OPENMP)
	    #warning expect 1 "iteration variable is unsigned" warning here
	    #pragma omp parallel for default(shared) reduction(+:rv)
	    for ( HiddenStateSpace::Non0IdxType stIdx = 0 ; stIdx < alpha_Tm1.size() ; ++stIdx )
	#else
	    for ( HiddenStateSpace::Non0IdxType stIdx = alpha_Tm1.size() ; stIdx-- != 0 ; )
	#endif
		rv += alpha_Tm1[ stIdx ];
    #endif

    rv = log(rv);
    #if HMM_OTF_RENORM
	rv += norm_log_sum_alpha;
    #endif


    //-----------------------------------------------------------------------------
    // This is strictly for debugging: consistency check.
    //-----------------------------------------------------------------------------
    #if CHECK_ALL_LIKELIHOODS
	while ( t-- != 0 ) // This could be parallelized.
	    {
	    double l_at_t = 0.0;
	    const ProbsAtLocusType & alpha_t = alpha[ t ];
	    const ProbsAtLocusType & beta_t  = beta [ t ];
	    gp_assert_eq( alpha_t.size(), beta_t.size() );
	    for ( HiddenStateSpace::Non0IdxType stIdx = alpha_t.size() ; stIdx-- != 0 ; )
		l_at_t += alpha_t[ stIdx ] * beta_t[ stIdx ];

	    #if 0 // Well, not _exactly_ equal
		gp_assert_eq( rv, l_at_t );
	    #else
		double error = ((rv + 1.0) / (l_at_t + 1.0)) - 1.0;
		if ( error < 0.0 )
		    error = - error;
		gp_assert_lt( error, 0.000001 );
	    #endif
	    }
    #endif


    // FIXME: under what circumstances does this happen?  Copied out of AdmixedIndividual.
    #if AGGRESSIZE_RANGE_CHECK
	if ( isnan(rv) )
	    throw std::runtime_error( "HMM returns log-likelihood as nan (not a number)" );
    #endif


    return rv;
    }



//-----------------------------------------------------------------------------
// getCondStateProbsAtLocus()
//-----------------------------------------------------------------------------

const cvector<double> & HiddenMarkovModel::getCondStateProbsAtLocus( SLocIdxType t ) const
    {

    if ( dirtyForwards || dirtyBackwards )
	computeForwardsBackwards();

    if ( dirtyCondStateProbs )
	for ( SLocIdxType sLoc = getNLoci() ; sLoc-- != 0 ; )
	    {

	    gp_assert_eq( getNLoci(), alpha.size() );
	    gp_assert_eq( getNLoci(), beta.size() );
	    gp_assert_eq( getNLoci(), condStateProbs.size() );

	    ProbsAtLocusType & locAlpha = alpha		[ sLoc ];
	    ProbsAtLocusType & locBeta	= beta		[ sLoc ];
	    ProbsAtLocusType & locProbs = condStateProbs[ sLoc ];

	    locProbs.resize( locAlpha.size() );

	    #if USE_LIBOIL

		double normSum;
		oil_multiply_f64( locProbs.data(), locAlpha.data(), locBeta.data(), locAlpha.size() );
		oil_sum_f64( locProbs.data(), &normSum, 1, locProbs.size() );
		normSum = 1.0 / normSum;
		oil_scalarmult_f64( locProbs.data(), 1, locProbs.data(), 1, &normSum, locProbs.size() );

	    #else

		double normSum = 0.0;
		for ( size_t stIdx = locAlpha.size() ; stIdx-- != 0 ; )
		    {
		    const double el = locAlpha[stIdx] * locBeta[stIdx];
		    locProbs[stIdx] = el;
		    normSum += el;
		    }
		normSum = 1.0 / normSum;
		for ( size_t stIdx = locProbs.size() ; stIdx-- != 0 ; )
		    locProbs.at_unsafe(stIdx) *= normSum;

	    #endif

	    }

    return condStateProbs[ t ];

    }



} // ---- end namespace genepi
