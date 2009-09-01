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
#define DEBUG_PRINT_ALPHA_SUMS	0
#define USE_LIBOIL		0


#if USE_LIBOIL
    #include <liboil/liboil.h>
#endif



namespace genepi { // ----



//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------

HiddenMarkovModel::HiddenMarkovModel( const Pedigree &	      _ped	,
				      const TransProbCache &  _tpCache  ,
				      const ThetaType *	      _theta	) :
	ped	      ( &_ped	  ) ,
	tpCache	      ( &_tpCache ) ,
	theta	      ( _theta	  ) ,
	dirtyForwards ( true	  ) ,
	dirtyBackwards( true	  )
    {
    }



//-----------------------------------------------------------------------------
// setTheta()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::setTheta( const ThetaType * nv )
    {
    theta = nv;
    thetaChanged();
    }



//-----------------------------------------------------------------------------
// thetaChanged()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::thetaChanged() const
    {
    dirtyForwards = true;
    dirtyBackwards = true;
    }



//-----------------------------------------------------------------------------
// computeForwardsBackwards()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::computeForwardsBackwards() const
    {

    // Verify that all input for the model has been specified:
    gp_assert( theta != 0 );


    const SLocIdxType T = getNLoci();


    #if HMM_PARALLELIZE_FWD_BKWD && defined(_OPENMP)
      //#pragma omp parallel default(shared) num_threads(2)
      #pragma omp parallel default(shared)
      { // begin parallel

      #pragma omp sections
      { // begin sections

      #pragma omp section
      { // begin compute-forwards section
    #endif



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

    alpha_0.resize( hss_0.getNStates() );

    // Probability of any given IV is 1/(2^M): perhaps this should be cached in
    // the pedigree or HSS?
    const double prob_of_each_iv = 1.0 / (1LL << getPed().getNMeiosis());

    // Iterator is more efficient than index here:
    for ( HiddenStateSpace::Iterator it( hss_0 ) ; it ; ++it )
	{
	double pi = prob_of_each_iv;
	const AncestryVector & av = it.getAV();
	for ( AncestryVector::IdxType idx = av.size() ; idx-- != 0 ; )
	    pi *= (*theta) [ AncestryVector::founderOf(idx) ] [ av.at(idx) ];
	const double val = pi * it.getEProb();
	alpha_0[ it.getIndex() ] = val;
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

	const HiddenStateSpace &	     hss_t	= getPed().getStateProbs( t );
	ProbsAtLocusType &		     alpha_t	= alpha[ t ];
	const HiddenStateSpace &	     hss_t_m1	= getPed().getStateProbs( t_m1 );
	ProbsAtLocusType &		     alpha_t_m1 = alpha[ t_m1 ];
	const HiddenStateSpace::StateIdxType n_st_t_m1	= hss_t_m1.getNStates();

	#if HMM_OTF_RENORM
	    norm_log_sum_alpha += log( normalize_sum );
	    gp_assert( normalize_sum != 0.0 );
	    normalize_sum = 1.0 / normalize_sum;
	    for ( HiddenStateSpace::StateIdxType j = hss_t_m1.getNStates() ; j-- != 0 ; )
		alpha_t_m1[ j ] *= normalize_sum;
	    normalize_sum = 0.0;
	#endif

	alpha_t.resize( hss_t.getNStates() );

	for ( HiddenStateSpace::StateIdxType j = hss_t.getNStates() ; j-- != 0 ; )
	    {
	    double val = 0.0;
	    for ( HiddenStateSpace::StateIdxType i = n_st_t_m1 ; i-- != 0 ; )
		val += alpha_t_m1[ i ] * tpCache->getProb( t_m1, i, j );
	    val *= hss_t.getEProb( j );
	    #if HMM_OTF_RENORM
		normalize_sum += val;
	    #endif
	    alpha_t[ j ] = val;
	    }
	}

    dirtyForwards = false;


    #if HMM_PARALLELIZE_FWD_BKWD && defined(_OPENMP)
      } // end compute-forwards section

      #pragma omp section
      { // begin compute-backwards section
    #endif



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

    beta_Tm1.resize( hss_Tm1.getNStates() );

    for ( HiddenStateSpace::StateIdxType i = hss_Tm1.getNStates() ; i-- != 0 ; )
	beta_Tm1[ i ] = 1.0;

    #if HMM_OTF_RENORM
	#if HMM_PARALLELIZE_FWD_BKWD && defined(_OPENMP)
	    double
	#endif
	normalize_sum = hss_Tm1.getNStates(); // 1.0 + 1.0 + ... + 1.0
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
	const HiddenStateSpace &	     hss_t	= getPed().getStateProbs( t );	  // HSS at locus t
	ProbsAtLocusType &		     beta_t	= beta[ t ];			  // beta at locus t
	const HiddenStateSpace &	     hss_t_p1	= getPed().getStateProbs( t_p1 ); // HSS at locus t+1
	ProbsAtLocusType &		     beta_t_p1	= beta[ t_p1 ];			  // beta at locus t+1
	const HiddenStateSpace::StateIdxType n_st_t_p1	= hss_t_p1.getNStates();	  // # of states at t+1

	// Normalize the probabilities for beta at t+1
	#if HMM_OTF_RENORM
	    norm_log_sum_beta += log( normalize_sum );
	    gp_assert( normalize_sum != 0.0 );
	    normalize_sum = 1.0 / normalize_sum;
	    for ( HiddenStateSpace::StateIdxType j = hss_t_p1.getNStates() ; j-- != 0 ; )
		beta_t_p1[ j ] *= normalize_sum;
	    normalize_sum = 0.0;
	#endif

	beta_t.resize( hss_t.getNStates() );

	for ( HiddenStateSpace::StateIdxType i = hss_t.getNStates() ; i-- != 0 ; )
	    {
	    double val = 0.0;

	    for ( HiddenStateSpace::StateIdxType j = n_st_t_p1 ; j-- != 0 ; )
		val += beta_t_p1[ j ] * tpCache->getProb( t, i, j ) * hss_t_p1.getEProb( j );

	    #if HMM_OTF_RENORM
	      normalize_sum += val;
	    #endif

	    beta_t[ i ] = val;
	    }

	t_p1 = t;
	}

    dirtyBackwards = false;


    #if HMM_PARALLELIZE_FWD_BKWD && defined(_OPENMP)
      } // end compute-backwards section
      } // end sections
      } // end parallel
    #endif



    }



//-----------------------------------------------------------------------------
// getLogLikelihood()
//-----------------------------------------------------------------------------

double HiddenMarkovModel::getLogLikelihood() const
    {
    if ( dirtyForwards )
	computeForwardsBackwards();


    #if DEBUG_PRINT_ALPHA_SUMS
	{
	printf( "\n\n ==== DEBUG ====\n" );
	for ( size_t t = 0 ; t < getNLoci() ; ++t )
	    {
	    double sum = 0.0;
	    const ProbsAtLocusType & a = alpha[ t ];
	    for ( HiddenStateSpace::StateIdxType stIdx = a.size() ; stIdx-- != 0 ; )
		sum += a[ stIdx ];
	    printf( "  Alpha-sum at locus %lu = %.16lf\n", t, sum );
	    }
	printf( "\n\n" );
	}
    #endif


    SLocIdxType t = getNLoci();
    const ProbsAtLocusType & alpha_Tm1 = alpha[ --t ];

    double rv = 0.0;

    // Parallelizing this probably gains little after overhead.
    #if defined(_OPENMP)
	#warning expect 1 "iteration variable is unsigned" warning here
	#pragma omp parallel for default(shared) reduction(+:rv)
	for ( HiddenStateSpace::StateIdxType stIdx = 0 ; stIdx < alpha_Tm1.size() ; ++stIdx )
    #else
	for ( HiddenStateSpace::StateIdxType stIdx = alpha_Tm1.size() ; stIdx-- != 0 ; )
    #endif
	    rv += alpha_Tm1[ stIdx ];

    rv = log(rv);
    #if HMM_OTF_RENORM
	rv += norm_log_sum_alpha;
    #endif


    //-----------------------------------------------------------------------------
    // This is strictly for debugging: consistency check.
    //-----------------------------------------------------------------------------
    #if CHECK_ALL_LIKELIHOODS
	// This could be parallelized.
	#if 0 && defined(_OPENMP)
	    #pragma omp parallel for default(shared)
	#endif
	while ( t-- != 0 )
	    {
	    double l_at_t = 0.0;
	    const ProbsAtLocusType & alpha_t = alpha[ t ];
	    const ProbsAtLocusType & beta_t  = beta [ t ];
	    gp_assert_eq( alpha_t.size(), beta_t.size() );
	    for ( HiddenStateSpace::StateIdxType stIdx = alpha_t.size() ; stIdx-- != 0 ; )
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


    return rv;
    }



//------------------------------------------------------------------------
// HMMBase-OVERRIDE METHODS:
//------------------------------------------------------------------------

#if HAVE_HMM_BASE

void HiddenMarkovModel::SetGenotypeProbs( const double * /*lambdain*/, const bool * /*missing*/ )
    {
    throw std::runtime_error( "SetGenotypeProbs() not implemented" );
    }

void HiddenMarkovModel::SetStateArrivalProbs( const double * /*Theta*/, int /*Mcol*/, bool /*isdiploid*/ )
    {
    throw std::runtime_error( "SetStateArrivalProbs() not implemented" );
    }

void HiddenMarkovModel::SampleHiddenStates( int * /*SStates*/, bool /*isdiploid*/ )
    {
    throw std::runtime_error( "SampleHiddenStates() not implemented" );
    }

const bclib::pvector<double> & HiddenMarkovModel::GetHiddenStateProbs( const bool /*isDiploid*/, int /*t*/ )
    {
    throw std::runtime_error( "GetHiddenStateProbs() not implemented" );
    }

double HiddenMarkovModel::getLogLikelihood( bool /*isdiploid*/ )
    {
    return getLikelihood();
    }

void HiddenMarkovModel::SampleJumpIndicators( const int * /*LocusAncestry*/, const unsigned int /*gametes*/,
		int * /*SumLocusAncestry*/, std::vector<unsigned> & /*SumNumArrivals*/,
		bool /*SampleArrivals*/, unsigned /*startlocus*/ ) const
    {
    throw std::runtime_error( "SampleJumpIndicators() not implemented" );
    }

#endif



} // ---- end namespace genepi
