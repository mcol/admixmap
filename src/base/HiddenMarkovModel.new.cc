//=============================================================================
//
// Copyright (C) 2009  David D. Favro
// Portions Copyright (C) 2010 Marco Colombo
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
/// Implementation of the genepi::HiddenMarkovModel class.
//=============================================================================

#include "HiddenMarkovModel.new.h"
#include "HVIterator.h"

#include <cmath> // log()


#define CHECK_ALL_LIKELIHOODS	0 ///< Validation check that dot-products of all alpha & beta are "equal"
				  ///< Which only works if we're _not_ renormalizing on-the-fly
#define USE_LIBOIL		0
#define DEBUG_TRANSRECURSION	0 ///< Debugging output for TransRecursion



#if USE_LIBOIL
    #include <liboil/liboil.h>
#endif



namespace genepi { // ----



//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------

HiddenMarkovModel::HiddenMarkovModel( const Pedigree &	_ped	 ,
				      TransProbCache &	_tpCache ,
				      IsXChromType	_isX	 ,
				      const AdmixType * _theta	 ) :
	ped		    ( &_ped	) ,
	tpCache		    ( &_tpCache ) ,
	isX		    ( _isX	) ,
	theta		    ( _theta	) ,
	dirtyForwards	    ( true	) ,
	dirtyBackwards	    ( true	) ,
	dirtyPi		    ( true	) ,
	dirtyCondStateProbs ( true	)
    {

    // Find the starting locus index, and number of loci.  We assume that if any
    // X-chromosome loci are modeled, it is the last chromosome, i.e. those loci
    // are at the end of the locus-list.
    // FIXME-XHMM: this should be enforced, or the data massaged so that it is true.
    const SimpleLocusArray & loci = getPed().getSLoci();
    const SLocIdxType nExtLoci = loci.size();
    if ( isX == CHR_IS_X )
	{
	firstLocus = INT_MAX;
	for ( SLocIdxType extLoc = 0 ; extLoc < nExtLoci ; ++extLoc )
	    if ( loci[extLoc].isXChrom() == CHR_IS_X )
		{
		firstLocus = extLoc;
		break;
		}
	nLoci = (firstLocus == INT_MAX) ? 0 : (nExtLoci - firstLocus);
	}
    else
	{
	firstLocus = 0;
	nLoci = 0;
	for ( SLocIdxType extLoc = nExtLoci ; extLoc-- != 0 ; )
	    if ( loci[extLoc].isXChrom() != CHR_IS_X )
		{
		nLoci = extLoc + 1;
		break;
		}
	}

    #if USE_LIBOIL
	oil_init();
    #endif

    //------------------------------------------------------------------------
    // Reserve memory space
    //------------------------------------------------------------------------

    alpha.resize( getNLoci() );
    beta.resize( getNLoci() );
    condStateProbs.resize( getNLoci() );

    }



//-----------------------------------------------------------------------------
// setTheta()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::setTheta( const AdmixType * nv )
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
    dirtyPi = true;
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
// assureNotDirty()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::assureNotDirty() const
    {
    if ( dirtyForwards || dirtyBackwards )
	computeForwardsBackwards();
    }



//-----------------------------------------------------------------------------
// computeForwardsBackwards()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::computeForwardsBackwards() const
    {

    #if HMM_PARALLELIZE_FWD_BKWD && defined(_OPENMP)


    // This is required by both computeForwards() and computeBackwards().  By
    // doing it here, prior to entering the parallel sections, we avoid a race
    // condition, without needing to create a mutex-protected critical-section.
    if ( dirtyPi )
	computeStationaryDistr( getPed().getHSS(0) );


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
// transRecursion() [private]
//
/// Recursive function to compute the product of the transition probability
/// matrix times alpha.
/// TransRecursion works with probability vectors of dimension corresponding
/// to the full state space size.
//
/// @param alpha  The probability vector at the "from" locus
/// @param res    The probability vector at the "to" locus
/// @param f      The probability of 0 arrivals for founder gametes
/// @param g      The probability of 0 arrivals for meiosis
/// @param h      Admixture proportions scaled by (1-f)
/// @param z      Hidden variables index corresponding to the recursion level
/// @param sz     Size of the probability vector having fixed the first z
///               hidden variables to a certain state
//-----------------------------------------------------------------------------

void HiddenMarkovModel::transRecursion(const double *alpha, double *res,
                                       double f, double g,
                                       const AdmixType& h,
                                       HVIterator z, size_t sz ) const {

  // Check if we've reached the bottom of the recursion
  if (z.leftToIterate() == 1) {

    double asum = 0.0;
    for (size_t i = 0; i < sz; ++i)
      asum += alpha[i];

    // Ancestry
    if (z.isOnAncestry()) {

      Pedigree::GameteType whichOne;
      const Pedigree::FounderIdx idx = ped->founderOfGameteIdx( z.getCurrentAncestry(), whichOne, isX );

      for (size_t i = 0; i < sz; ++i)
        res[i] = f * alpha[i] + h[idx][i] * asum;
    }

    // Meiosis
    else {
      for (size_t i = 0; i < sz; ++i)
        res[i] = g * alpha[i] + (1 - g) * 0.5 * asum;
    }
  }

  // Recursive computations
  else {

    // Advance the iterator that will enter the recursion
    HVIterator z_p1( z );
    z_p1.advance();

    const size_t nValues = z.nValues();
    sz /= nValues;

    cvector<double> asum(sz);

    // Enter the recursion
    double *alpha_tilde = res;
    for (size_t i = 0; i < nValues; ++i) {

      // Divide the alpha array in nValues parts of size sz, each corresponding
      // to one of the possible states for the current hidden variable
      const double *alpha_i = &alpha[i * sz];
      double *alpha_tilde_i = &alpha_tilde[i * sz];
      transRecursion( alpha_i, alpha_tilde_i, f, g, h, z_p1, sz );

      // Accumulate the sum the alpha_tilde_i arrays
      for (size_t j = 0; j < sz; ++j)
        asum[j] += alpha_tilde_i[j];
    }

    // Ancestry
    if (z.isOnAncestry()) {

      Pedigree::GameteType whichOne;
      const Pedigree::FounderIdx idx = ped->founderOfGameteIdx( z.getCurrentAncestry(), whichOne, isX );

      for (size_t i = 0; i < sz * nValues; ++i)
        res[i] = f * alpha_tilde[i] + h[idx][i / sz] * asum[i % sz];
    }

    // Meiosis
    else {
      for (size_t i = 0; i < sz * nValues; ++i)
        res[i] = g * alpha_tilde[i] + (1 - g) * 0.5 * asum[i % sz];
    }
  }
}




//------------------------------------------------------------------------
//
// computeStationaryDistr() [private]
//
/// Computation of the stationary distribution (pi).
//
//------------------------------------------------------------------------

void HiddenMarkovModel::computeStationaryDistr( const HiddenStateSpace & hss_0 ) const
    {

    // Vectors for the stationary distribution and its inverse, of dimension
    // equal to the full state-space size.
    Pi	   .resize( hss_0.getNStates( isX ) );
    piInv  .resize( hss_0.getNStates( isX ) );

    // Probability of any given IV is 1/(2^M): perhaps this should be cached in
    // the pedigree or HSS?
    const double prob_of_each_iv = 1.0 / (1LL << getPed().getNMeiosis(isX));


    const AdmixType & th = *theta;


    // Here we compute the stationary distribution (Pi) and store its inverse
    // in piInv. To ensure that we compute it for all elements in the hidden
    // state space, we force the iterator to not skip the elements for which
    // the emission probability is zero by setting the third parameter of
    // the constructor to false.

    for ( HiddenStateSpace::Iterator it( hss_0, isX, false ) ; it ; ++it )
	{
	double pi = prob_of_each_iv;
	const AncestryVector & av = it.getAV();
	for ( AncestryVector::FGIdx idx = av.size(isX) ; idx-- != 0 ; )
	    pi *= th[ av.founderOf(idx,isX) ] [ av.at(idx,isX) ];
	Pi   [ it.getOverallIndex() ] = pi;
	piInv[ it.getOverallIndex() ] = 1 / pi;
	}


    dirtyPi = false;

    }



//-----------------------------------------------------------------------------
// recursionProbs()
//
/// TransRecursion works with probability vectors of dimension equal to the
/// full state space size.  Therefore, the probability vectors at the "from"
/// and "to" loci are copied to and from dense vectors which explicitly contain
/// the zeros in the positions corresponding to the elements for which the
/// emission probability is zero.
//-----------------------------------------------------------------------------

void HiddenMarkovModel::recursionProbs( double	 		 f	 ,
					double 	 		 g	 ,
					const HiddenStateSpace & fr_hss	 ,
					const HiddenStateSpace & to_hss	 ,
					const ProbsAtLocusType & frProbs ,
					ProbsAtLocusType &	 toProbs ) const
    {

    const size_t size = fr_hss.getNStates( isX );

    ProbsAtLocusType frProbsDense(size);
    ProbsAtLocusType toProbsDense(size);

    // Copy the nonzero elements of the sparse frProbs into frProbsDense
    for ( HiddenStateSpace::Iterator it( fr_hss, isX ) ; it ; ++it )
	frProbsDense[ it.getOverallIndex() ] = frProbs[ it.getNon0Index() ];

    AdmixType h( *theta );
    h.scaleBy(1 - f);

    HVIterator z( getPed(), isX );
    transRecursion(frProbsDense.data_unsafe(), toProbsDense.data_unsafe(),
		   f, g, h, z, size );

    // Ensure that the probability vector has been resized to the size of
    // the hidden-state-space
    toProbs.resize( to_hss.getNNon0() );

    // Copy the nonzero elements of toProbsDense into the sparse toProbs
    for ( HiddenStateSpace::Iterator it( to_hss, isX ) ; it ; ++it )
	toProbs[ it.getNon0Index() ] = toProbsDense[ it.getOverallIndex() ];

    }



//-----------------------------------------------------------------------------
// initStaticAlpha()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::initStaticAlpha( ProbsAtLocusType &	  alpha0	,
					 const HiddenStateSpace & hss		,
					 double &		  normalize_sum ) const
    {

    alpha0.resize( hss.getNNon0() );

    for ( HiddenStateSpace::Iterator it( hss, isX ) ; it ; ++it )
	{
	const double val = it.getEProb() * Pi[ it.getOverallIndex() ];
	alpha0[ it.getNon0Index() ] = val;
	normalize_sum += val;
	}

    }



//-----------------------------------------------------------------------------
// computeForwards()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::computeForwards() const
    {

    // Verify that all input for the model has been specified:
    gp_assert( theta != 0 );


    const HiddenStateSpace & hss_0 = getHSS( 0 );

    if ( dirtyPi )
	computeStationaryDistr( hss_0 );


    //------------------------------------------------------------------------
    // Compute alpha_0 (forward probabilities for locus #0):
    //------------------------------------------------------------------------

    ProbsAtLocusType & alpha_0 = alpha[ 0 ];

    double normalize_sum = 0.0;

    initStaticAlpha( alpha_0, hss_0, normalize_sum );



    //------------------------------------------------------------------------
    // Fill in the alpha array for the rest of the loci:
    //------------------------------------------------------------------------

    norm_log_sum_alpha = 0.0;

    for ( IntSLocIdx t = 1 ; t < getNLoci() ; ++t )
	{

	const IntSLocIdx t_m1 = t - 1;

	const HiddenStateSpace & hss_t	 = getHSS( t );
	ProbsAtLocusType &	 alpha_t = alpha[ t ];

	const HiddenStateSpace & hss_t_m1   = getHSS( t_m1 );
	ProbsAtLocusType &	 alpha_t_m1 = alpha[ t_m1 ];

	norm_log_sum_alpha += log( normalize_sum );
	gp_assert( normalize_sum != 0.0 );
	normalize_sum = 1.0 / normalize_sum;
	for ( HiddenStateSpace::Non0IdxType j = hss_t_m1.getNNon0() ; j-- != 0 ; )
	    alpha_t_m1[ j ] *= normalize_sum;
	normalize_sum = 0.0;

	#if DEBUG_TRANSRECURSION
	    cout << "\n* alpha[" << t - 1 << "] norm\n";
	    for (size_t i = 0; i < alpha_t_m1.size(); ++i)
	      cout << alpha_t_m1[i] << endl;
	#endif

	// This is the main recursion: compute alpha[t] from alpha[t-1]

	// These factors should be pre-computed (as arrays) and cached; they
	// only need be updated when theta changes (f and g) or when rho
	// changes (f).
	const double f = TransProbCache::computeF( getPed().getSLoci(), extLoc(t_m1), tpCache->getRho() );
	const double g = TransProbCache::computeG( getPed().getSLoci(), extLoc(t_m1) );

	// Do the matrix multiplication by the transition probabilities.
	recursionProbs( f, g, hss_t_m1, hss_t, alpha_t_m1, alpha_t );

	#if DEBUG_TRANSRECURSION
	    cout << "\n* alpha[" << t - 1 << "] after\n";
	    for (size_t i = 0; i < alpha_t.size(); ++i)
	      cout << alpha_t[i] << endl;
	#endif

	// Next multiply by the emission probabilities, simultaneously
	// accumulating the normalization factor for the next iteration:
	for ( HiddenStateSpace::Iterator to_it( hss_t, isX ) ; to_it ; ++to_it )
	    {
	    double & val = alpha_t[ to_it->getNon0Index() ];
	    val *= to_it->getEProb();
	    normalize_sum += val;
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



    //------------------------------------------------------------------------
    // Compute beta_T-1 (backwards probabilities for locus # T-1):
    //------------------------------------------------------------------------

    gp_assert( getNLoci() >= 1 );
    IntSLocIdx t = getNLoci() - 1;
    const HiddenStateSpace & hss_Tm1  = getHSS( t );
    ProbsAtLocusType &	     beta_Tm1 = beta[ t ];

    beta_Tm1.resize( hss_Tm1.getNNon0() );

    for ( HiddenStateSpace::Non0IdxType i = hss_Tm1.getNNon0() ; i-- != 0 ; )
	beta_Tm1[ i ] = 1.0;


    if ( dirtyPi )
	computeStationaryDistr( hss_Tm1 );


    //------------------------------------------------------------------------
    // Fill in the beta array for the rest of the loci:
    //------------------------------------------------------------------------

    norm_log_sum_beta = 0.0;

    IntSLocIdx t_p1 = t; // t plus 1
    while ( t-- != 0 )
	{

	const HiddenStateSpace & hss_t	   = getHSS( t );    // HSS at locus t
	ProbsAtLocusType &	 beta_t	   = beta[ t ];	     // beta at locus t
	const HiddenStateSpace & hss_t_p1  = getHSS( t_p1 ); // HSS at locus t+1
	ProbsAtLocusType &	 beta_t_p1 = beta[ t_p1 ];   // beta at locus t+1


	// Normalize the probabilities for beta at t+1
        double normalize_sum = 0.0;
        for ( HiddenStateSpace::Non0IdxType j = hss_t_p1.getNNon0() ; j-- != 0 ; )
            normalize_sum += beta_t_p1[ j ];

        norm_log_sum_beta += log( normalize_sum );
        gp_assert( normalize_sum != 0.0 );
        normalize_sum = 1.0 / normalize_sum;
        for ( HiddenStateSpace::Non0IdxType j = hss_t_p1.getNNon0() ; j-- != 0 ; )
            beta_t_p1[ j ] *= normalize_sum;


	// These factors should be pre-computed (as arrays) and cached; they
	// only need be updated when theta changes (f and g) or when rho
	// changes (f).
	const double f = TransProbCache::computeF( getPed().getSLoci(), extLoc(t), tpCache->getRho() );
	const double g = TransProbCache::computeG( getPed().getSLoci(), extLoc(t) );

	// Pre-multiply beta[t+1] by the emission probabilities at t+1 (making a copy)
	ProbsAtLocusType beta_t_p1_mult( beta_t_p1 );
	for ( HiddenStateSpace::Iterator fr_it( hss_t_p1, isX ) ; fr_it ; ++fr_it )
	    beta_t_p1_mult[ fr_it->getNon0Index() ] *= fr_it->getEProb() * Pi[ fr_it->getOverallIndex() ];

	// Do the matrix multiplication by the transition probabilities.
	recursionProbs( f, g, hss_t_p1, hss_t, beta_t_p1_mult, beta_t );

	for ( HiddenStateSpace::Iterator fr_it( hss_t, isX ) ; fr_it ; ++fr_it )
	    beta_t[ fr_it->getNon0Index() ] *= piInv[ fr_it->getOverallIndex() ];

	#if DEBUG_TRANSRECURSION
	    cout << "\n* beta[" << t << "] after\n";
	    for (size_t i = 0; i < beta_t.size(); ++i)
	      cout << beta_t[i] << endl;
	#endif

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

    // This happens when an HMM is instantiated for the X chromosome, but the
    // model's input data contained no loci on the X chromosome:
    if ( nLoci == 0 )
	return 0.0;

    if ( dirtyForwards )
	computeForwards();


    IntSLocIdx t = getNLoci();
    const ProbsAtLocusType & alpha_Tm1 = alpha[ --t ];

    #if USE_LIBOIL
	double rv;
	oil_sum_f64( &rv, alpha_Tm1.data_unsafe(), sizeof(double), alpha_Tm1.size() );
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
    rv += norm_log_sum_alpha;


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
	{

	gp_assert_eq( getNLoci(), alpha.size() );
	gp_assert_eq( getNLoci(), beta.size() );
	gp_assert_eq( getNLoci(), condStateProbs.size() );

	for ( IntSLocIdx sLoc = getNLoci() ; sLoc-- != 0 ; )
	    {

	    ProbsAtLocusType & locAlpha = alpha		[ sLoc ];
	    ProbsAtLocusType & locBeta	= beta		[ sLoc ];
	    ProbsAtLocusType & locProbs = condStateProbs[ sLoc ];
	    gp_assert_eq( locAlpha.size(), locBeta.size() );

	    locProbs.resize( locAlpha.size() );

	    #if USE_LIBOIL

		double normSum;
		oil_multiply_f64( locProbs.data_unsafe(), locAlpha.data_unsafe(), locBeta.data_unsafe(), locAlpha.size() );
		oil_sum_f64( &normSum, locProbs.data_unsafe(), sizeof(double), locProbs.size() );
		normSum = 1.0 / normSum;
		oil_scalarmult_f64( locProbs.data_unsafe(), sizeof(double),
				    locProbs.data_unsafe(), sizeof(double),
				    &normSum, locProbs.size() );

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
 	}

    dirtyCondStateProbs = false;

    return condStateProbs[ intLoc(t) ];

    }



} // ---- end namespace genepi
