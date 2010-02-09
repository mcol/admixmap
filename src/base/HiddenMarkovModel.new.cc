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
#include "HVIterator.h"

#include <cmath> // log()


#define CHECK_ALL_LIKELIHOODS	0 ///< Validation check that dot-products of all alpha & beta are "equal"
				  ///< Which only works if we're _not_ renormalizing on-the-fly
#define USE_LIBOIL		0
#define DEBUG_TRANSRECURSION	0 ///< Debugging output for TransRecursion
#define ALPHA_REMAPPING 	1 ///< Reorder the vectors for TransRecursion



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
                                       const ThetaType& h,
                                       HVIterator z, size_t sz) const {

  // Check if we've reached the bottom of the recursion
  if (z.leftToIterate() == 1) {

    double asum = 0.0;
    for (size_t i = 0; i < sz; ++i)
      asum += alpha[i];

    // Ancestry
    if (z.isOnAncestry()) {

      Pedigree::GameteType whichOne;
      const Pedigree::FounderIdx idx = ped->founderOfGameteIdx(z.getCurrentAncestry(), whichOne);

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

    // Enter the recursion
    double *alpha_tilde = res;
    for (size_t i = 0; i < nValues; ++i) {

      // Divide the alpha array in nValues parts of size sz, each corresponding
      // to one of the possible states for the current hidden variable
      const double *alpha_i = &alpha[i * sz];
      double *alpha_tilde_i = &alpha_tilde[i * sz];
      transRecursion(alpha_i, alpha_tilde_i, f, g, h, z_p1, sz);
    }

    // Sum the alpha_tilde_i arrays along the rows
    double *asum = new double[sz];
    for (size_t i = 0; i < sz; ++i) {
      asum[i] = 0.0;
      for (size_t j = 0; j < nValues; ++j)
        asum[i] += alpha_tilde[j * sz + i];
    }

    // Ancestry
    if (z.isOnAncestry()) {
      Pedigree::GameteType whichOne;
      const Pedigree::FounderIdx idx = ped->founderOfGameteIdx(z.getCurrentAncestry(), whichOne);

      for (size_t i = 0; i < sz * nValues; ++i)
        res[i] = f * alpha_tilde[i] + h[idx][i / sz] * asum[i % sz];
    }

    // Meiosis
    else {
      for (size_t i = 0; i < sz * nValues; ++i)
        res[i] = g * alpha_tilde[i] + (1 - g) * 0.5 * asum[i % sz];
    }

    delete[] asum;
  }
}


//-----------------------------------------------------------------------------
// recursionProbs()
//-----------------------------------------------------------------------------

void HiddenMarkovModel::recursionProbs( double	 		 f	   ,
					double 	 		 g	   ,
					const ThetaType &	 h	   ,
					const HiddenStateSpace & fr_hss	   ,
					const HiddenStateSpace & to_hss	   ,
					const ProbsAtLocusType & frProbs   ,
					ProbsAtLocusType &	 toProbs) const
    {
    gp_assert(frProbs.size() == toProbs.size());
  
    // create a flat array for alpha
    double *alphaflat = new double[frProbs.size()];
    double *res = new double[toProbs.size()];

#if ALPHA_REMAPPING
    memset(alphaflat, 0, frProbs.size() * sizeof(double));
    for ( HiddenStateSpace::Iterator it( fr_hss ) ; it ; ++it ) {
      const AncestryVector & av = it.getAV();
      const InheritanceVector & iv = it.getIV();
      alphaflat[fr_hss.idxMapping(av, iv)] = frProbs[it.getOverallIndex()];
    }
#else
    copy(frProbs.begin(), frProbs.end(), alphaflat);
#endif

    HVIterator z(getPed());
    transRecursion(alphaflat, res, f, g, h, z, frProbs.size());

#if ALPHA_REMAPPING
    // store res into alpha[t+1]
    for ( HiddenStateSpace::Iterator it( to_hss ) ; it ; ++it ) {
      const AncestryVector & av = it.getAV();
      const InheritanceVector & iv = it.getIV();
      toProbs[it.getOverallIndex()] = res[to_hss.idxMapping(av, iv)];
    }
#else
    copy(res, res + toProbs.size(), toProbs.begin());
#endif

    delete[] alphaflat;
    delete[] res;
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

    alpha_0.resize( hss_0.getNStates() );

    // Probability of any given IV is 1/(2^M): perhaps this should be cached in
    // the pedigree or HSS?
    const double prob_of_each_iv = 1.0 / (1LL << getPed().getNMeiosis());

    const ThetaType & th = *theta;
    for ( HiddenStateSpace::Iterator it( hss_0 ) ; it ; ++it )
	{
	double pi = prob_of_each_iv;
	const AncestryVector & av = it.getAV();
	for ( AncestryVector::FGIdx idx = av.size() ; idx-- != 0 ; )
	    pi *= th[ AncestryVector::founderOf(idx) ] [ av.at(idx) ];
	const double val = pi * it.getEProb();
	alpha_0[ it.getOverallIndex() ] = val;
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
	    for ( HiddenStateSpace::Non0IdxType j = hss_t_m1.getNStates() ; j-- != 0 ; )
		alpha_t_m1[ j ] *= normalize_sum;
	    normalize_sum = 0.0;
	#endif

#if DEBUG_TRANSRECURSION
        cout << "\n* alpha[" << t - 1 << "] norm\n";
        for (size_t i = 0; i < alpha_t_m1.size(); ++i)
          cout << alpha_t_m1[i] << endl;
#endif

	// This is the main recursion: compute alpha[t] from alpha[t-1]

	// These factors should be pre-computed (as arrays) and cached; they
	// only need be updated when theta changes (f,g,h) or when rho changes
	// (f and h).
	const double f = TransProbCache::computeF( getPed().getSLoci(), t_m1, tpCache->getRho() );
	const double g = TransProbCache::computeG( getPed().getSLoci(), t_m1 );
	ThetaType h( th );
	h *= (1 - f);

	// Do the matrix multiplication by the transition probabilities.  We
	// assure that alpha[t] has been resized to the size of the
	// hidden-state-space at locus t so that recursionProbs() knows how
	// large the space is.
	alpha_t.resize( hss_t.getNStates() );
	recursionProbs( f, g, h, hss_t_m1, hss_t, alpha_t_m1, alpha_t );

#if DEBUG_TRANSRECURSION
        cout << "\n* alpha[" << t - 1 << "] after\n";
        for (size_t i = 0; i < alpha_t.size(); ++i)
          cout << alpha_t[i] << endl;
#endif

	// Next multiply by the emission probabilities, simultaneously
	// accumulating the normalization factor for the next iteration:
	for ( HiddenStateSpace::Iterator to_it( hss_t ) ; to_it ; ++to_it )
	    {
	    double & val = alpha_t[ to_it->getOverallIndex() ];
	    val *= to_it->getEProb();
	    #if HMM_OTF_RENORM
		normalize_sum += val;
	    #endif
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
    // Compute beta_T-1 (backwards probabilities for locus # T-1):
    //------------------------------------------------------------------------

    gp_assert( T >= 1 );
    SLocIdxType t = T - 1;
    const HiddenStateSpace & hss_Tm1  = getPed().getStateProbs( t );
    ProbsAtLocusType &	     beta_Tm1 = beta[ t ];

    beta_Tm1.resize( hss_Tm1.getNStates() );

    for ( HiddenStateSpace::Non0IdxType i = hss_Tm1.getNStates() ; i-- != 0 ; )
	beta_Tm1[ i ] = 1.0;


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

	    double normalize_sum = 0.0;
	    for ( HiddenStateSpace::Non0IdxType j = hss_t_p1.getNStates() ; j-- != 0 ; )
		normalize_sum += beta_t_p1[ j ];

	    norm_log_sum_beta += log( normalize_sum );
	    gp_assert( normalize_sum != 0.0 );
	    normalize_sum = 1.0 / normalize_sum;
	    for ( HiddenStateSpace::Non0IdxType j = hss_t_p1.getNStates() ; j-- != 0 ; )
		beta_t_p1[ j ] *= normalize_sum;

	#endif


	// These factors should be pre-computed (as arrays) and cached; they
	// only need be updated when theta changes (f,g,h) or when rho changes
	// (f and h).
	const double f = TransProbCache::computeF( getPed().getSLoci(), t, tpCache->getRho() );
	const double g = TransProbCache::computeG( getPed().getSLoci(), t );
	ThetaType h( *theta );
	h *= (1 - f);


	// Pre-multiply beta[t+1] by the emission probabilities at t+1 (making a copy)
	ProbsAtLocusType beta_t_p1_mult( beta_t_p1 );
	for ( HiddenStateSpace::Iterator fr_it( hss_t_p1 ) ; fr_it ; ++fr_it )
	    beta_t_p1_mult[ fr_it->getOverallIndex() ] *= fr_it->getEProb();


	// Do the matrix multiplication by the transition probabilities.  We
	// assure that beta[t] has been resized to the size of the
	// hidden-state-space at locus t so that recursionProbs() knows how
	// large the space is.
	beta_t.resize( hss_t.getNStates() );
	recursionProbs( f, g, h, hss_t_p1, hss_t, beta_t_p1_mult, beta_t );

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
    if ( dirtyForwards )
	computeForwards();


    SLocIdxType t = getNLoci();
    const ProbsAtLocusType & alpha_Tm1 = alpha[ --t ];

    #if USE_LIBOIL && 0 // disabled, as it doesn't produce the same results
	double rv;
	oil_sum_f64( &rv, alpha_Tm1.data_unsafe(), 1, alpha_Tm1.size() );
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
		oil_multiply_f64( locProbs.data_unsafe(), locAlpha.data_unsafe(), locBeta.data_unsafe(), locAlpha.size() );
		oil_sum_f64( locProbs.data_unsafe(), &normSum, 1, locProbs.size() );
		normSum = 1.0 / normSum;
		oil_scalarmult_f64( locProbs.data_unsafe(), 1, locProbs.data_unsafe(), 1, &normSum, locProbs.size() );

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
