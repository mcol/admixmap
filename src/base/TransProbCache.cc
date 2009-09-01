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
/// \file TransProbCache.cc
/// Implementation of the TransProbCache class.
//=============================================================================

#pragma implementation
#include "TransProbCache.h"


#define DEBUG_PRINT_TP_DETAIL	0


#include <cmath>


namespace genepi { // ----



//-----------------------------------------------------------------------------
// computeFactors() [protected]
// Recompute the cached "f" values (called initially and when rho changes).
//-----------------------------------------------------------------------------

void TransProbCache::computeFactors()
    {
    // No need to recompute g since it doesn not depend on rho
    for ( int t = factors.size() ; t-- != 0 ; )
	factors[ t ].f = computeF( loci, t, rho );
    }



//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------

TransProbCache::TransProbCache( const Pedigree & _pedigree, double _rho, const MuType & _mu ) :
	pedigree ( _pedigree		) ,
	loci	 ( _pedigree.getSLoci() ) ,
	rho	 ( _rho			) ,
	mu	 ( &_mu			)
    {

    #if TPC_CACHE_MODEL == TPC_NO_CACHE

	muChanged();

	gp_assert( loci.size() != 0 ); // Unnecessary, but a bug did manifest itself this way once

	factors.resize( loci.size() - 1 );
	for ( int t = factors.size() ; t-- != 0 ; )
	    factors[ t ].g = computeG( loci, t	);
	computeFactors();

    #elif TPC_CACHE_MODEL == TPC_BIG_CACHE

	const size_t T_minus_1 = loci.size() - 1;

	// Cache f & g by locus:
	double fa[ loci.size() ];
	double ga[ loci.size() ];
	for ( SLocIdxType t = T_minus_1 ; t-- != 0 ; )
	    {
	    fa[ t ] = computeF( loci, t, rho );
	    ga[ t ] = computeG( loci, t );
	    #if DEBUG_PRINT_TP_DETAIL
		printf( "At locus %lu, f=%lf g=%lf\n", t, fa[t], ga[t] );
	    #endif
	    }


	// Assume that all loci have the same sized hidden-state-space:
	const size_t n_states = pedigree.getStateProbs(0).getNStates();
	pCache.stMultiplier  = n_states;
	pCache.locMultiplier = n_states * n_states;
	size_t n_probs = pCache.locMultiplier * T_minus_1;
	#if AGGRESSIVE_RANGE_CHECK
	    // Keep the number of probabilities (array-size) for later range checks:
	    pCache.n_probs = n_probs;
	#endif

	#if 0 // MEM-DEBUG
	    fprintf( stderr, "T-1: %lu; N-states: %lu; N-probs: %lu; mem request: %lu\n",
		    T_minus_1, n_states, n_probs, n_probs * sizeof(CacheProbType) );
	#endif
	float * const probs = new CacheProbType[ n_probs ];
	pCache.probs = probs;

	for ( SLocIdxType t = T_minus_1 ; t-- != 0 ; )
	    {
	    const double f = fa[t];
	    const double g = ga[t];

	    const HiddenStateSpace & frSpace = pedigree.getStateProbs( t     );
	    const HiddenStateSpace & toSpace = pedigree.getStateProbs( t + 1 );

	    gp_assert_eq( n_states, frSpace.getNStates() );
	    gp_assert_eq( n_states, toSpace.getNStates() );

	    const size_t locOffset = t * pCache.locMultiplier;

	    for ( HiddenStateSpace::StateIdxType frIdx = n_states ; frIdx-- != 0 ; )
		{
		const size_t stOffset = frIdx * pCache.stMultiplier;
		for ( HiddenStateSpace::StateIdxType toIdx = n_states ; toIdx-- != 0 ; )
		    probs[ locOffset + stOffset + toIdx ] =
			computeProb( frSpace, frIdx, toSpace, toIdx, f, g, getMu() );
		}
	    }

    #endif

    }



//-----------------------------------------------------------------------------
// Destructor
//-----------------------------------------------------------------------------

TransProbCache::~TransProbCache()
    {
    #if TPC_CACHE_MODEL == TPC_BIG_CACHE
	delete[] probs;
    #endif
    }



void TransProbCache::setRho( double nv )
    {
    // Needs conversion to TP_BIG_CACHE
    rho = nv;
    computeFactors();
    }


void TransProbCache::setMu( const MuType & nv )
    {
    mu = &nv;
    muChanged();
    }


void TransProbCache::muChanged()
    {
    // Needs conversion to TP_BIG_CACHE
    gp_assert_eq( pedigree.getNFounders(), getMu().size() );
    gp_assert_eq( pedigree.getK(), getMu()[0].size() );
    }



//=============================================================================
// Transition probabilities:
//=============================================================================

//-----------------------------------------------------------------------------
// computeF() [static]
/// \f$f = (1 - e^{-\rho x})\f$
//-----------------------------------------------------------------------------

inline double TransProbCache::computeF( const SimpleLocusArray & loci, SLocIdxType t, double rho )
    {
    // The distance from one chromosome to another is effectively infinity;
    // therefore e^-x is 0.  We make an explicit check for this here.

    const SimpleLocus & locTPlusOne = loci[ t + 1 ];

    return locTPlusOne.startsNewChromosome()			  ?
	1.0							  :
	1.0 - exp( - rho * locTPlusOne.getDistance().inMorgans() );
    }



//-----------------------------------------------------------------------------
// computeG() [static]
/// \f$g = e^-x\f$
//-----------------------------------------------------------------------------

inline double TransProbCache::computeG( const SimpleLocusArray & loci, SLocIdxType t )
    {
    // The distance from one chromosome to another is effectively infinity;
    // therefore e^-x is 0.  We make an explicit check for this here.

    const SimpleLocus & locTPlusOne = loci[ t + 1 ];

    return locTPlusOne.startsNewChromosome()		?
	0.0						:
	exp( - locTPlusOne.getDistance().inMorgans() )	;
    }



#if 0
    //-----------------------------------------------------------------------------
    // bitsInULong() [static helper]
    // Could use for more sophisticated cache strategy
    //-----------------------------------------------------------------------------
    static inline size_t bitsInULong( unsigned long x )
	{
	#define T_SIZE (1U << 16)
	static size_t tbl16[ T_SIZE ];
	static bool initialized = false;
	if ( ! initialized )
	    {
	    for ( size_t i = T_SIZE ; i-- != 0 ; )
		tbl16 ofoo here
	    }
	}
#endif



//-----------------------------------------------------------------------------
// computeProb() [static]
//-----------------------------------------------------------------------------

double TransProbCache::computeProb( const HiddenStateSpace & frSpace, HiddenStateSpace::StateIdxType frIdx,
				    const HiddenStateSpace & toSpace, HiddenStateSpace::StateIdxType toIdx,
				    double f, double g, const MuType & _mu )
    {
    const Pedigree & ped = frSpace.getPed();
    gp_assert( &ped == &(toSpace.getPed()) );

    const HiddenStateSpace::State & A = frSpace.stateAtIdx( frIdx );
    const HiddenStateSpace::State & B = toSpace.stateAtIdx( toIdx );

    gp_assert( A.av.size() == B.av.size() ); // DEBUG

    #if 0 // overkill
	gp_assert( &A.av.getPed() == &ped );
	gp_assert( &A.iv.getPed() == &ped );
	gp_assert( &B.av.getPed() == &ped );
	gp_assert( &B.iv.getPed() == &ped );
    #endif

    double rv = 1.0;

    const size_t n_meiosis = (ped.getNNonFndrs() << 1);

    unsigned long ivTransBits = A.iv.to_ulong() ^ B.iv.to_ulong();
    for ( size_t ctr = n_meiosis ; ctr-- != 0 ; )
	{
	const bool notEqual = (ivTransBits & 1);
	ivTransBits >>= 1;
	rv *= (notEqual ? (g + 1) : (1 - g)) * 0.5;
	#if DEBUG_PRINT_TP_DETAIL
	    printf( " (%lu,%lu) IV%lu=%d: %lf %lf\n", frIdx, toIdx, ctr, notEqual, (notEqual ? (g + 1) : (1 - g)) * 0.5, rv );
	#endif
	}

    for ( AncestryVector::IdxType avIdx = A.av.size() ; avIdx-- != 0 ; )
	{
	const PopIdx i = A.av.at( avIdx );
	const PopIdx j = B.av.at( avIdx );

	// The ancestry-vector is indexed on founder-gamete; mu (the ancestry
	// proportions) is indexed on founder (since we are assuming that the
	// ancestry proportions are the same for both of a founder's gametes) To
	// find mu, we need to translate the founder-gamete index into a
	// founder-index.
	const Pedigree::FounderIdx fIdx = avIdx >> 1; // NOTE *X117*!!!!
	double factor = (1 - f) * _mu[fIdx][j];
	if ( i == j )
	    factor += f;
	rv *= factor;
	#if DEBUG_PRINT_TP_DETAIL
	    printf( " (%lu,%lu) AV%lu = %lu,%lu{%d}: %lf %lf\n", frIdx, toIdx, avIdx, i, j, (i == j), factor, rv );
	#endif
	}

    #if DEBUG_PRINT_TP_DETAIL
	printf( " (%lu,%lu) trans-prob=%lf\n", frIdx, toIdx, rv );
    #endif

    return rv;
    }



//-----------------------------------------------------------------------------
// getProb()
//-----------------------------------------------------------------------------

#if TPC_CACHE_MODEL == TPC_NO_CACHE

    double TransProbCache::getProb( SLocIdxType t,
				    HiddenStateSpace::StateIdxType i,
				    HiddenStateSpace::StateIdxType j ) const
	{
	// t must be less than T-1 since we are computing the probability from t to t+1:
	gp_assert_lt( t, (loci.size() - 1) );

	const SLocFacts & fact = factors[ t ];

	return computeProb( pedigree.getStateProbs( t	), i,
			    pedigree.getStateProbs( t+1 ), j,
			    fact.f, fact.g, getMu() );
	}

#endif



//-----------------------------------------------------------------------------
// PedCache::getNProbs()
//-----------------------------------------------------------------------------

#if TPC_CACHE_MODEL == TPC_BIG_CACHE
    size_t TransProbCache::PedCache::getNProbs() const
	{
	#if AGGRESSIVE_RANGE_CHECK
	    return n_probs;
	#else
	    return 0;//locMultiplier * T_minus_1;
	#endif
	}
#endif



} // ---- end namespace genepi
