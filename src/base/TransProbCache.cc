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
#define DEBUG_TPC_MEMORY	0


#include <cmath> // exp()


namespace genepi { // ----



//-----------------------------------------------------------------------------
// reComputeFactors() [protected]
// Recompute the cached "f" values (called initially and when rho changes).
//-----------------------------------------------------------------------------

void TransProbCache::reComputeFactors()
    {
    // No need to recompute g since it doesn not depend on rho
    for ( int t = factors.size() ; t-- != 0 ; )
	factors[ t ].f = computeF( loci, t, rho );
    }



//-----------------------------------------------------------------------------
// LocusTP::init()
//-----------------------------------------------------------------------------

#if TPC_CACHE_MODEL == TPC_BIG_CACHE

    void TransProbCache::LocusTP::init( size_t from_n_states, size_t to_n_states )
	{
	stMultiplier = to_n_states;

	const size_t npr = from_n_states * to_n_states;

	gp_assert( probs == 0 );
	probs = new CacheProbType[ npr ];

	#if AGGRESSIVE_RANGE_CHECK
	    // Keep the number of probabilities (array-size) for later range checks:
	    n_probs = npr;
	#endif
	}

#endif



//-----------------------------------------------------------------------------
//
// recomputeBigCache() [protected] [big-cache model only]
//
/// Recompute the cached transition probabilities for one whole pedigree
/// (i.e. one matrix for each locus).  Called when constructed and again
/// whenever rho or mu (aka theta) changes.
//
//-----------------------------------------------------------------------------

#if TPC_CACHE_MODEL == TPC_BIG_CACHE

    void TransProbCache::recomputeBigCache()
	{
	const size_t T_minus_1 = loci.size() - 1;

	// Let's try using the same code that caches the F and G factors for the
	// no-cache model, even though they don't need to be persistent for the
	// big-model (they take up so much less space than the TPs themselves, it's
	// insignificant).
	#define BIG_MODEL_TEMP_FACTORS 0
	#if BIG_MODEL_TEMP_FACTORS
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
	#endif


	for ( SLocIdxType t = T_minus_1 ; t-- != 0 ; )
	    {

	    LocusTP &	 lProbs	      = probs[ t ];
	    const size_t stMultiplier = lProbs.stMultiplier;

	    const HiddenStateSpace & frSpace = pedigree.getStateProbs( t     );
	    const HiddenStateSpace & toSpace = pedigree.getStateProbs( t + 1 );

	    #if AGGRESSIVE_RANGE_CHECK
		gp_assert_eq( lProbs.n_probs, frSpace.getNNon0() * toSpace.getNNon0() );
		gp_assert_eq( stMultiplier, toSpace.getNNon0() );
	    #endif

	    #if BIG_MODEL_TEMP_FACTORS
		const double f = fa[t];
		const double g = ga[t];
	    #endif

	    for ( HiddenStateSpace::Iterator fr_it( frSpace ) ; fr_it ; ++fr_it )
		{
		const size_t	  stOffset = fr_it.getNon0Index() * stMultiplier;
		const SLocFacts & fact	   = factors[ t ];
		for ( HiddenStateSpace::Iterator to_it( toSpace ) ; to_it ; ++to_it )
		    // lProbs.getProbs( fr_it.getNon0Index(), to_it.getNon0Index() ) =
		    lProbs.probs[ stOffset + to_it.getNon0Index() ] =
			computeProb( fr_it, to_it, fact, getMu() );
		}
	    }

	cacheIsDirty = false;

	}


#endif // TPC_CACHE_MODEL == TPC_BIG_CACHE



//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------

TransProbCache::TransProbCache( const Pedigree & _pedigree, double _rho, const MuType & _mu ) :
	#if TPC_CACHE_MODEL == TPC_BIG_CACHE
	    cacheIsDirty( true ) ,
	#endif
	pedigree ( _pedigree		) ,
	loci	 ( _pedigree.getSLoci() ) ,
	rho	 ( _rho			) ,
	mu	 ( &_mu			)
    {

    #if (TPC_CACHE_MODEL == TPC_NO_CACHE) || ((TPC_CACHE_MODEL == TPC_BIG_CACHE) && (! BIG_MODEL_TEMP_FACTORS))

	muChanged();

	gp_assert( loci.size() != 0 ); // Unnecessary, but a bug did manifest itself this way once


	// Pre-compute and cache the table of factors that are based on the
	// number of non-zero bits in the XOR of the from and to IVs.  As
	// explained in the header-file, the double-lookup is condensed to a
	// single table-lookup, but when computing the table, we first compute a
	// table (preCompTab) containing the values of the factor, indexed on
	// the number of non-0 bits; then create the lookup table (gtab) that is
	// indexed on the value of the XOR of the bitmaps of the two IVs.
	//
	// The value of the g-factor for a transition in which there are n
	// differing bits between the two IVs is (here "^" means "to the power
	// of"):
	//	[0.5*(1-g)]^n * [0.5*(g+1)]^k
	// Where k is the number of bits which do _not_ differ between the two
	// IVs, and is therefore equal to (M-n) where M is the total number of
	// meiosis in the pedigree at that locus.


	factors.resize( loci.size() - 1 );
	for ( size_t t = factors.size() ; t-- != 0 ; )
	    {

	    const double g = computeG( loci, t );
	    factors[ t ].g = g;

	    cvector<CachedIVFactor> & gtab = factors[ t ].iv_factors;
	    const size_t n_meiosis = _pedigree.getNMeiosis();
	    gtab.resize( 1 << n_meiosis ); // 2^n_meiosis

	    if ( n_meiosis == 0 )
		gtab[0] = 1.0; // When the pedigree contains no meiosis, "null" IV
	    else
		{

		// The number of non-0 bis can be from 0 to M, so we need a
		// lookup-table of size M+1:
		double preCompTab[ n_meiosis + 1 ];

		// First, compute the first half of the factor in each slot,
		// [0.5*(1-g)]^n
		const double oneMinusG = (1 - g) * 0.5;
		double oneMinusGTerm = 1.0;
		for ( size_t ctr = 0 ; ctr <= n_meiosis ; ++ctr )
		    {
		    preCompTab[ ctr ] = oneMinusGTerm;
		    oneMinusGTerm *= oneMinusG;
		    }

		// Next multiply the value in each slot by the second half of
		// the factor, [0.5*(g+1)]^k
		const double gPlusOne = (g + 1) * 0.5;
		double gPlusOneTerm = gPlusOne;
		for ( size_t ctr = n_meiosis ; ctr-- != 0 ; )
		    {
		    preCompTab[ ctr ] *= gPlusOneTerm;
		    gPlusOneTerm *= gPlusOne;
		    }

		// Now fill in the values in the main table (gtab), indexed on
		// the result of the XOR, with the corresponding values.
		for ( size_t idx = 0 ; idx < gtab.size() ; ++idx )
		    {
		    unsigned int n_non_zero_bits = 0;
		    size_t ivTransBits = idx;
		    for ( size_t ctr = n_meiosis ; ctr-- != 0 ; )
			{
			n_non_zero_bits += (ivTransBits & 1);
			ivTransBits >>= 1;
			}
		    gp_assert_le( n_non_zero_bits, n_meiosis );
		    gtab[ idx ] = preCompTab[ n_non_zero_bits ];
		    }

		}
	    }

	reComputeFactors();

    #endif



    #if TPC_CACHE_MODEL == TPC_BIG_CACHE

	#if DEBUG_TPC_MEMORY
	    size_t tot_nprob = 0;
	    size_t tot_mem   = 0;
	#endif

	const SLocIdxType T_minus_1 = loci.size() - 1;

	probs.resize( T_minus_1 );

	for ( SLocIdxType loc = T_minus_1 ; loc-- != 0 ; )
	    {

	    // We no longer assume that all loci have the same state-space size,
	    // i.e. we only allocate for non-0 states.
	    const size_t from_n_states = pedigree.getStateProbs( loc	 ).getNNon0();
	    const size_t to_n_states   = pedigree.getStateProbs( loc + 1 ).getNNon0();

	    probs[ loc ].init( from_n_states, to_n_states );

	    #if DEBUG_TPC_MEMORY
		#if 0
		  if ( pedigree.getId() == "H010007" )
		     fprintf( stderr, "  TPC for %s fr-loc %zd %s: from-nst=%zd/to-nst=%zd n-tprobs=%zd -> %zd bytes\n",
			pedigree.getId().c_str(), loc, loci[loc].getName().c_str(),
			from_n_states, to_n_states, probs[loc].n_probs, probs[loc].n_probs*sizeof(CacheProbType) );
		#endif
		tot_nprob += probs[loc].n_probs;
		tot_mem += probs[loc].n_probs*sizeof(CacheProbType);
	    #endif
	    }

	#if DEBUG_TPC_MEMORY
	    fprintf( stderr, "Construct TPC for %s (%zd members): n-tprobs=%zd -> %zd bytes\n",
			pedigree.getId().c_str(), pedigree.getNMembers(), tot_nprob, tot_mem );
	    #if 0 && (_BSD_SOURCE || _SVID_SOURCE || _XOPEN_SOURCE >= 500)
		const void * const brkVal = sbrk(0);
		fprintf( stderr, "  break=%p\n", brkVal );
	    #endif
	#endif

    #endif

    }



//-----------------------------------------------------------------------------
// Destructor
//-----------------------------------------------------------------------------

TransProbCache::~TransProbCache()
    {
    }



void TransProbCache::setRho( double nv )
    {
    // Needs conversion to TP_BIG_CACHE
    if ( nv != rho )
	{
	rho = nv;
	reComputeFactors();
	#if (TPC_CACHE_MODEL == TPC_BIG_CACHE)
	    cacheIsDirty = true;
	#endif
	}
    }


void TransProbCache::setMu( const MuType & nv )
    {
    mu = &nv;
    muChanged();
    }


void TransProbCache::muChanged()
    {
    gp_assert_eq( pedigree.getNFounders(), getMu().size() );
    gp_assert_eq( pedigree.getK(), getMu()[0].size() );
    #if (TPC_CACHE_MODEL == TPC_BIG_CACHE)
	cacheIsDirty = true;
    #endif
    }



//=============================================================================
// Transition probabilities:
//=============================================================================

//-----------------------------------------------------------------------------
// computeF() [static]
/// \f$f = e^{-\rho x}\f$
//-----------------------------------------------------------------------------

double TransProbCache::computeF( const SimpleLocusArray & loci, SLocIdxType t, double rho )
    {
    // The distance from one chromosome to another is effectively infinity;
    // therefore e^-x is 0.  We make an explicit check for this here.

    const SimpleLocus & locTPlusOne = loci[ t + 1 ];

    return locTPlusOne.startsNewChromosome()		     ?
	0.0						     :
	exp( - rho * locTPlusOne.getDistance().inMorgans() ) ;
    }



//-----------------------------------------------------------------------------
// computeG() [static]
/// \f$g = e^-x\f$
//-----------------------------------------------------------------------------

double TransProbCache::computeG( const SimpleLocusArray & loci, SLocIdxType t )
    {
    // The distance from one chromosome to another is effectively infinity;
    // therefore e^-x is 0.  We make an explicit check for this here.

    const SimpleLocus & locTPlusOne = loci[ t + 1 ];

    return locTPlusOne.startsNewChromosome()		?
	0.0						:
	exp( - locTPlusOne.getDistance().inMorgans() )	;
    }



//-----------------------------------------------------------------------------
// computeProb() [static]
//-----------------------------------------------------------------------------

double TransProbCache::computeProb( const HiddenStateSpace::Iterator & frState,
				    const HiddenStateSpace::Iterator & toState,
				    const SLocFacts & facts, const MuType & _mu )
    {
    const AncestryVector & frAV = frState.getAV();
    const AncestryVector & toAV = toState.getAV();

    const Pedigree & ped = frState.getSpace().getPed();
    gp_assert( &(frState.getSpace().getPed()) == &(toState.getSpace().getPed()) );

    gp_assert_eq( frAV.size(), toAV.size() ); // DEBUG

    double rv = facts.iv_factors[ frState.getIV() ^ toState.getIV() ];

    const double f = facts.f;
    const double one_minus_f = 1 - f;

    #if DEBUG_PRINT_TP_DETAIL
	const unsigned long ivXor = frState.getIV() ^ toState.getIV();
	cout << "Compute trans-prob from-state=" << frState.getNon0Index() << ", to-state=" << toState.getNon0Index()
	     << ", from-IV=" << frState.getIV() << ", to-IV=" << toState.getIV() << ", xor=" << ivXor << ", f=" << f << '\n';
    #endif

    for ( Pedigree::FGameteIdx fgIdx = ped.getNFounderGametes() ; fgIdx-- != 0 ; )
	{

	Pedigree::GameteType whichOne;
	const Pedigree::FounderIdx fIdx = ped.founderOfGameteIdx( fgIdx, whichOne ); // , t.isXChromosome() );

	// The ancestry-vector is indexed on founder-gamete; mu (the ancestry
	// proportions) is indexed on founder (since we are assuming that the
	// ancestry proportions are the same for both of a founder's gametes) To
	// find mu, we need to translate the founder-gamete index into a
	// founder-index.
	const cvector<double> & mu_of_fIdx = _mu[ fIdx ];

	const PopIdx frAncestry = frAV.at( fgIdx );
	const PopIdx toAncestry = toAV.at( fgIdx );

	double factor = one_minus_f * mu_of_fIdx[ toAncestry ];
	if ( frAncestry == toAncestry )
	    factor += f;
	rv *= factor;

	#if DEBUG_PRINT_TP_DETAIL
	    printf( "  fndr:%zu(gmete:%zu)  = %zu,%zu: mu=%.9lf %.9lf %.9lf\n", fIdx, fgIdx, frAncestry, toAncestry, mu_of_fIdx[toAncestry], factor, rv );
	#endif

	}

    #if DEBUG_PRINT_TP_DETAIL
	printf( " (%lu,%lu) trans-prob=%lf\n", frState.getNon0Index(), toState.getNon0Index(), rv );
    #endif

    return rv;
    }



#if TPC_CACHE_MODEL == TPC_BIG_CACHE

    //-----------------------------------------------------------------------------
    // LocusTP::debugOut()
    //-----------------------------------------------------------------------------

    void TransProbCache::LocusTP::debugOut(
				    std::ostream &	os    ,
				    const Pedigree &	ped   ,
				    SLocIdxType		frLoc ) const
	{
	os << "Transition-probabilities for ped #" << ped.getMyNumber() << " (" << ped.getId()
	    << ") from locus #" << frLoc << " to " << (frLoc+1) << ":\n";

	const HiddenStateSpace & frSpace = ped.getStateProbs( frLoc );
	const HiddenStateSpace & toSpace = ped.getStateProbs( frLoc + 1 );

	for ( HiddenStateSpace::Iterator to_it(toSpace) ; to_it ; ++to_it )
	    os << '\t' << to_it.getOverallIndex() << '=' << to_it.getNon0Index();
	os << '\n';

	for ( HiddenStateSpace::Iterator fr_it(frSpace) ; fr_it ; ++fr_it )
	    {
	    os << fr_it.getOverallIndex() << '=' << fr_it.getNon0Index();
	    for ( HiddenStateSpace::Iterator to_it(toSpace) ; to_it ; ++to_it )
		os << '\t' << getProb(fr_it.getNon0Index(),to_it.getNon0Index());
	    os << '\n';
	    }
	}



    //-----------------------------------------------------------------------------
    // debugOut( std::ostream & ostr, SLocIdxType frLoc ) const;
    //-----------------------------------------------------------------------------

    void TransProbCache::debugOut( std::ostream & ostr, SLocIdxType frLoc ) const
	{
	probs[frLoc].debugOut( ostr, pedigree, frLoc );
	}

#endif



} // ---- end namespace genepi
