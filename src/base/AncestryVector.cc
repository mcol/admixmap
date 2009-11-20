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
/// \file AncestryVector.cc
/// Implementation of the AncestryVector class
//=============================================================================

#include "AncestryVector.h"


#if AV_OSTREAM
    #include <iostream>
#endif


#define MULTIPLE_PARMS	0

#if MULTIPLE_PARMS && defined(_OPENMP)
    #define THREAD_PRIVATE_CACHE	0
    #if THREAD_PRIVATE_CACHE
	#include <omp.h>
    #endif
#endif





namespace genepi { // ----



//=============================================================================
// Index cache:
//=============================================================================

static inline unsigned long K_to_the_F( unsigned int K, unsigned int F )
    {
    if	    ( K == 2 )
	return (1U << F);
    else if ( K == 4 )
	return (1U << (F << 1));
    else if ( F == 2 )
	return (K * K);
    else
	{
	unsigned long rv = 1;
	while ( F-- != 0 )
	    rv *= K;
	return rv;
	}
    }


namespace
{  // Begin anonymous namespace to give IdxCache internal linkage


//-----------------------------------------------------------------------------
/// The IdxCache data-structure caches the "space" of ancestry vectors for a
/// given pedigree structure, allowing an integer index into the space to
/// quickly be translated into an AncestryVector by table-lookup.  This is
/// important for performance because it is a very frequent operation in the
/// simulations and on-the-fly generation of the AV from the index, while
/// possible, is non-trivial.  The number of vectors in the space is K^F where K
/// is the number of populations and F is the number of founder-gametes in the
/// Pedigree.  Since the larger space of AVs of a larger pedigree with a larger
/// number of founder-gametes is a superset of smaller spaces, we can just keep
/// one table with as many elements as required for the largest possible number
/// of founder-gametes, and use the lower-indexed elements for smaller
/// pedigrees.  It is a "static linkage" class which is only used internally by
/// the implementation of AncestryVector.  While it could (and should) be used
/// in a "AncestryVectorSpace" class thus allowing multiple instances, as
/// currently used by AncestryVector there is only one global table.  This is
/// nonetheless relatively thread-safe because it is computed once at
/// initialization time (which is protected as a critical section for OpenMP)
/// and thereafter used in a read-only manner.  The initial computation is done
/// by set_parms(), which then fixes the number of populations and maximum
/// number of founder gametes to single global values.
//-----------------------------------------------------------------------------

class IdxCache
    {
    public:
	typedef AncestryVector::DType DType;

    private:

	PopIdx K; ///< Number of populations
	struct El
	    {
	    size_t  F	  ; ///< Number of founder gametes
	    size_t  n_els ; ///< K^F, size of val vect, cached for range-checking
	    DType * val   ;

	    El() : val ( 0 ) { }
	    El( const El & rhs ) : F( rhs.F ), n_els( rhs.n_els ), val( rhs.val ) { const_cast<El&>(rhs).val = 0; }
	    El& operator=( const El & rhs )
		{ F = rhs.F; val = rhs.val; const_cast<El&>(rhs).val = 0; return *this; }
	    ~El() { delete[] val; }

	    void generate( size_t K, size_t F );
	    };

	#if MULTIPLE_PARMS
	    cvector<El> els;
	#else
	    El el; ///< One element with F = maximum-F
	#endif

    public:

	#if MULTIPLE_PARMS
	    IdxCache( PopIdx _K ) : K( _K ) {}
	    void reserve( size_t max_F ) { els.reserve( max_F ); }
	#else
	    IdxCache( PopIdx _K, size_t max_F ) : K( _K )
		{
		el.generate( _K, max_F );
		}

	    /// Reverse-lookup:
	    const DType & getByIndex( size_t idx ) const { return el.val[idx]; }
	#endif


	PopIdx getK() const { return K; }

	const DType & lookup( PopIdx _K, size_t F, unsigned long index ) const;
    };



void IdxCache::El::generate( size_t _K, size_t _F )
    {
    F = _F;
    const size_t size = K_to_the_F( _K, _F );
    n_els = size;

    val = new DType[ size ];
    memset( val, '\0', size );

    #if 0
        DType * last = val + (size-1); // **** DEBUG ****
    #endif
    DType * ptr = val;
    for ( size_t x = size ; --x != 0 ; )
        {
        DType & prev = *ptr;
        DType & next = *++ptr;
        next = prev;
        #if 0
            gp_assert_eq( prev.size(), F ); // **** DEBUG ****
            gp_assert_eq( next.size(), F ); // **** DEBUG ****
            fprintf( stderr, "-> %lu %lu %p %p %p\n", size, x, &prev, &next, last ); // **** DEBUG ****
        #endif

        // ---- begin "increment" operation ----
        size_t curEl = 0;
        PopIdx val = next.at_unsafe( curEl );

        while ( ++val == _K )
            {
            next.setAt_unsafe( curEl++, 0 );
            gp_assert_lt( curEl, F );
            val = next.at_unsafe( curEl );
            }
        next.setAt_unsafe( curEl, val );
        // ---- end "increment" operation ----
        }
    }



#if MULTIPLE_PARMS

    // This implementation can maintain separate caches for multiple values of K
    // and F.  Not currently used.

    const IdxCache::DType & IdxCache::lookup( PopIdx _K, size_t F, unsigned long index ) const
	{
	gp_assert_eq( _K, K );

	#if 0 // DEBUG, will be checked below for cheaper in one of two places (see *!*!*)
	    gp_assert_lt( index, K_to_the_F(_K,F) );
	#endif

	for ( std::vector<El>::iterator it = els.begin(); it != els.end(); ++it )
	    if ( F == it->F )
		{
		gp_assert_lt( index, it->n_els ); // (*!*!*)
		return it->val[ index ]; // **** RETURN HERE ****
		}



	// We didn't find an existing entry in the cache with number of
	// founder-gametes "F", so we'll create one and populate it.  Since this
	// involves modifying the (shared) cache, we mark this section as critical.

	const DType * rv;

	#if defined(_OPENMP) && (! THREAD_PRIVATE_CACHE)
	  #pragma omp critical
	    {
	#endif

	#define MOD_IN_PLACE 0
	#if MOD_IN_PLACE
	    els.resize( els.size() + 1 );
	    El & newEl = els.back();
	#else
	    El newEl;
	#endif

	newEl.generate( K, F );

	gp_assert_lt( index, size ); // (*!*!*)
	#if 0
	    gp_assert_eq( newEl.val[index].size(), F ); // **** DEBUG ****
	#endif
	rv = &( newEl.val[ index ] );

	#if ! MOD_IN_PLACE
	    els.push_back( newEl );
	#endif

	#if defined(_OPENMP) && (! THREAD_PRIVATE_CACHE)
	  } // End critical section
	#endif

	return *rv;
	}

#else

    const IdxCache::DType & IdxCache::lookup( PopIdx _K, size_t F, unsigned long index ) const
	{

	#if AGGRESSIVE_RANGE_CHECK
	    gp_assert_eq( _K, K );
	    gp_assert_le( F, el.F );
	    #if 0 // DEBUG, will be checked below for cheaper in one of two places (see *!*!*)
		gp_assert_lt( index, K_to_the_F(_K,F) );
	    #endif
	#else
	    if ( _K ) {;} // suppress compiler warning
	    if ( F  ) {;} // suppress compiler warning
	#endif

	return el.val[ index ];
	}

#endif


} // End anonymous namespace to give IdxCache internal linkage



//-----------------------------------------------------------------------------
// set_parms()
//-----------------------------------------------------------------------------

#if THREAD_PRIVATE_CACHE

    void AncestryVector::set_parms( PopIdx /*K*/, size_t /*maxF*/ )
	{
	// Explicitly turn off dynamic threads (to preserve the value of
	// cache within a thread from one parallel task to another) (is this
	// required?)
	omp_set_dynamic( false );
	}

#else

    static IdxCache * cache = 0;

    void AncestryVector::set_parms( PopIdx K, size_t maxF )
	{
	gp_assert( cache == 0 );
	cache = new IdxCache( K, maxF ); // Could use AV_MAX_FOUNDER_GAMETES
	#if ! THREAD_PRIVATE_CACHE
	    #if MULTIPLE_PARMS
		for ( size_t F = 1 ; F < 6 ; ++F )
		    cache->lookup( K, F, 0 );
	    #endif
	#endif
	}

#endif


static const IdxCache & get_cache()
    {
    #if THREAD_PRIVATE_CACHE
	#pragma omp threadprivate( cache )
	if ( cache == 0 )
	    cache = new IdxCache( K );
    #endif

    return *cache;
    }



//-----------------------------------------------------------------------------
// set_ulong()
//-----------------------------------------------------------------------------

void AncestryVector::set_ulong( unsigned long nv )
    {
    memcpy( data.bytes, get_cache().lookup(K, ped.getNFounderGametes(), nv).bytes, AV_MC_BYTES() );
    }



//-----------------------------------------------------------------------------
// to_ulong()
//-----------------------------------------------------------------------------

unsigned long AncestryVector::to_ulong() const
    {
    unsigned long rv = 0;

    if ( K == 2 )	// Optimized version for K==2:
	for ( FGIdx idx = size() ; idx-- != 0 ; )
	    rv = (rv << 1) + at(idx);

    else if ( K == 3 )	// Optimized version for K==3:
	for ( FGIdx idx = size() ; idx-- != 0 ; )
	    rv = (rv << 1) + rv + at(idx);

    else
	for ( FGIdx idx = size() ; idx-- != 0 ; )
	    rv = (rv * K) + at(idx);


    // We need to cache the limit here:
    //gp_assert_lt( rv, K^F );


    return rv;
    }



//-----------------------------------------------------------------------------
// Iterator
//-----------------------------------------------------------------------------

AncestryVector::Iterator::Iterator( const Pedigree & _ped ) :
	idx	 ( 0 ) ,
	idxLimit ( K_to_the_F( get_cache().getK(), _ped.getNFounderGametes() ) ) ,
	av	 ( _ped, get_cache().getK(), get_cache().getByIndex(0) )
    {
    gp_assert( ! MULTIPLE_PARMS );
    }


bool AncestryVector::Iterator::advance()
    {
    const bool rv = (++idx < idxLimit);

    if ( rv )
	memcpy( av.data.bytes, get_cache().getByIndex(idx).bytes, AV_MC_BYTES() );

    return rv;
    }



//-----------------------------------------------------------------------------
// Print an AncestryVector to an ostream
//-----------------------------------------------------------------------------

#if AV_OSTREAM

    static AVOutputStyle style = AV_NUMERIC;

    void setAVOutputStyle( AVOutputStyle nv )
	{
	style = nv;
	}

    inline static std::ostream & output( std::ostream & os, PopIdx pop )
	{
	#if 0
	    if ( style == AV_NUMERIC )
		return os << popVect[pop];
	    else
	#endif
		return os << size_t(pop);
	}

    std::ostream & output( std::ostream & os,
			   const AncestryVector & av,
			   AncestryVector::FGIdx nValid )
	{
	AncestryVector::FGIdx idx;

	gp_assert( av.size() != 0 );
	gp_assert( nValid <= av.size() );

	os << "AV(";
	for ( idx = 0; idx < nValid; ++idx )
	    output( os << ((idx==0)?"":","), av.at(idx) );
	for ( ; idx < av.size(); ++idx )
	    os << ((idx==0) ? "X" : ",X");
	os << ')';

	return os;
	}

    std::ostream & operator<<( std::ostream & os, const AncestryVector & av )
	{
	return output( os, av, av.size() );
	}

#endif // AV_OSTREAM



} // ---- end namespace genepi
