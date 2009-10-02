//=============================================================================
//
// Copyright (C) 2009  David D. Favro  gpl-copyright@meta-dynamic.com
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



//-----------------------------------------------------------------------------
// Index cache:
//-----------------------------------------------------------------------------

static unsigned long K_to_the_F( unsigned int K, unsigned int F )
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


// "Static" class used internally to the implementation of AncestryVector:
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
	#endif


	const DType & lookup( PopIdx _K, size_t F, unsigned long index );
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

    const IdxCache::DType & IdxCache::lookup( PopIdx _K, size_t F, unsigned long index )
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

    const IdxCache::DType & IdxCache::lookup( PopIdx _K, size_t F, unsigned long index )
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
// set_ulong()
//-----------------------------------------------------------------------------

#if defined(_OPENMP)

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

#endif


void AncestryVector::set_ulong( unsigned long nv )
    {
    #if ! defined(_OPENMP)
	static IdxCache _cache( K, AV_MAX_FOUNDER_GAMETES );
	#define cache (&_cache)
    #elif THREAD_PRIVATE_CACHE
	#pragma omp threadprivate( cache )
	if ( cache == 0 )
	    cache = new IdxCache( K );
    #endif

    memcpy( data.bytes, cache->lookup(K, ped.getNFounderGametes(), nv).bytes, AV_MC_BYTES() );
    }



//-----------------------------------------------------------------------------
// to_ulong()
//-----------------------------------------------------------------------------

unsigned long AncestryVector::to_ulong() const
    {
    unsigned long rv = 0;

    if ( K == 2 )	// Optimized version for K==2:
	for ( IdxType idx = size() ; idx-- != 0 ; )
	    rv = (rv << 1) + at(idx);

    else if ( K == 3 )	// Optimized version for K==3:
	for ( IdxType idx = size() ; idx-- != 0 ; )
	    rv = (rv << 1) + rv + at(idx);

    else
	for ( IdxType idx = size() ; idx-- != 0 ; )
	    rv = (rv * K) + at(idx);


    // We need to cache the limit here:
    //gp_assert_lt( rv, K^F );


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

    std::ostream & operator<<( std::ostream & os, const AncestryVector & av )
	{
	gp_assert( av.size() != 0 );

	os << "AV(";
	const AncestryVector::IdxType limit = av.size() - 1;
	for ( AncestryVector::IdxType idx = 0 ; idx < limit ; ++idx )
	    output( os, av.at(idx) ) << ',';
	os << av.at(limit) << ')';

	return os;
	}

#endif // AV_OSTREAM



} // ---- end namespace genepi
