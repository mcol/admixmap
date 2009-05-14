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
/// \file HiddenStateSpace.cc
/// Implementation of the HiddenStateSpace class.
//=============================================================================

#pragma implementation
#include "HiddenStateSpace.h"

#include <cmath> // exp()


#define DEBUG_DENSITY	0



namespace genepi { // ----



//-----------------------------------------------------------------------------
/// Static helper function, returns K^F (optimized for certain common cases)
//-----------------------------------------------------------------------------

static inline size_t k_pow_f( PopIdx K , const Pedigree & ped )
    {
    size_t F = ped.getNFounderGametes();

    size_t rv;

    if ( K == 2 )
	rv = (size_t(1) << F);
    else if ( K == 4 )
	rv = (size_t(1) << (F<<1));
    else
	{
	rv = 1;
	while ( F-- != 0 )
	    rv *= K;
	}

    return rv;
    }



//-----------------------------------------------------------------------------
// Constructors:
//-----------------------------------------------------------------------------

HiddenStateSpace::HiddenStateSpace( const Pedigree & _ped, PopIdx _K ) :
	ped   ( &_ped			 ) ,
	K     ( _K			 ) ,
	N_IVs ( 1U << _ped.getNMeiosis() ) , // 2^M
	N_AVs ( k_pow_f( _K, _ped )	 ) , // K^F
	probs ( new ProbType [ aSize() ] )
    {
    }


HiddenStateSpace::HiddenStateSpace() :
	ped( 0 )
    {
    }


void HiddenStateSpace::init( const Pedigree & _ped, PopIdx _K )
    {
    ped	  = &_ped		     ;
    K	  = _K			     ;
    N_IVs = 1U << _ped.getNMeiosis() ; // 2^M
    N_AVs = k_pow_f( _K, _ped )	     ; // K^F
    probs = new ProbType [ aSize() ] ;
    #if DEBUG_DENSITY
	if ( aSize() > 70000 )
	    fprintf( stderr, "Allocate(%s): %lu = %lu\n", _ped.getId().c_str(),
		aSize(), aSize()*sizeof(ProbType) );
    #endif
    }



//-----------------------------------------------------------------------------
// Destructor:
//-----------------------------------------------------------------------------

HiddenStateSpace::~HiddenStateSpace()
    {
    delete[] probs;
    }



//-----------------------------------------------------------------------------
// resetEmProbsToZero()
//-----------------------------------------------------------------------------

void HiddenStateSpace::resetEmProbsToZero()
    {
    for ( size_t idx = aSize() ; idx-- != 0 ; )
	probs[ idx ] = 0.0;
    }



//=============================================================================
// Iterator:
//=============================================================================


HiddenStateSpace::Iterator::Iterator( const HiddenStateSpace & sp ) :
	space	 ( sp ) ,
	av	 ( sp.getPed(), sp.getK() ) ,
	iv	 ( sp.getPed() ) ,
	idx	 ( 0 ) ,
	finished ( false )
    {
    for ( size_t aIdx = av.size() ; aIdx-- != 0 ; )
	av.setAt( aIdx, 0 );
    iv.set_ulong( 0 );

    if ( sp.probs[ idx ] == 0.0 )
	advance();
    }



bool HiddenStateSpace::Iterator::advance()
    {

    #if ! HSS_AV_MOST_SIG
	#error Must reimplement HiddenSpaceState::Iterator::advance() for IV-most-significant order
    #endif


    do
	{
	// In theory this is vulnerable to overflow, but in practice we would use a
	// different (sparse) implementation of HiddenStateSpace for any pedigree
	// with enough members to even get close:
	const unsigned long nextIV = iv.to_ulong() + 1;

	if ( nextIV == space.N_IVs )
	    {
	    finished = true;
	    iv.set_ulong( 0 );

	    for ( size_t aIdx = 0 ; aIdx < av.size() ; ++aIdx )
		{
		PopIdx p = av[aIdx] + 1;
		if ( p == space.getK() )
		    {
		    p = 0;
		    av.setAt( aIdx, p );
		    }
		else
		    {
		    av.setAt( aIdx, p );
		    finished = false;
		    break;
		    }
		}

	    if ( ! finished )
		{
		++idx;
		gp_assert( idx == (av.to_ulong() * space.N_IVs) );
		}
	    }
	else
	    {
	    iv.set_ulong( nextIV );
	    ++idx;
	    //gp_assert( idx == (av.to_ulong() * space.N_IVs) + nextIV );
	    }

	} while ( (! finished) && (space.probs[ idx ] == 0.0) );

    return finished;
    }



HiddenStateSpace::State HiddenStateSpace::Iterator::getState() const
    {
    gp_assert( ! finished );

    const State rv = { av, iv, space.probs[idx] };
    return rv;
    }



HiddenStateSpace::State HiddenStateSpace::stateAtIdx( StateIdxType idx ) const
    {
    const size_t iv_idx = idx % N_IVs;
    const size_t av_idx = idx / N_IVs;
    State rv = { AncestryVector(getPed(),K), InheritanceVector(getPed()), probs[idx] };
    gp_assert( rv.av.size() != 0 );
    rv.av.set_ulong( av_idx );
    gp_assert( rv.av.size() != 0 );
    rv.iv.set_ulong( iv_idx );
    return rv;
    }



//-----------------------------------------------------------------------------
// Compatibility methods for "old" (individual-based) HMM:
//-----------------------------------------------------------------------------

void HiddenStateSpace::lambdaAsArrayOfDouble( double * lambda ) const
    {
    StateIdxType limit = getNStates();
    for ( StateIdxType idx = 0 ; idx < limit ; ++idx )
	lambda[ idx ] = probs[ idx ];
    }



} // ---- end namespace genepi
