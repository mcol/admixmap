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
/// \file HiddenStateSpace.cc
/// Implementation of the HiddenStateSpace class.
//=============================================================================

#pragma implementation
#include "HiddenStateSpace.h"

#include <cmath>  // exp()
#include <limits> // numeric_limits::min()



namespace genepi { // ----



#if TRACK_UNVISITED_STATES
    const HiddenStateSpace::ProbType HiddenStateSpace::State::NOT_VISITED = std::numeric_limits<ProbType>::min();
#endif



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
	nNon0 ( 0			 ) ,
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
    nNon0 = 0			     ;
    probs = new ProbType [ aSize() ] ;
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
	#if TRACK_UNVISITED_STATES
	    probs[ idx ] = State::NOT_VISITED;
	#else
	    probs[ idx ] = 0.0;
	#endif
    }



//=============================================================================
// Iterator:
//=============================================================================


HiddenStateSpace::Iterator::Iterator( const HiddenStateSpace & sp ) :
	space	 ( sp	       ) ,
	av_it	 ( sp.getPed() ) ,
	iv_it	 ( sp.getPed() ) ,
	sIdx	 ( 0	       ) ,
	non0Idx	 ( 0	       ) ,
	finished ( false       )
    {
    if ( sp.getEProb( sIdx ) == 0.0 )
	advance();
    }



bool HiddenStateSpace::Iterator::advance()
    {

    #if ! HSS_AV_MOST_SIG
	#error Must reimplement HiddenSpaceState::Iterator::advance() for IV-most-significant order
    #endif


    if ( isFinished() )
	throw runtime_error(
	    string("Advance iterator when already past last state in pedigree ") +
	    space.getPed().getId() );


    if ( (space.getEProb(sIdx) != 0.0) &&
	    (++non0Idx == space.getNNon0()) )	// Saves skipping past the last segment of 0-EPs
	return (finished = true);		// **** RETURN HERE ****

    do
	{

	if ( iv_it.is_on_last_el() )
	    {
	    finished = ! av_it.advance();

	    if ( ! finished )
		{
		iv_it.reset();
		++sIdx;
		gp_assert_eq( sIdx, (av_it.to_ulong() * space.N_IVs) );
		}
	    }
	else
	    {
	    iv_it.advance();
	    ++sIdx;
	    }

	} while ( (! finished) && (space.getEProb(sIdx) == 0.0) );


    return finished;
    }



HiddenStateSpace::State HiddenStateSpace::Iterator::getState() const
    {
    gp_assert( ! finished );
    State rv = { av_it.getAV(), iv_it.getIV(), space.getEProb(sIdx) };
    return rv;
    }



HiddenStateSpace::State HiddenStateSpace::stateAtIdx( StateIdxType idx ) const
    {
    const size_t iv_idx = idx % N_IVs;
    const size_t av_idx = idx / N_IVs;
    State rv = { AncestryVector(getPed(),K), InheritanceVector(getPed()), getEProb(idx) };
    gp_assert( rv.av.size() != 0 );
    rv.av.set_ulong( av_idx );
    gp_assert( rv.av.size() != 0 );
    rv.iv.set_ulong( iv_idx );
    return rv;
    }



} // ---- end namespace genepi
