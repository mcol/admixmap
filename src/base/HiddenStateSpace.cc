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
/// Implementation of the genepi::HiddenStateSpace class.
//=============================================================================

#pragma implementation
#include "HiddenStateSpace.h"

#include <cmath>
#include <stdexcept> // for runtime_error()
#include <string>



namespace genepi { // ----



#if TRACK_UNVISITED_STATES
    #include <limits>
    const HiddenStateSpace::ProbType HiddenStateSpace::State::NOT_VISITED = std::numeric_limits<ProbType>::min();
#endif



//-----------------------------------------------------------------------------
/// Static helper function, returns K^F (optimized for certain common cases)
//-----------------------------------------------------------------------------

static inline size_t k_pow_f( PopIdx K , size_t F )
    {
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
	ped	   ( &_ped				    ) ,
	K	   ( _K					    ) ,
	N_IVs_X	   ( 1U << _ped.getNMeiosis(CHR_IS_X)	    ) , // 2^M
	N_IVs_notX ( 1U << _ped.getNMeiosis(CHR_IS_NOT_X)   ) , // 2^M
	N_AVs_X	   ( k_pow_f ( _K, _ped.getNFounderGametes(CHR_IS_X    ) ) ) , // K^F
	N_AVs_notX ( k_pow_f ( _K, _ped.getNFounderGametes(CHR_IS_NOT_X) ) ) , // K^F
	nNon0	   ( 0					    ) ,
	probs	   ( new ProbType [ aSize(CHR_IS_NOT_X) ]   )
    {
    }


HiddenStateSpace::HiddenStateSpace() :
	ped  ( 0 ) ,
	probs( 0 )
    {
    }


void HiddenStateSpace::init( const Pedigree & _ped, PopIdx _K )
    {

    gp_assert( ped == 0 );

    ped = &_ped;
    K	= _K;

    N_IVs_X    = 1U << _ped.getNMeiosis(CHR_IS_X    ); // 2^M
    N_IVs_notX = 1U << _ped.getNMeiosis(CHR_IS_NOT_X); // 2^M
    N_AVs_X    = k_pow_f ( _K, _ped.getNFounderGametes(CHR_IS_X	   ) ) ; // K^F
    N_AVs_notX = k_pow_f ( _K, _ped.getNFounderGametes(CHR_IS_NOT_X) ) ; // K^F

    nNon0 = 0;

    gp_assert( probs == 0 );
    probs = new ProbType [ aSize(CHR_IS_NOT_X) ] ;

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
    for ( size_t idx = aSize(CHR_IS_NOT_X) ; idx-- != 0 ; )
	#if TRACK_UNVISITED_STATES
	    probs[ idx ] = State::NOT_VISITED;
	#else
	    probs[ idx ] = 0.0;
	#endif
    }



//=============================================================================
// Iterator:
//=============================================================================


HiddenStateSpace::Iterator::Iterator( const HiddenStateSpace & sp, IsXChromType isX, bool sparse ) :
	space	 ( sp		    ) ,
	av_it	 ( sp.getPed(), isX ) ,
	iv_it	 ( sp.getPed(), isX ) ,
	sIdx	 ( 0		    ) ,
	non0Idx	 ( 0		    ) ,
	isXChrom ( isX		    ) ,
	skipNon0 ( sparse	    ) ,
	finished ( false	    )
    {
    if ( skipNon0 && sp.getEProb( sIdx ) == 0.0 )
	advance();
    }



HiddenStateSpace::Iterator::Iterator( const Iterator & rhs ) :
	space	 ( rhs.space	),
	av_it	 ( rhs.av_it	),
	iv_it	 ( rhs.iv_it	),
	sIdx	 ( rhs.sIdx	),
	non0Idx  ( rhs.non0Idx  ),
	finished ( rhs.finished )
    {
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


    if ( skipNon0 && (space.getEProb(sIdx) != 0.0) &&
	    (++non0Idx == space.getNNon0()) )	// Saves skipping past the last segment of 0-EPs
	return (finished = true);		// **** RETURN HERE ****

    do
	{

	if ( iv_it.is_on_last_el() )
	    {
	    finished = ! av_it.advance();

	    if ( ! finished )
		{
		++sIdx;
		iv_it.reset( isX() );
		gp_assert_eq( sIdx, (av_it.to_ulong() * space.N_IVs(isX())) );
		}
	    }
	else
	    {
	    iv_it.advance();
	    ++sIdx;
	    }

	} while ( skipNon0 && (! finished) && (space.getEProb(sIdx) == 0.0) );


    return finished;
    }



HiddenStateSpace::State HiddenStateSpace::Iterator::getState() const
    {
    gp_assert( ! finished );
    State rv = { av_it.getAV(), iv_it.getIV(), space.getEProb(sIdx) };
    return rv;
    }



HiddenStateSpace::State HiddenStateSpace::stateAtIdx( StateIdxType idx, IsXChromType isX ) const
    {
    const size_t iv_idx = idx % N_IVs( isX );
    const size_t av_idx = idx / N_IVs( isX );
    State rv = { AncestryVector(getPed(),K), InheritanceVector(getPed(),isX), getEProb(idx) };
    gp_assert( rv.av.size(isX) != 0 );
    rv.av.set_ulong( av_idx, isX );
    gp_assert( rv.av.size(isX) != 0 );
    rv.iv.set_ulong( iv_idx );
    return rv;
    }



} // ---- end namespace genepi
