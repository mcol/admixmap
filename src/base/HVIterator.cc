//=============================================================================
//
// Copyright (C) 2009  David D. Favro
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License, version 2 or later, as published by
// the Free Software Foundation.
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
/// \file HVIterator.cc
/// Implementation of the genepi::HVIterator class.
//=============================================================================


#include "HVIterator.h"



namespace genepi { // ----




//-----------------------------------------------------------------------------
// Constructors
//-----------------------------------------------------------------------------

HVIterator::HVIterator( const Pedigree & ped, IsXChromType isX ) :
	n_ancestries( ped.getNFounderGametes(isX) ) ,
	n_meiosis   ( ped.getNMeiosis( isX )	  ) ,
	K	    ( ped.getK()		  ) ,
	cur_ancestry( 0 ) ,
	cur_meiosis ( 0 )
    {
    }



HVIterator::HVIterator( const HVIterator & rhs ) :
	n_ancestries( rhs.n_ancestries ),
	n_meiosis   ( rhs.n_meiosis    ),
	K	    ( rhs.K	       ),
	cur_ancestry( rhs.cur_ancestry ),
	cur_meiosis ( rhs.cur_meiosis  )
    {
    }



//-----------------------------------------------------------------------------
// advance()
//-----------------------------------------------------------------------------

void HVIterator::advance()
    {
    gp_assert( ! isFinished() ); // Don't allow advance after finished

    // NOTE *1*:
    if ( isOnMeiosis() )
	++cur_meiosis;
    else
	++cur_ancestry;
    }



//-----------------------------------------------------------------------------
// nValues()
//-----------------------------------------------------------------------------

int HVIterator::nValues() const
    {
    return (isOnAncestry() ? K : 2);
    }



} // ---- end namespace genepi
