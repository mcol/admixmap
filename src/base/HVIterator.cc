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
/// Implementation of the HVIterator class.
//=============================================================================


#include "HVIterator.h"



namespace genepi { // ----




#if 0
    HVIterator::HVIterator( const HiddenStateSpace & hss ) :
	    n_ancestries( hss.getPed().getNFounderGametes() ) ,
	    n_meiosis	( hss.getPed().getNMeiosis()	    ) ,
	    K		( hss.getPed().getK()		    ) ,
	    cur_ancestry( 0 ) ,
	    cur_meiosis ( 0 )
	{
	}
#endif



HVIterator::HVIterator( const Pedigree & ped ) :
	n_ancestries( ped.getNFounderGametes()	) ,
	n_meiosis   ( ped.getNMeiosis()		) ,
	K	    ( ped.getK()		) ,
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



void HVIterator::advance()
    {
    gp_assert( ! isFinished() ); // Don't allow advance after finished

    // NOTE *1*:
    if ( isOnMeiosis() || (++cur_ancestry == n_ancestries) )
	++cur_meiosis;
    }



int HVIterator::nValues() const
    {
    return (isOnAncestry() ? K : 2);
    }



} // ---- end namespace genepi
