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
/// \file HVIterator.h
/// Definition of the HVIterator class.
//=============================================================================



#ifndef __base_HVIterator_h
#define __base_HVIterator_h



#include "Pedigree.h"



namespace genepi { // ----




//-----------------------------------------------------------------------------
//
/// HVIterator class -- iterates across the hidden-variables of a
/// hidden-state-space, in this case the founder-gametes and the meiosis of a
/// pedigree.
///
/// Use like this:
///	const Pedigree & ped;
///	for ( HVIterator it( ped ) ; ! it.isFinished() ; ++it )
///	    {
///	    // Do something...
///	    }
///
/// NOTE *1*: We arbitrarily choose to do the ancestries first, then the meiosis.
///
/// NOTE *2*: We carry attributes of the pedigree in each iterator object; we
///	could instead just keep a reference to the Pedigree object.  This is
///	especially egregious in the case of K, which is fixed for the whole
///	simulation (i.e. global/static).
//
//-----------------------------------------------------------------------------

class HVIterator
    {

    private:
	// See NOTE *2* about all of these:
	const size_t n_ancestries ; // Number of founder-gametes
	const size_t n_meiosis    ; // Number of meiosis
	const size_t K		  ; // Number of populations

	size_t cur_ancestry;
	size_t cur_meiosis;

    public:

	HVIterator( const HiddenStateSpace & hss );
	HVIterator( const Pedigree &	     ped );
	HVIterator( const HVIterator &	     rhs );

	bool isOnMeiosis () const { return (cur_ancestry < n_ancestries); } // NOTE *1*
	bool isOnAncestry() const { return (! isOnMeiosis()); }

	bool isFinished() const { return (cur_ancestry == n_ancestries) && (cur_meiosis == n_meiosis); }
	operator bool() const { return (! isFinished()); }

	void advance();

	// The number of discrete values that the "pointed-to" hidden-variable can take on.
	int nValues() const;

    };



} // ---- end namespace genepi




#endif // ! __base_HVIterator_h
