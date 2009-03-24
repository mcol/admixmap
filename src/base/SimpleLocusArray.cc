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
/// \file SimpleLocusArray.cc
//=============================================================================

#include "SimpleLocusArray.h"

#include <stdexcept>



namespace genepi { // ----



//-----------------------------------------------------------------------------
// throwRange() [protected]
//-----------------------------------------------------------------------------

void SimpleLocusArray::throwRange( size_t idx ) const
    {
    throw std::runtime_error( estr("SimpleLocus idx (") + idx +
		") out-of-range (" + size() + ')' );
    }



//-----------------------------------------------------------------------------
// Constructor:
//-----------------------------------------------------------------------------

SimpleLocusArray::SimpleLocusArray() :
	gdu	     ( Morgans	) ,
	nComposite   ( 0	) ,
	nChromosomes ( 0	)
    {
    }



//-----------------------------------------------------------------------------
// Destructor:
//-----------------------------------------------------------------------------

SimpleLocusArray::~SimpleLocusArray()
    {
    }



//-----------------------------------------------------------------------------
// Locate locus-index based on name.
//
/// For small sets of loci, this is implemented as a linear search; for larger
/// sets we should implement a more efficient algorithm such as a hash-table or
/// sorted-index.
//-----------------------------------------------------------------------------

SLocIdxType SimpleLocusArray::findIndexOf( const std::string & locus ) const
    {
    for ( SLocIdxType idx = size() ; idx-- != 0 ; )
	if ( atUnsafe(idx).getName() == locus )
	    return idx;

    throw std::runtime_error( estr("Unknown locus \"") + locus + '"' );
    }



} // ---- end namespace genepi
