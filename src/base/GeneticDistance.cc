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
/// \file GeneticDistance.cc
/// Implementation of the GeneticDistance class.
//=============================================================================

#include "GeneticDistance.h"

#include <stdexcept>



namespace genepi { // ----



//-----------------------------------------------------------------------------
// set()
//-----------------------------------------------------------------------------

void GeneticDistance::set( GeneticDistanceUnit unit, double val )
    {
    switch ( unit )
	{
	case basepairs:
	    distance = val / 1e6;
	    isBP = true;
	    break;

	case kilobases:
	    distance = val / 1e3;
	    isBP = true;
	    break;

	case megabases:
	    distance = val;
	    isBP = true;
	    break;

	case centimorgans:
	    distance = val / 100.0;
	    isBP = false;
	    break;

	case Morgans:
	    distance = val;
	    isBP = false;
	    break;

	default:
	    throw std::runtime_error( "Invalid GDU" );
	}
    }



} // ---- end namespace genepi
