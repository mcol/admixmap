//=============================================================================
//
// Copyright (C) 2009  David D. Favro  gpl@meta-dynamic.com
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
/// \file Genotype.cc
/// Implementation of the Genotype class.
//=============================================================================

#include "Genotype.h"



namespace genepi { // ----



//-----------------------------------------------------------------------------
// desc()
//-----------------------------------------------------------------------------

estr Genotype::desc() const
    {
    const char LOC_DELIM = ','; // FIXME: gtypeDelim

    estr rv;

    if ( val1 != MISSING_VAL )
	rv = val1;

    if ( isDiploid() )
	{
	rv << LOC_DELIM;

	if ( val2 != MISSING_VAL )
	    rv += val2;
	}

    return rv;
    }



} // ---- end namespace genepi
