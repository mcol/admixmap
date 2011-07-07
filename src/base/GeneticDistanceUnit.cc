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
/// \file GeneticDistanceUnit.cc
/// Implementation of the GeneticDistanceUnit methods.
//=============================================================================

#pragma implementation
#include "GeneticDistanceUnit.h"

#include <stdexcept>
#include <cstring>  // strcasecmp()
#include "config.h" // AGGRESSIVE_RANGE_CHECK



static const char * const descs[] = { "bp", "kb", "Mb", "cM", "M" };

static const int N_DESCS = sizeof(descs) / sizeof(*descs);



//-----------------------------------------------------------------------------
// gduAsString()
//-----------------------------------------------------------------------------

const char * gduAsString( GeneticDistanceUnit u )
    {
    #if AGGRESSIVE_RANGE_CHECK
	if ( N_DESCS != N_GDUS )
	    throw std::runtime_error( "internal GDU desync" );

	if ( size_t(u) >= N_GDUS )
	    throw std::runtime_error( "invalid GDU" );
    #endif

    return descs[ u ];
    }



//-----------------------------------------------------------------------------
// gduFromString()
//-----------------------------------------------------------------------------

GeneticDistanceUnit gduFromString( const char * desc )
    {
    #if AGGRESSIVE_RANGE_CHECK
	if ( N_DESCS != N_GDUS )
	    throw std::runtime_error( "internal GDU desync" );
    #endif

    for ( size_t idx = N_DESCS ; idx-- != 0 ; )
	if ( strcasecmp( desc, descs[idx] ) == 0 )
	    return static_cast<GeneticDistanceUnit>( idx );

    throw std::runtime_error( std::string("invalid genetic distance unit: \"") + desc + '"' );
    }
