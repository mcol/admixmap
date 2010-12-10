//=============================================================================
//
// Copyright (C) 2009  David D. Favro  amm-header@meta-dynamic.com
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
/// \file NRand.cc
/// Implementation of the genepi::NRand class.
//=============================================================================


#include "bclib/NRand.h"


#include <stdexcept>



namespace genepi { // ----
//namespace bclib { // ----



NRand::NRand( const gsl_rng_type * type, unsigned long init_seed ) :
	gen( gsl_rng_alloc( type ) )
    {
    if ( gen == 0 )
	throw std::runtime_error( "unable to allocate GSL PRNG (presumably insufficient memory)" );

    if ( init_seed != gsl_rng_default_seed )
	gsl_rng_set( gen, init_seed );
    }



NRand::NRand( unsigned long init_seed ) :
	gen( gsl_rng_alloc( gsl_rng_taus ) )
    {
    if ( gen == 0 )
	throw std::runtime_error( "unable to allocate GSL PRNG (presumably insufficient memory)" );

    if ( init_seed != gsl_rng_default_seed )
	gsl_rng_set( gen, init_seed );
    }



NRand::~NRand()
    {
    gsl_rng_free( gen );
    }



void NRand::seed( unsigned long new_seed )
    {
    gsl_rng_set( gen, new_seed );
    }



//} // ---- end namespace bclib
} // ---- end namespace genepi
