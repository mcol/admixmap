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
/// \file GenotypeIterator.cc
/// Implementation of the GenotypeIterator class
//=============================================================================

#include "GenotypeIterator.h"

#include <stdexcept>
#include "estr.h"



using namespace genepi;



#if USE_GENOTYPE_PARSER
    unsigned short AdmixGenotypeIterator::operator() ( unsigned int j , unsigned int g ) const
	{
	const Genotype & gtype = v.at( j );

	if ( g == 0 )
	    return gtype.getVal1( 0 ); // 0 is the secret undocumented missing value
	else if ( g == 1 )
	    return gtype.getVal2( 0 ); // 0 is the secret undocumented missing value
	else
	    throw std::invalid_argument( estr("Invalid allele-index: ") + g );
	}
#endif
