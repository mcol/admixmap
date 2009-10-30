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
/// \file InheritanceVector.cc
/// Implementation of the InheritanceVector class
//=============================================================================

#include "InheritanceVector.h"

#include "Pedigree.h"



#if IV_OSTREAM
    #include <iostream>
#endif



namespace genepi { // ----



const size_t InheritanceVector::MAX_ORGANISMS;



#if (! IV_KEEP_PED_REF) && IV_PEDIGREE_FORWARD

    InheritanceVector::InheritanceVector( const Pedigree & p ) :
	    nFounders( p.getNFounders() ) ,
	    nMembers ( p.getNMembers () )
	{
	gp_assert_le( p.getNNonFndrs(), MAX_ORGANISMS );
	}

#endif



//-----------------------------------------------------------------------------
// Iterator
//-----------------------------------------------------------------------------

void InheritanceVector::Iterator::reset()
    {
    cur_val = 0;
    pattern.set_ulong( cur_val );
    }


InheritanceVector::Iterator::Iterator( const Pedigree & _ped ) :
	pattern( _ped ),
	max_val( (1UL << pattern.n_meiosis()) - 1 )
    {
    reset();
    }



//-----------------------------------------------------------------------------
// Print an InheritanceVector to an ostream
//-----------------------------------------------------------------------------

#if IV_OSTREAM

    static IVOutputStyle style = IV_BINARY;

    void setIVOutputStyle( IVOutputStyle nv )
	{
	style = nv;
	}

    inline static char si2char( const InheritanceVector::SegInd & si )
	{
	return (si == InheritanceVector::SI_PATERNAL) ? 'p' : 'm';
	}

    inline static char si2digit( const InheritanceVector::SegInd & si )
	{
	return (si == InheritanceVector::SI_PATERNAL) ? '0' : '1';
	}

    inline static std::ostream & operator<<( std::ostream & os, const InheritanceVector::Bits & b )
	{
	if ( style == IV_BINARY )
	    return os << si2digit(b.paternal()) << ',' << si2digit(b.maternal());
	else
	    return os << char(toupper(si2char(b.paternal()))) << si2char(b.maternal());
	}

    std::ostream & operator<<( std::ostream & os, const InheritanceVector & iv )
	{
	os << "IV(";

	if ( iv.n_meiosis() == 0 )
	    os << "-no-child-)";
	else
	    {
	    const size_t limit = iv.getNMembers() - 1;
	    for ( size_t sib = iv.getNFounders() ; sib < limit ; ++sib )
		os << iv.getMember(sib) << ';';
	    os << iv.getMember(limit) << ')';
	    }

	return os;
	}

#endif // IV_OSTREAM



} // ---- end namespace genepi
