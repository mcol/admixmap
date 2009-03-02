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
/// \file InheritanceVector.h
/// Implementation of the InheritanceVector class.
//=============================================================================

#ifndef __base_InheritanceVector_h
#define __base_InheritanceVector_h



#include <bitset>

#include "exceptions.h"


#define IV_KEEP_PED_REF		0
#define IV_PEDIGREE_FORWARD	1
#define IV_MAX_BITS		128 // 64 max non-founders



namespace genepi { // ----



/** \addtogroup base
 * @{ */



#if IV_PEDIGREE_FORWARD
    class Pedigree;
#endif



//-----------------------------------------------------------------------------
///
/// Inheritance vector for a Pedigree group.  This may require generalization
/// for polyploidal organisms.
///
/// This is effectively a bitmap with two bits for each non-founder member of
/// the pedigree.
///
/// <A name="note-1"></A>
/// <TABLE STYLE="border: groove 3pt aqua;">
///
///  <TR>
///	<TD><B>NOTE *1*</B></TD>
///	<TD>
///	Legacy, nothing here.
///	</TD>
///  </TR>
///
///  <TR>
///	<TD><B>NOTE *2*</B></TD>
///	<TD>
///	We could keep a reference to the Pedigree object in this class, but due
///	to a desire to implement inline operations in this header (thus
///	precluding forward-references), and a dizzying circular dependency
///	between this class and Pedigree (as they are currently implemented,
///	because the hidden-state enumeration code is in Pedigree itself), it now
///	just caches the relevant parameters from Pedigree (number of members,
///	number of founders, etc.), which are immutable once the Pedigree object
///	is constructed.
///	</TD>
///  </TR>
///
/// </TABLE>
///
//-----------------------------------------------------------------------------

class InheritanceVector
    {
    public:
	// Segregation Indicator for one meosis:
	enum SegInd
	    {
	    SI_PATERNAL ,
	    SI_MATERNAL
	    };

	// Complete segregation indicator (both maternal and paternal) for
	// diploid organism's inheritance:
	class Bits
	    {
	    friend class InheritanceVector;

	    private:
		SegInd pBit : 8;
		SegInd mBit : 8;

	    protected:
		Bits( unsigned int x ) : pBit( (x & 1) != 0 ), mBit( (x & 2) != 0 ) {}

	    public:
		Bits( SegInd p, SegInd m ) : pBit(p), mBit(m) {}

		SegInd paternal() const { return pBit; } ///< Did father contribute his father's or his mother's?
		SegInd maternal() const { return mBit; } ///< Did mother contribute her father's or her mother's?
	    };

	/// This is the offset from the first non-founder; it is safer to access
	/// by Individual and let the InheritanceVector calculate the index.
	#if IV_PEDIGREE_FORWARD
	    typedef size_t SibIdx;
	#else
	    typedef Pedigree::SibIdx SibIdx;
	#endif

	static const size_t MAX_ORGANISMS = (IV_MAX_BITS >> 1);

    private:
	#if IV_KEEP_PED_REF
	    const Pedigree & pedigree ;
	#else
	    unsigned short nFounders;
	    unsigned short nMembers ;
	#endif

	std::bitset<IV_MAX_BITS> bits;

    public:

	#if IV_KEEP_PED_REF
	    const Pedigree & getPedigree() const { return pedigree; }
	    size_t getNFounders() const { return pedigree.getNFounders(); }
	    size_t getNMembers () const { return pedigree.getNMembers (); }
	#else
	    size_t getNFounders() const { return nFounders; }
	    size_t getNMembers () const { return nMembers ; }
	#endif

	size_t getNSibs() const { return (getNMembers() - getNFounders()); }

	#if IV_KEEP_PED_REF
	    InheritanceVector( const Pedigree & p ) :
		pedigree( p ) {}

	    InheritanceVector( const InheritanceVector & rhs ) :
		pedigree( rhs.pedigree ) ,
		bits	( rhs.bits     ) {}

	    InheritanceVector & operator=( const InheritanceVector & rhs )
		{
		gp_assert( &rhs.pedigree == &pedigree );
		bits = rhs.bits;
		return *this;
		}
	#else
	    InheritanceVector( const Pedigree & p );

	    InheritanceVector( const InheritanceVector & rhs ) :
		    nFounders( rhs.nFounders ) ,
		    nMembers ( rhs.nMembers  ) ,
		    bits     ( rhs.bits	     )
		{
		gp_assert_lt( getNSibs(), MAX_ORGANISMS );
		}

	    InheritanceVector & operator=( const InheritanceVector & rhs )
		{
		nFounders = rhs.nFounders ;
		nMembers  = rhs.nMembers  ;
		bits	  = rhs.bits	  ;
		return *this;
		}
	#endif


	/// Extract both inheritance bits for a given member of the pedigree
	Bits getMember( SibIdx mIdx ) const
	    {
	    gp_assert_lt( mIdx, getNMembers () );
	    gp_assert_ge( mIdx, getNFounders() );
	    const size_t shift = ((mIdx - getNFounders()) << 1);
	    return Bits( bits[shift    ] ? SI_MATERNAL : SI_PATERNAL ,
			 bits[shift + 1] ? SI_MATERNAL : SI_PATERNAL );
	    }


	#if ! IV_PEDIGREE_FORWARD
	  Bits getMember( const Pedigree::Member & organism ) const
	    {
	    #if 0 // See *TD1*
		gp_assert( pedigree.isAMember( organism );
	    #endif
	    gp_assert( ! organism.isFounder() );
	    return getMember( organism.getPIdx() );
	    }
	#endif


	/// Set both inheritance bits for a given member of the pedigree:
	void setMember( SibIdx mIdx, const Bits & nv )
	    {
	    gp_assert_lt( mIdx, getNMembers () );
	    gp_assert_ge( mIdx, getNFounders() );
	    const size_t bit_shift = ((mIdx - getNFounders()) << 1);
	    bits[ bit_shift	] = nv.paternal() == SI_PATERNAL ? 0 : 1;
	    bits[ bit_shift + 1 ] = nv.maternal() == SI_PATERNAL ? 0 : 1;
	    }

	// Would be nice to make a member-reference type so we can implement this:
	// BitsReference & operator[]( SibIdx mIdx );


	/// Convenience: did father contribute his father's or his mother's?
	SegInd paternal( SibIdx n ) const { return getMember(n).paternal(); }
	#if ! IV_PEDIGREE_FORWARD
	  SegInd paternal( const Pedigree::Member & organism ) const
			{ return getMember(organism).paternal(); }
	#endif


	/// Convenience: did mother contribute her father's or her mother's?
	SegInd maternal( SibIdx n ) const { return getMember(n).maternal(); }
	#if ! IV_PEDIGREE_FORWARD
	  SegInd maternal( const Pedigree::Member & organism ) const
			{ return getMember(organism).maternal(); }
	#endif
    };



#if (! IV_KEEP_PED_REF) && (! IV_PEDIGREE_FORWARD)

    inline InheritanceVector::InheritanceVector( const Pedigree & p ) :
	    nFounders( p.getNFounders() ) ,
	    nMembers ( p.getNMembers () )
	{
	gp_assert_le( p.getNSibs(), MAX_ORGANISMS );
	}

#endif



} // ---- end namespace genepi



/** @} */



#endif // ! __base_InheritanceVector_h
