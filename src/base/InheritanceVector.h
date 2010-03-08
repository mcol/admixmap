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
/// \file InheritanceVector.h
/// Implementation of the InheritanceVector class.
//=============================================================================

#ifndef __base_InheritanceVector_h
#define __base_InheritanceVector_h


#include <bitset>

#include "bclib/exceptions.h"
#include "bclib/cvector.h"


#define IV_OSTREAM		1 /// Should we compile ostream support for IVs?

#define IV_PEDIGREE_FORWARD	1
//#define IV_MAX_BITS		128 /// 64 max non-founders
#define IV_MAX_BITS		32  /// 16 max non-founders



#if IV_OSTREAM
    #include <iosfwd>
#endif



namespace genepi { // ----



/** \addtogroup base
 * @{ */



#if IV_PEDIGREE_FORWARD
    class Pedigree;
#endif

class InheritanceSpace;



//-----------------------------------------------------------------------------
//
/// Inheritance vector for a Pedigree group.  This is effectively a bitmap with
/// one bit for each meiosis; this would typically be two bits for each
/// non-founder member of a pedigree of diploid organisms: one for the father's
/// meiosis and one for the mother's; however, if we are only modeling a single
/// gamete for one of the parents, there will be no meiosis for that parent.
///
/// Internally, each meiosis is represented as a single bit, 1 for the maternal
/// gamete (parent's mother's allele received by child) and 0 for the paternal
/// gamete (parent's father's allele received by child).  Externally, we
/// represent a meiosis via the InheritanceVector::Bits class, which contains
/// both meiosis for a given child; if one of the parents had no meiosis
/// (haploid or single-gamete-model), the corresponding meiosis will have value
/// @a SI_NONE.
///
/// <A name="note-1"></A>
/// <TABLE STYLE="border: groove 3pt aqua;">
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
	/// Segregation Indicator for one meiosis:
	enum SegInd
	    {
	    SI_PATERNAL , /// Child received the parent's father's allele
	    SI_MATERNAL , /// Child received the parent's mother's allele
	    SI_NONE	  /// No meiosis: parent had haploid genotype or we only model one gamete
	    };

	/// Complete segregation indicator (both maternal and paternal) for
	/// diploid organism's inheritance.
	class Bits
	    {
	    friend class InheritanceVector;

	    private:
		SegInd pBit : 8;
		SegInd mBit : 8;

	    protected:
		Bits( unsigned int x ) :
			pBit( SegInd((x & 1) != 0) ) ,
			mBit( SegInd((x & 2) != 0) ) {}

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

	const InheritanceSpace & space;
	std::bitset<IV_MAX_BITS> bits;


    protected:
	// Access as unsigned long: use with caution.  We befriend
	// HiddenStateSpace only so that it can use these methods.
	friend class HiddenStateSpace;
	unsigned long to_ulong  () const { return bits.to_ulong(); }
	void	      set_ulong ( unsigned long nv ) { bits = nv; }

    public:

	const InheritanceSpace & getSpace() const { return space; }

	size_t getNFounders() const;
	size_t getNMembers () const;
	size_t getNNonFndrs() const { return (getNMembers() - getNFounders()); }


	/// Return the number of meiosis (i.e. the number of bits) in this vector.
	size_t getNMeiosis() const;


	InheritanceVector( const Pedigree & p );

	InheritanceVector( const InheritanceVector & rhs ) :
		    space ( rhs.space ) ,
		    bits  ( rhs.bits  )
	    {
	    gp_assert_lt( getNNonFndrs(), MAX_ORGANISMS );
	    }

	InheritanceVector & operator=( const InheritanceVector & rhs )
	    {
	    gp_assert( &space == &rhs.space );
	    bits = rhs.bits;
	    return *this;
	    }


	/// Extract both inheritance bits for a given member of the pedigree
	Bits getMember( SibIdx mIdx ) const;


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
	void setMember( SibIdx mIdx, const Bits & nv );

	// Would be nice to make a member-reference type so we can implement this:
	// BitsReference & operator[]( SibIdx mIdx );


	/// Convenience: did father contribute his father's or his mother's?
	SegInd paternal( SibIdx n ) const { return getMember(n).paternal(); }

	/// Convenience: did mother contribute her father's or her mother's?
	SegInd maternal( SibIdx n ) const { return getMember(n).maternal(); }



	//-----------------------------------------------------------------------
	/// How many meiosis differ between two InheritanceVector's?  This turns
	/// out to be useful for computation of transition-probabilities, so we
	/// would like a convenient and efficient way to compute it.  The two
	/// IVs must come from the same space (i.e. be for the same Pedigree).
	/// At present, this returns the raw XOR of the two IVs when treated as
	/// bitmaps (rather than the number of non-zero bits in that number),
	/// and thus is a highly specialized method that links the internal
	/// representation of InheritanceVector and the internal algorithms of
	/// TransProbCache.  However, for the purposes of the
	/// transition-probability computation as implemented in TransProbCache,
	/// the order of the bits in the bitmap of InheritanceVector is
	/// irrelevant; only the number of non-zero bits in the XOR is used in
	/// the computation.  This could, however, be made cleaner and more
	/// generic by implementing option (C) of the note at
	/// TPC_PRE_CALC_IV_FACTORS in TransProbCache.h and thus returning the
	/// true number of differing meiosis, resulting in a slight
	/// runtime/memory tradeoff.
	///
	/// See also: TransProbCache.h
	//-----------------------------------------------------------------------
	unsigned long operator^( const InheritanceVector & rhs ) const
	    {
	    return to_ulong() ^ rhs.to_ulong();
	    }


	// The nested Iterator class must be forward-defined because it contains
	// a reference to InheritanceVector.
	class Iterator;

    };



//-----------------------------------------------------------------------------
// Iterator
/// Allow iteration across the "space" of all possible inheritance vectors for a
/// given pedigree.  In fact, this space is the same for all pedigrees sharing a
/// common graph; except that it also differs for human males at loci on the X
/// chromosome.  As currently implemented, it is instantiated for a given
/// Pedigree.  No particular ordering of the space is guaranteed.  Note that in
/// the special case of a pedigree with no meiosis (typically a single unrelated
/// individual), we create a "null" inheritance vector and thus a space of size 1.
//-----------------------------------------------------------------------------

class InheritanceVector::Iterator
    {
    private:
	// For efficiency, keep a copy of the inheritance vector that
	// will be used when dereferencing the iterator.
	InheritanceVector pattern;

	unsigned long cur_val;
	unsigned long max_val;

    public:
	Iterator( const Pedigree & _ped );

	bool is_on_last_el() const { return (cur_val == max_val); }

	/// Advances to the next IV in the space.  Returns true if did advanced;
	/// returns false if was _already_ on the last element before the call.
	bool advance() { if ( is_on_last_el() ) return false; pattern.set_ulong( ++cur_val ); return true; }

	const InheritanceVector & getIV	    () const { return pattern; }
	const InheritanceVector & operator* () const { return getIV(); }
	const Iterator *	  operator->() const { return this; }


	/// For efficiency and convenience, we will allow an iterator to be
	/// reset to the the beginning of the space, starting a new complete
	/// iteration, rather than forcing the construction of a new iterator.
	void reset();
    };




//-----------------------------------------------------------------------------
//
// InheritanceSpace
//
/// Consider maintaining two for each Pedigree, one for the X chromosome and one
/// for the non-X chromosomes, and eliminating the IsXChrom parameters.
//
//-----------------------------------------------------------------------------

class InheritanceSpace
    {

    public:

	enum IsXChrom
	    {
	    IS_X_CHROM	= 0 ,
	    NOT_X_CHROM = 1
	    };
	typedef size_t NonFounderIdx;

	struct WhichMeiosis
	    {
	    bool haveMaternalMeiosis;
	    bool havePaternalMeiosis;
	    };


    private:

	friend class InheritanceVector;

	struct Entry
	    {
	    static const size_t HAPLOID_PARENT = size_t(-1);

	    WhichMeiosis wm;
	    size_t	 p_bit;
	    size_t	 m_bit;
	    };

	cvector<Entry>	non_x_vals;
	cvector<Entry>	x_vals;
	size_t		nFounders;
	size_t		nMembers;
	size_t		nMeiosis_X;
	size_t		nMeiosis_nonX;

	// Used by friend InheritanceVector.  Remove the default value for
	/// @a isX when X chromosomes are fully implemented for pedigrees.
	const Entry & getEntry( InheritanceVector::SibIdx mIdx, IsXChrom isX = NOT_X_CHROM ) const
	    {
	    gp_assert_ge( mIdx, nFounders );
	    gp_assert_lt( mIdx, nMembers  );
	    const size_t ent_idx = mIdx - nFounders;
	    return (isX == IS_X_CHROM) ? x_vals[ent_idx] : non_x_vals[ent_idx];
	    }

    public:

	InheritanceSpace( const Pedigree & _ped );

	size_t getNFounders() const { return nFounders; }
	size_t getNMembers () const { return nMembers ; }

	/// The number of meiosis in the pedigree.  Remove the default value for
	/// @a isX when X chromosomes are fully implemented for pedigrees.
	size_t getNMeiosis( IsXChrom isX = NOT_X_CHROM ) const
	    {
	    return (isX == IS_X_CHROM) ? nMeiosis_X : nMeiosis_nonX;
	    }

	/// Given a member of a pedigree, report which meiosis contributed to its genotype.
	WhichMeiosis getWhichMeiosis( NonFounderIdx nonFounderIdx, IsXChrom x_chrom ) const;

    };



//-----------------------------------------------------------------------------
// Inline InheritanceVector methods that use InheritanceSpace's methods:
//-----------------------------------------------------------------------------


inline size_t InheritanceVector::getNFounders() const
    {
    return space.getNFounders();
    }


inline size_t InheritanceVector::getNMembers() const
    {
    return space.getNMembers();
    }


inline size_t InheritanceVector::getNMeiosis() const
    {
    return getSpace().getNMeiosis();
    }



inline InheritanceVector::Bits InheritanceVector::getMember( SibIdx mIdx ) const
    {

    gp_assert_lt( mIdx, getNMembers() );
    gp_assert_ge( mIdx, getNFounders() );

    const InheritanceSpace::Entry & entry = space.getEntry( mIdx );

    SegInd patSI;
    if ( entry.p_bit == InheritanceSpace::Entry::HAPLOID_PARENT )
	patSI = SI_NONE;
    else
	patSI = bits[ entry.p_bit ] ? SI_MATERNAL : SI_PATERNAL;

    SegInd matSI;
    if ( entry.m_bit == InheritanceSpace::Entry::HAPLOID_PARENT )
	matSI = SI_NONE;
    else
	matSI = bits[ entry.m_bit ] ? SI_MATERNAL : SI_PATERNAL;

    return Bits( patSI, matSI );

    }



inline void InheritanceVector::setMember( SibIdx mIdx, const Bits & nv )
    {
    gp_assert_lt( mIdx, getNMembers()  );
    gp_assert_ge( mIdx, getNFounders() );

    const InheritanceSpace::Entry & entry = space.getEntry( mIdx );

    if ( nv.paternal() == SI_NONE )
	{
	gp_assert( entry.p_bit == InheritanceSpace::Entry::HAPLOID_PARENT );
	}
    else
	{
	gp_assert( entry.p_bit != InheritanceSpace::Entry::HAPLOID_PARENT );
	bits[ entry.p_bit ] = (nv.paternal() == SI_PATERNAL) ? 0 : 1;
	}

    if ( nv.maternal() == SI_NONE )
	{
	gp_assert( entry.m_bit == InheritanceSpace::Entry::HAPLOID_PARENT );
	}
    else
	{
	gp_assert( entry.m_bit != InheritanceSpace::Entry::HAPLOID_PARENT );
	bits[ entry.m_bit ] = (nv.maternal() == SI_PATERNAL) ? 0 : 1;
	}
    }



//-----------------------------------------------------------------------------
// Other inline methods:
//-----------------------------------------------------------------------------


#if IV_OSTREAM
    enum IVOutputStyle { IV_BINARY , IV_ALPHA };
    void setIVOutputStyle( IVOutputStyle );
    /// Useful for debugging: output an IV to an ostream.  Each segregation
    /// indicator is output as a '0' for paternal and a '1' for maternal.  For
    /// diploid inheritance, the meiosis on the father's side is output first
    /// and the mother's meiosis is second.
    std::ostream & operator<<( std::ostream & os, const InheritanceVector & iv );
#endif // IV_OSTREAM




} // ---- end namespace genepi



/** @} */



#endif // ! __base_InheritanceVector_h
