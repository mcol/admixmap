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
/// \file AncestryVector.h
/// Definition of the AncestryVector class.
//=============================================================================

#ifndef __base_AncestryVector_h
#define __base_AncestryVector_h



#include "Pedigree.h"
#include "Organism.h"	// PopIdx
#include "bclib/exceptions.h"

#include <cstring> // memcpy()


#define AV_OSTREAM		1 ///< Should we compile ostream support for AVs?

#if AV_OSTREAM
    #include <iosfwd>
#endif



#define AV_MAX_FOUNDER_GAMETES	16 ///< Should be even
#define AV_MAX_K		16



namespace genepi { // ----



//-----------------------------------------------------------------------------
//
/// Stores an ancestry vector for a pedigree: one element for each founder
/// gamete; each element is a population (PopIdx), a number from 0 to K-1.
///
/// Since we may be storing many of these vectors in memory, we'll optimize the
/// internal storage of the elements as "unsigned char", leading to a maximum of
/// 256 populations; but externally we present them as PopIdx (typically
/// size_t).  This leads to private inheritance from std::vector.
///
/// <A name="note-1"></A>
/// NOTE *1*: We need to index both based on the founder-index and which gamete
/// (i.e. using {Pedigree::FounderIdx,MaternalPaternalType}) for inheritance,
/// and by a simple integral index (IdxType) for comparison and probability
/// computation.  We might consider an iterator rather than the index externally.
//
//-----------------------------------------------------------------------------

//template < size_t MAX_FNDR_GAMETES >
class AncestryVector
    {
    public:
	/// Integer index into the elements, i.e. founder-gamete-index:
	typedef size_t IdxType;

    protected:

	typedef unsigned char WordType;
	static const size_t MAX_BYTES = (AV_MAX_FOUNDER_GAMETES + 1) >> 1;

    public:
	/// We publicize this internal component only so that we can use it from
	/// the IdxCache class that is internal to the implementation, without
	/// making IdxCache a friend/member of AncestryVector.  DType is not
	/// usable from outside AncestryVector.
	struct DType
	    {

	    WordType bytes[ MAX_BYTES ];

	    PopIdx at_unsafe( IdxType el ) const
		{
		const WordType byte = bytes[ el >> 1 ];
		return ((el&1) != 0) ? (byte >> 4) : (byte & 0x0F);
		}

	    void setAt_unsafe( IdxType el, PopIdx val )
		{
		WordType & byte = bytes[ el >> 1 ];

		if ( (el & 1) == 0 )
		    { // Lower-order nibble (even index, maternal)
		    byte &= 0xF0;
		    byte |= val;
		    }
		else
		    { // Higher-order nibble (odd index, paternal)
		    byte &= 0x0F;
		    byte |= (val << 4);
		    }
		}

	    PopIdx getBoth( const Pedigree::FounderIdx & f, PopIdx & maternalAncestry ) const
		{
		const WordType byte = bytes[ f ];
		maternalAncestry = byte & 0x0F;
		return byte >> 4;
		}

	    int nCopiesFromKAtFounder( const Pedigree::FounderIdx & f, const PopIdx & k ) const
		{
		const WordType byte = bytes[ f ];
		const PopIdx matAnc = byte & 0x0F;
		const PopIdx patAnc = byte >> 4;
		return (matAnc == k) + (patAnc == k);
		}

	    bool isHetrozygousForPop( const Pedigree::FounderIdx & f, const PopIdx & k ) const
		{
		const WordType byte = bytes[ f ];
		const PopIdx matAnc = byte & 0x0F;
		const PopIdx patAnc = byte >> 4;
		return	((matAnc == k) && (patAnc != k)) ||
			((patAnc == k) && (matAnc != k));
		}

	    };


    private:
	const Pedigree & ped;
	const PopIdx	 K;
	DType		 data;

	/// Maybe we should cache this somewhere.
	size_t nBytes() const { return ((size()+1) >> 1); }
	#if 1
	    #define AV_MC_BYTES() MAX_BYTES
	#else
	    #define AV_MC_BYTES() nBytes()
	#endif

    protected:
	/// Converts between (founder-index,which-gamete) and
	/// founder-gamete-index (see <A HREF="#note-1">NOTE *1*</A>).
	static IdxType idxOf( Pedigree::FounderIdx f, bool isPaternalGamete )
	    {
	    IdxType idx = f << 1;
	    if ( isPaternalGamete )
		++idx;
	    return idx;
	    }

	AncestryVector( const Pedigree & _ped, PopIdx _K, const DType & val ) :
		ped   ( _ped ) ,
		K     ( _K   )
	    {
	    memcpy( data.bytes, val.bytes, AV_MC_BYTES() );
	    }

    public:
	AncestryVector( const Pedigree & _ped, PopIdx _K ) :
		ped ( _ped ) ,
		K   ( _K   ) {}

	AncestryVector( const AncestryVector & rhs ) :
		ped ( rhs.ped ) ,
		K   ( rhs.K   )
	    {
	    memcpy( data.bytes, rhs.data.bytes, AV_MC_BYTES() );
	    }



	/// Converts between founder-gamete-index and founder-index.
	static Pedigree::FounderIdx founderOf( IdxType idx )
	    {
	    return idx >> 1;
	    }


	/// For use as an index into arrays; in the range of (0,K^F)
	unsigned long to_ulong() const;

	/// @warning set_parms() must be called prior to this method.
	void set_ulong( unsigned long nv );



	//---------------------------------------------------------------------
	// Index (founder-gamete-index) based read-access to the elements:
	//---------------------------------------------------------------------

	/// Number of elements (founder-gamete-indexed).
	/// Maybe we should cache this somewhere.
	IdxType size() const { return ped.getNFounderGametes(); }
	PopIdx at_unsafe( IdxType el ) const		///< at(IdxType), but not range-checked
	    { return data.at_unsafe(el); }
	PopIdx at( IdxType el ) const			///< Element access via founder-gamete-index
	    {
	    gp_assert_lt( el, size() );
	    return at_unsafe( el );
	    }
	PopIdx operator[]( IdxType el ) const { return at(el); } ///< equivalent to at(IdxType)

	/// See (see <A HREF="#note-1">NOTE *1*</A>) and at(IdxType)
	PopIdx at( Pedigree::FounderIdx f, bool isPaternalGamete ) const
	    {
	    return at( idxOf(f,isPaternalGamete) );
	    }


	//---------------------------------------------------------------------
	// Index (founder-gamete-index) based write-access to the elements:
	//---------------------------------------------------------------------

	/// Like @link setAt(IdxType,PopIdx) setAt() @endlink, but not range-checked.
	void setAt_unsafe( IdxType el, PopIdx val )
	    { return data.setAt_unsafe( el, val ); }

	/// Set the ancestry at position @a el to @a val.
	void setAt( IdxType el, PopIdx val )
	    {
	    gp_assert_lt( el , size() );
	    gp_assert_lt( val, K );
	    setAt_unsafe( el, val );
	    }

	/// See (see <A HREF="#note-1">NOTE *1*</A>) and setAt(IdxType,PopIdx).
	void setAt( Pedigree::FounderIdx f, bool isPaternalGamete, PopIdx val )
	    {
	    setAt( idxOf(f,isPaternalGamete), val );
	    }


	/// Retrieve two ancestry values at once (for the maternal and paternal
	/// gametes of one individual).  Pass in the founder index in @a f;
	/// returns the paternal gamete's ancestry, and puts the maternal
	/// gamete's ancestry into @a maternalAncestry.  Implementation note:
	/// this depends on the internal implementation of DType in a slightly
	/// inelegant way; but getBoth() can improve performance in
	/// certain situations, and it does not expose the implementation
	/// details outside of the AncestryVector class.
	PopIdx getBoth( const Pedigree::FounderIdx & f, PopIdx & maternalAncestry ) const
	    {
	    return data.getBoth( f, maternalAncestry );
	    }

	bool isHetrozygousForPop( const Pedigree::FounderIdx & f, const PopIdx & k ) const
	    {
	    return data.isHetrozygousForPop( f, k );
	    }

	int nCopiesFromKAtFounder( const Pedigree::FounderIdx & f, const PopIdx & k ) const
	    {
	    return data.nCopiesFromKAtFounder( f, k );
	    }



	//-------------------------------------------------------------------------
	/// AncestryVector iterator class -- while there is no "container", it
	/// behooves us to have a class that iterates over a "space" defined as
	/// all AncestryVectors possible for @a numFndr founders, @a K populations.
	//-------------------------------------------------------------------------

	class Iterator;


	/// Hack!  Must be called prior to first call to set_ulong().
	/// Number of populations and maximum number of founder-gametes.
	static void set_parms( PopIdx K, Pedigree::FounderIdx maxF );

    };


// C++ bizarrely forbids defining this nested class within the surrounding class.
class AncestryVector::Iterator
    {
    private:
	unsigned long  idx;
	unsigned long  idxLimit;
	AncestryVector av;

    public:

	Iterator( const Pedigree & _ped );

	/// Returns true if did _not_ advance _past_ the end of the
	/// space.  If returns false, the iterator is _no_longer_ valid.
	bool advance();
	bool isValid() const { return (idx < idxLimit); }

	unsigned long to_ulong() const { gp_assert(isValid()); return idx; }

	const AncestryVector & getAV() const
	    {
	    gp_assert( isValid() );
	    return av;
	    }
    };



//typedef AncestryVectorT<AV_MAX_FOUNDER_GAMETES> AncestryVector;



#if AV_OSTREAM
    /// Alpha output is not currently supported since the vector keeps no
    /// reference to the population-names vector.
    enum AVOutputStyle { AV_NUMERIC , AV_ALPHA };
    void setAVOutputStyle( AVOutputStyle );
    /// Useful for debugging: output an AV to an ostream
    std::ostream & operator<<( std::ostream & os, const AncestryVector & av );
#endif // AV_OSTREAM



} // ---- end namespace genepi



#endif // ! __base_AncestryVector_h
