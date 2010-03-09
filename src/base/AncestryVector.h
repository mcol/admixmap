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
/// and by a simple integral index (FGIdx) for comparison and probability
/// computation.  We might consider an iterator rather than the index externally.
//
//-----------------------------------------------------------------------------

//template < size_t MAX_FNDR_GAMETES >
class AncestryVector
    {

    public:
	/// Integer index into the elements, i.e. founder-gamete-index:
	typedef Pedigree::FGameteIdx FGIdx;

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

	    PopIdx at_unsafe( FGIdx el ) const
		{
		const WordType byte = bytes[ el >> 1 ];
		return ((el&1) != 0) ? (byte >> 4) : (byte & 0x0F);
		}

	    void setAt_unsafe( FGIdx el, PopIdx val )
		{
		WordType & byte = bytes[ el >> 1 ];

		if ( (el & 1) == 0 )
		    { // Lower-order nibble (even/smaller index)
		    byte &= 0xF0;
		    byte |= val;
		    }
		else
		    { // Higher-order nibble (odd/larger index)
		    byte &= 0x0F;
		    byte |= (val << 4);
		    }
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

	AncestryVector( const Pedigree & _ped, PopIdx _K, const DType & val ) :
		ped ( _ped ) ,
		K   ( _K   )
	    {
	    memcpy( data.bytes, val.bytes, AV_MC_BYTES() );
	    }


	//----------------------------------------------------------------------
	// Access as unsigned long: use with caution.  We befriend
	// HiddenStateSpace only so that it can use these methods, which it does
	// to compute the index into its table of emission-probabilities when
	// passed a specific AV and IV during state-space-generation
	// (HiddenStateSpace::getEProb() called by
	// Pedigree::accumStateInArray().
	//----------------------------------------------------------------------

	friend class HiddenStateSpace;

	/// For use as an index into arrays; in the range of (0,K^F)
	unsigned long to_ulong() const;

	/// @warning set_parms() must be called prior to this method.
	void set_ulong( unsigned long nv );



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
	static Pedigree::FounderIdx founderOf( FGIdx idx )
	    {
	    return idx >> 1;
	    }



	//---------------------------------------------------------------------
	// Index (founder-gamete-index) based read-access to the elements:
	//---------------------------------------------------------------------

	/// Number of elements (founder-gamete-indexed).
	/// Maybe we should cache this somewhere.
	FGIdx size() const { return ped.getNFounderGametes(); }
	PopIdx at_unsafe( FGIdx el ) const		///< at(FGIdx), but not range-checked
	    { return data.at_unsafe(el); }
	PopIdx at( FGIdx el ) const			///< Element access via founder-gamete-index
	    {
	    gp_assert_lt( el, size() );
	    return at_unsafe( el );
	    }
	PopIdx operator[]( FGIdx el ) const { return at(el); } ///< equivalent to at(FGIdx)

	/// See (see <A HREF="#note-1">NOTE *1*</A>) and at(FGIdx)
	PopIdx at( Pedigree::FounderIdx f, Pedigree::GameteType whichOne ) const
	    {
	    return at( ped.founderGameteOfFounder( f, whichOne ) );
	    }


	//---------------------------------------------------------------------
	// Index (founder-gamete-index) based write-access to the elements:
	//---------------------------------------------------------------------

	/// Like @link setAt(FGIdx,PopIdx) setAt() @endlink, but not range-checked.
	void setAt_unsafe( FGIdx el, PopIdx val )
	    { return data.setAt_unsafe( el, val ); }

	/// Set the ancestry at position @a el to @a val.
	void setAt( FGIdx el, PopIdx val )
	    {
	    gp_assert_lt( el , size() );
	    gp_assert_lt( val, K );
	    setAt_unsafe( el, val );
	    }

	/// See (see <A HREF="#note-1">NOTE *1*</A>) and setAt(FGIdx,PopIdx).
	void setAt( Pedigree::FounderIdx f, Pedigree::GameteType whichOne, PopIdx val )
	    {
	    setAt( ped.founderGameteOfFounder(f,whichOne), val );
	    }


	/// Returns true if founder is 2-gamete model, one gamete's ancestry is
	/// from population @a k, the other's is not.  False otherwise,
	/// including if @a f is modeled as a single gamete.  Used by
	/// AdmixPedigree::accumAOScore().
	bool isHetrozygousForPop( const Pedigree::FounderIdx & f, const PopIdx & k ) const
	    {
	    bool rv;

	    const Organism & founder = ped.founderAt( f );

	    if ( ! founder.isHaploid() )
		{
		const bool patAncIsK = at_unsafe(ped.founderGameteOfFounder(f,Pedigree::GT_PATERNAL)) == k;
		const bool matAncIsK = at_unsafe(ped.founderGameteOfFounder(f,Pedigree::GT_MATERNAL)) == k;
		rv = (matAncIsK && (! patAncIsK)) || (patAncIsK && (! matAncIsK));
		}
	    else
		rv = false;

	    return rv;
	    }

	/// The number of gametes of founder @a f with ancestry from population
	/// @a k.  Used by AdmixPedigree::accumAOScore().  Consider passing in
	/// is-haploid flag, eliminating look-up of founder here.
	int nCopiesFromKAtFounder( const Pedigree::FounderIdx & f, const PopIdx & k ) const
	    {
	    int rv;
	    const Organism & founder = ped.founderAt( f );
	    if ( founder.isHaploid() )
		rv = at_unsafe( ped.founderGameteOfFounder(f,Pedigree::GT_SINGLE) ) == k;
	    else
		rv = ( at_unsafe( ped.founderGameteOfFounder(f,Pedigree::GT_PATERNAL) ) == k ) +
		     ( at_unsafe( ped.founderGameteOfFounder(f,Pedigree::GT_MATERNAL) ) == k );
	    return rv;
	    }



	// The nested Iterator class must be forward-defined because it contains
	// a reference to AncestryVector.
	class Iterator;


	/// Hack!  Must be called prior to first call to set_ulong().
	/// Number of populations and maximum number of founder-gametes.
	static void set_parms( PopIdx K, Pedigree::FounderIdx maxF );

    };



//-------------------------------------------------------------------------
/// AncestryVector iterator class -- while there is no "container", it
/// behooves us to have a class that iterates over a "space" defined as
/// all AncestryVectors possible for @a numFndr founders, @a K populations.
//-------------------------------------------------------------------------

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

    /// Useful for debugging: output an AV to an ostream; only the first @a
    /// nValid elements are valid.
    std::ostream & output( std::ostream & os,
			   const AncestryVector & av,
			   AncestryVector::FGIdx nValid );

    /// Useful for debugging: output an AV to an ostream.  Prints one digit for
    /// each founder-gamete, showing the ancestry-population [0..K-1] for that
    /// gamete.  Note that there will only be one gamete for single-gamete
    /// founders; in the case of two-gamete founders, the paternal gamete will
    /// be output first in the pair (this is controlled by
    /// Pedigree::founderGameteOfFounder(), which gives the maternal gamete the
    /// higher index).
    std::ostream & operator<<( std::ostream & os, const AncestryVector & av );

#endif // AV_OSTREAM



} // ---- end namespace genepi



#endif // ! __base_AncestryVector_h
