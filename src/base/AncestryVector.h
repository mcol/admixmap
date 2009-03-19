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
#include "exceptions.h"
#include <vector>



#define AV_OSTREAM		1 /// Should we compile ostream support for AVs?

#if AV_OSTREAM
    #include <iosfwd>
#endif



#define AV_STORAGE unsigned char // PopIdx



namespace genepi { // ----



//-----------------------------------------------------------------------------
//
/// Since we may be storing many of these vectors in memory, we'll optimize the
/// internal storage of the elements as "unsigned char", leading to a maximum of
/// 256 populations; but externally we present them as PopIdx (typically
/// size_t).  This leads to private inheritance from std::vector.
///
/// NOTE *2*: we could consider exclusively indexing based on the founder-index
///	and which gamete; i.e. using (Pedigree::FounderIdx,MaternalPaternalType)
///	to reference the elements.
//
//-----------------------------------------------------------------------------

class AncestryVector : private std::vector<AV_STORAGE>
    {
    private:
	const Pedigree & ped;
	const PopIdx	 K;

	typedef std::vector<AV_STORAGE> SUPER;


    protected:
	/// See NOTE *2*:
	static size_t idxOf( Pedigree::FounderIdx f, bool isPaternalGamete )
	    {
	    size_t idx = f << 1;
	    if ( isPaternalGamete )
		++idx;
	    return idx;
	    }

    public:
	AncestryVector( const Pedigree & _ped, PopIdx _K ) :
	    SUPER ( (_ped.getNFounders() << 1)  ) ,
	    ped	  ( _ped			) ,
	    K	  ( _K				) {}

	AncestryVector( const AncestryVector & rhs ) :
	    SUPER ( rhs		) ,
	    ped	  ( rhs.ped	) ,
	    K	  ( rhs.K	) {}


	/// For use as an index into arrays; in the range of (0,K^F)
	unsigned long to_ulong() const;


	size_t size() const { return SUPER::size(); } ///< Publicized vector method
	PopIdx at( size_t el ) const		      ///< Element access (publicized vector method)
	    {
	    gp_assert_lt( el, size() );
	    return SUPER::operator[](el);
	    }
	PopIdx operator[]( size_t el ) const { return at(el); } ///< equivalent to at(size_t)

	/// See NOTE *2*.  @see { at(size_t) }
	PopIdx at( Pedigree::FounderIdx f, bool isPaternalGamete ) const
	    {
	    return at( idxOf(f,isPaternalGamete) );
	    }


	/// Set the ancestry at position @a el to @a val
	void setAt( size_t el, PopIdx val )
	    {
	    gp_assert_lt( el , size() );
	    gp_assert_lt( val, K      );
	    SUPER::operator[](el) = val;
	    }

	/// See NOTE *2*
	void setAt( Pedigree::FounderIdx f, bool isPaternalGamete, PopIdx val )
	    {
	    setAt( idxOf(f,isPaternalGamete), val );
	    }
    };



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
