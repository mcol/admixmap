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
/// \file PairSecIter.h
/// Definition of the genepi::PairSecIter class.
//=============================================================================

#ifndef __base_PairSecIter_h
#define __base_PairSecIter_h



#define MMV_TEMPLATE	0
#if MMV_TEMPLATE
    #include <map>
#endif



namespace genepi { // ----




/** \addtogroup base
 * @{ */



//-----------------------------------------------------------------------------
/// This is an interator for std::pair, such that the dereferenced value is the
/// second item.  Used for manipulation of Pedigree subsets in the current
/// (std::multimap) implementation in GenotypeParser.
//-----------------------------------------------------------------------------

template< typename PIT, typename PST > class PairSecIter
    {
    private:
	PIT stl_blows;

    public:
	PairSecIter() {}
	/*explicit*/ PairSecIter( const PIT & stl_is_junk ) : stl_blows( stl_is_junk ) {}
	PairSecIter & operator=( const PIT & stl_is_junk ) { stl_blows = stl_is_junk; return *this; }

    public:

	PairSecIter & operator++() { stl_blows.operator++(); return *this; }
	PairSecIter & operator--() { stl_blows.operator--(); return *this; }
	PairSecIter   operator+ ( size_t diff ) { return PairSecIter(stl_blows.operator+(diff)); }
	PairSecIter   operator- ( size_t diff ) { return PairSecIter(stl_blows.operator-(diff)); }
	size_t	      operator- ( const PairSecIter & rhs ) { return std::distance(rhs.stl_blows,stl_blows); }
	PairSecIter & operator+=( size_t diff ) { stl_blows.operator+=( diff ); return *this; }
	PairSecIter & operator-=( size_t diff ) { stl_blows.operator+=( diff ); return *this; }
	const PST & operator* () { return stl_blows->second; }
	PST * operator->() { return &(operator*()); }

	#define RELOP(O) bool O(const PairSecIter & rhs) { return stl_blows.O(rhs.stl_blows); }
	    // bool operator>( const PairSecIter & rhs ) { return (stl_blows < rhs.stl_blows); }
	    RELOP(operator>)
	    RELOP(operator>=)
	    RELOP(operator<)
	    RELOP(operator<=)
	    RELOP(operator==)
	    RELOP(operator!=)
	#undef RELOP

	#define RELOP(O) bool O(const PST & rhs) { return stl_blows.O(rhs); }
	    RELOP(operator>)
	    RELOP(operator>=)
	    RELOP(operator<)
	    RELOP(operator<=)
	    RELOP(operator==)
	    RELOP(operator!=)
	#undef RELOP
    };



#if MMV_TEMPLATE
    template< typename KT, typename VT > class MMVIter :
	    public PairSecIter< std::multimap<KT,VT>::iterator , VT >
	{
	public:
	    MMVIter( std::multimap<KT,VT>::iterator & i ) : PairSecIter( i ) {}
	};
#else
    #define MMVIter(KT,VT)  PairSecIter<std::multimap<KT,VT>::iterator	    ,VT>
    #define MMVCIter(KT,VT) PairSecIter<std::multimap<KT,VT>::const_iterator,VT>
#endif



} // ---- end namespace genepi



/** @} */



#endif // ! __base_PairSecIter_h
