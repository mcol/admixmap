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
/// \file SimpleLocusArray.h
/// Implementation of the SimpleLocusArray class
//=============================================================================


#ifndef __base_SimpleLocusArray_h
#define __base_SimpleLocusArray_h



#include <cstddef>	// size_t
#include <vector>

#include "SimpleLocus.h"
#include "GeneticDistanceUnit.h"
#include "config.h"	// AGGRESSIVE_RANGE_CHECK



namespace genepi { // ----



/** \addtogroup base
 * @{ */



typedef size_t SLocIdxType;



//-----------------------------------------------------------------------------
//
/// Container for the SimpleLocus records.  Also holds some summary and
/// aggregate information such as the number of composite loci, number of
/// chromosomes, and the genetic distance unit.
///
/// <A name="note-1"></A>
/// <TABLE STYLE="border: groove 3pt aqua;">
///
///  <TR>
///	<TD><B>NOTE *1*</B></TD>
///	<TD>
///	We currently use private inheritance from std::vector so that we can
///	impose stricter range-checking; public inheritance would work fine as
///	well.
///	</TD>
///  </TR>
///
/// </TABLE>
//
//-----------------------------------------------------------------------------


class SimpleLocusArray : private std::vector<SimpleLocus> // see NOTE *1*
    {
    typedef std::vector<SimpleLocus> SUPER;
    friend class SimpleLocusParser;


    public:

	typedef iterator	Iter	  ;
	typedef const_iterator	ConstIter ;


    private:

	GeneticDistanceUnit gdu		;
	size_t		    nComposite	; ///< number of composite loci
	size_t		    nChromosomes; ///< number of chromosomes


    protected:

	void throwRange( SLocIdxType idx ) const __attribute__((noreturn));

	/// Range-check index
	void rc( SLocIdxType idx ) const
	    {
	    #if AGGRESSIVE_RANGE_CHECK
		if ( idx >= size() )
		    throwRange( idx );
	    #endif
	    }


    public:

	GeneticDistanceUnit getGDU() const { return gdu; }

	size_t getNComposite   () const { return nComposite  ; } ///< number of composite loci
	size_t getNChromosomes () const { return nChromosomes; } ///< number of chromosomes


	//-------------------------------------------------------------------
	// Vector-style access (see NOTE *1*):
	//-------------------------------------------------------------------

	SimpleLocus &	    operator[]( SLocIdxType idx )	{ rc(idx); return SUPER::operator[](idx); }
	const SimpleLocus & operator[]( SLocIdxType idx ) const { rc(idx); return SUPER::operator[](idx); }
	const SimpleLocus & atUnsafe  ( SLocIdxType idx ) const {	   return SUPER::operator[](idx); }
	SLocIdxType	    size()			  const { return SUPER::size(); }

	/// Convenience method, in case someone is used to calling STL's
	/// range-checked at() (be aware that this overrides a (privately)
	/// inherited non-virtual method).
	const SimpleLocus & at( SLocIdxType idx ) const { return operator[](idx); }


	const SimpleLocus & last() const { rc(0); return SUPER::back(); }

	Iter	  begin()	{ return SUPER::begin(); }
	Iter	  end  ()	{ return SUPER::end  (); }
	ConstIter begin() const { return SUPER::begin(); }
	ConstIter end  () const { return SUPER::end  (); }


	/// Locate the index of a locus in the array, based on its name.  Throws
	/// an exception if not found.
	SLocIdxType findIndexOf( const std::string & locus ) const;


	SimpleLocusArray();
	~SimpleLocusArray();

    };




/** @} */



} // ---- end namespace genepi



#endif // ! __base_SimpleLocusArray_h
