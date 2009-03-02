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



//-----------------------------------------------------------------------------
//
/// Container for the SimpleLocus records.  Also holds some summary and
/// aggregate information such as the number of composite loci and the genetic
/// distance unit.
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

	GeneticDistanceUnit gdu	       ;
	size_t		    nComposite ; ///< number of composite loci


    protected:

	void throwRange( size_t idx ) const __attribute__((noreturn));

	/// Range-check index
	void rc( size_t idx ) const
	    {
	    #if AGGRESSIVE_RANGE_CHECK
		if ( idx >= size() )
		    throwRange( idx );
	    #endif
	    }


    public:

	GeneticDistanceUnit getGDU	  () const { return gdu	      ; }
	size_t		    getNComposite () const { return nComposite; } ///< number of composite loci


	//-------------------------------------------------------------------
	// Vector-style access (see NOTE *1*):
	//-------------------------------------------------------------------

	SimpleLocus &	    operator[]( size_t idx )	   { rc(idx); return SUPER::operator[](idx); }
	const SimpleLocus & operator[]( size_t idx ) const { rc(idx); return SUPER::operator[](idx); }
	size_t		    size()		     const { return SUPER::size(); }

	Iter	  begin()	{ return SUPER::begin(); }
	Iter	  end  ()	{ return SUPER::end  (); }
	ConstIter begin() const { return SUPER::begin(); }
	ConstIter end  () const { return SUPER::end  (); }


	SimpleLocusArray();
	~SimpleLocusArray();

    };




/** @} */



} // ---- end namespace genepi



#endif // ! __base_SimpleLocusArray_h
