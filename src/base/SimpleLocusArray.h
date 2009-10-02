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

#include "GeneticDistanceUnit.h"
#include "SimpleLocus.h"
#include "config.h"	// AGGRESSIVE_RANGE_CHECK
#include <bclib/cvector.h>



namespace genepi { // ----



/** \addtogroup base
 * @{ */



typedef size_t SLocIdxType;



//-----------------------------------------------------------------------------
//
/// Container for the SimpleLocus records.  Also holds some summary and
/// aggregate information such as the number of composite loci, number of
/// chromosomes, and the genetic distance unit.
//
//-----------------------------------------------------------------------------


class SimpleLocusArray : public cvector<SimpleLocus>
    {
    typedef cvector<SimpleLocus> SUPER;
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
	    #else
		if ( idx ) {;} // Suppress compiler warning
	    #endif
	    }


    public:

	GeneticDistanceUnit getGDU() const { return gdu; }

	size_t getNComposite   () const { return nComposite  ; } ///< number of composite loci
	size_t getNChromosomes () const { return nChromosomes; } ///< number of chromosomes

	const SimpleLocus & atUnsafe( SLocIdxType idx ) const { return getVector_unsafe().operator[](idx); }


	/// Locate the index of a locus in the array, based on its name.  Throws
	/// an exception if not found.
	SLocIdxType findIndexOf( const std::string & locus ) const;


	SimpleLocusArray();
	~SimpleLocusArray();

    };




/** @} */



} // ---- end namespace genepi



#endif // ! __base_SimpleLocusArray_h
