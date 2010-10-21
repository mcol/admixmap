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
/// \file SimpleLocusParser.h
/// Definition of the genepi::SimpleLocusParser class.
//=============================================================================

#ifndef __base_SimpleLocusParser_h
#define __base_SimpleLocusParser_h



#include "SimpleLocusArray.h"
#include "GFileLexer.h"



namespace genepi { // ----



/** \addtogroup base
 * @{ */



//-----------------------------------------------------------------------------
/// Class to read and parse a locus file into a SimpleLocusArray
//-----------------------------------------------------------------------------

class SimpleLocusParser : public GFileLexer
    {
    private:
	static const size_t MAX_COLS = 4;

	std::string	   headers[ MAX_COLS ] ;
	SimpleLocusArray & loci		       ;

    public:

	const SimpleLocusArray & getLoci() const { return loci; }
	SimpleLocusArray &	 getLoci()	 { return loci; }


	//---------------------------------------------------------------
	// Constructors/destructor:
	//---------------------------------------------------------------

	/// Constructs the object, but doesn't yet open the file (use parse() )
	SimpleLocusParser( const char * fileName, SimpleLocusArray & output );

	/// Opens the file, parses the data, closes the file
	void parse();

	/// Static convenience method to construct the parser object and parse the file.
	static void parse( const char * fileName, SimpleLocusArray & output )
	    { SimpleLocusParser( fileName, output ) . parse(); }
    };



} // ---- end namespace genepi



/** @} */



#endif // ! __base_SimpleLocusParser_h
