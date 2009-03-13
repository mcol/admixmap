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
/// \file SimpleLocusParser.cc
//=============================================================================

#include "SimpleLocusParser.h"



#define STATUS_TO_COUT	1


#if STATUS_TO_COUT
    #include <iostream>
#endif



static const size_t DISTANCE_COL = 2; // 0-based



namespace genepi { // ----



//-----------------------------------------------------------------------------
// Constructor:
//-----------------------------------------------------------------------------

SimpleLocusParser::SimpleLocusParser( const char * fileName, SimpleLocusArray & output ) :
	GFileLexer ( fileName ) ,
	loci	   ( output   )
    {
    }



//-----------------------------------------------------------------------------
// parse()
//-----------------------------------------------------------------------------

void SimpleLocusParser::parse()
    {

    #if STATUS_TO_COUT
	std::cout << "Loading " << getFileName() << "... ";
    #endif


    // Set up a try-block so that we can catch exceptions that aren't
    // DataValidError's, and wrap them in a DataValidError exception so that
    // filename/line-number are remembered and reported back to the user.
    try
	{

	//------------------------------------------------------------------
	// Skip blank lines, comments, etc:
	//------------------------------------------------------------------
	if ( ! skipToToken() )
	    throwError( "empty locus file" );


	//------------------------------------------------------------------
	// Read the header line:
	//------------------------------------------------------------------

	for ( size_t col = 0 ; ; ++col )
	    {
	    const Token & tok = lexToken();
	    if ( tok.isToken() )
		{
		if ( col == MAX_COLS )
		    throwError( estr("too many columns in header line (max ") + MAX_COLS + ')' );
		headers[ col ] = tok.asString();
		}
	    else
		{
		if ( col < DISTANCE_COL )
		    throwError( estr("Too few (") + col + ") columns in header line"
				    " (require at least " + DISTANCE_COL + ')' );
		break;
		}
	    }

	try
	    {
	    loci.gdu = gduFromString( headers[DISTANCE_COL].c_str() );
	    }
	catch ( std::runtime_error & e )
	    {
	    loci.gdu = Morgans; // Default, usual for admixture mapping
	    }


	//------------------------------------------------------------------
	// Prepare a "buffer" SimpleLocus object, to read each row into before
	// inserting in the vector.
	//------------------------------------------------------------------

	SimpleLocus row;
	row.inComposite = false;


	//------------------------------------------------------------------
	// Read the data lines, skipping blank lines and comments:
	//------------------------------------------------------------------

	while ( skipToToken() )
	    {

	    // Column 1: locus name
	    row.name = lexString();

	    // Column 2: number of alleles
	    const long nAlleles = lexInteger();
	    if ( nAlleles < 1 )
		throwError( "number of alleles is less than 1" );
	    row.numAlleles = nAlleles;

	    // Column 3: distance to this locus from previous locus
	    #if 0 // This would be way too simple:
		row.distance = lexFloat();
	    #else
		{ // Begin scope
		const Token & tok = lexToken();
		if ( tok.isType(T_STRING) && (tok.getStrVal() == "NA") )
		    row.makeUnlinkedToPrevious();
		else
		    row.distance = tok.asFloat();
		} // End scope
	    #endif

	    // Column 4 (optional): chromosome label:
	    const Token tok = lexToken();
	    if ( ! tok.isToken() ) // EOL or EOF reached
		row.chromNum = -1; // The no-chromosome-label code
	    else
		{
		if ( tok.isType( T_INTEGER ) )
		    row.chromNum = tok.getIntVal();
		else if ( tok.isType( T_STRING ) && (tok.getStrVal().size() == 1) &&
			    ((tok.getStrVal()[0] == 'x') || (tok.getStrVal()[0] == 'X')) )
		    row.chromNum = -2; // The X-chromosome code
		else
		    throwError( estr("invalid chromosome label ") + tok.asString() );

		if ( lexToken().isToken() ) // Not EOL or EOF
		    throwError( "spurious garbage at end of line (too many fields)" );
		}


	    if ( loci.empty() )
		{
		if ( row.isCompositeWithPrevious() )
		    throwError( "first locus in file is composite with previous" );
		++loci.nComposite;
		}
	    else
		{

		// Check for duplicate locus-name:
		for ( SimpleLocusArray::Iter iter = loci.begin() ; iter != loci.end() ; ++iter )
		    if ( iter->name.equalsCaseInsens( row.name ) )
			throwError( estr("Locus name ") + row.name +
				    " occurs multiple times in locus-file" );

		// Check for composite locus:
		if ( row.isCompositeWithPrevious() )
		    {
		    // Increment the number of composite loci if appropriate:
		    if ( ! loci.back().isInComposite() )
			loci.back().makeInComposite();
		    row.makeInComposite();
		    }
		else
		    ++loci.nComposite;
		}


	    // Add the newly-parsed row to the container of simple loci.
	    loci.push_back( row );

	    }
	}

    // Catch and rethrow DataValidError's just so that we can isolate other exceptions below:
    catch ( DataValidError & err )
	{
	throw;
	}

    // Wrap other exceptions in a DataValidError exception so that
    // filename/line-number are remembered and reported back to the user.
    catch ( std::exception & err )
	{
	throwError( err.what() );
	}

    #if STATUS_TO_COUT
	std::cout << loci.size() << " simple-loci.\n";
    #endif

    }



} // ---- end namespace genepi
