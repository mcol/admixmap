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
/// \file SimpleLocusParser.cc
/// Implementation of genepi::SimpleLocusParser class
//=============================================================================


#include "SimpleLocusParser.h"


#include <cctype>	// isblank()


#define STATUS_TO_COUT	1


#if STATUS_TO_COUT
    #include <iostream>
#endif



static const size_t DISTANCE_COL = 2; // 0-based



namespace genepi { // ----



const size_t SimpleLocusParser::MAX_COLS;



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
//
/// If the file contains the optional chromosome labels, they are used.  If not,
/// they will be assigned default labels (numbered sequentially starting with
/// 1).  Explicitly-labeled and default-labeled chromosomes may @b not be mixed
/// in the same locus-file.
//
//-----------------------------------------------------------------------------

void SimpleLocusParser::parse()
    {

    #if STATUS_TO_COUT
	std::cout << "Loading " << getFileName() << "... ";
    #endif


    int chromLabelCtr = 0; // Only used for default chromosome labeling
    bool explicit_chrom_labels = false; // Initialize to suppress compiler warning


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


	    // --- Column 1: locus name ---
	    row.name = lexString();

	    // It's hard to understand why this would happen and if so why we
	    // would want to strip it, but strip leading/trailing whitespace:
	    stripWhitespace( row.name );


	    // --- Column 2: number of alleles ---
	    const long nAlleles = lexInteger();
	    if ( nAlleles < 2 )
		throwError( "number of alleles is less than 2" );
	    row.numAlleles = nAlleles;


	    // --- Column 3: distance to this locus from previous locus ---
	    { // Begin scope
	    const Token & tok = lexToken();
	    if ( ! tok.isToken() )
		throwError( "too few fields (no distance)" );
	    if ( tok.isType(T_STRING) && (tok.getStrVal() == "NA") )
		{
		row.makeUnlinkedToPrevious();
		row.distance.makeNewChromosome(); // Contradictory?
		}
	    else
		{
		const double dist = tok.asFloat( "distance" );

		if ( dist < 0.0 )
		    throwError( "distance from previous locus is negative" );

		row.distance.set( loci.getGDU(), dist );

		if ( (! row.distance.isUnlinked()) && row.distance.exceedsThreshold() )
		    {
		    if ( dist != 100.0 ) // Backwards-compatibility: no warning
			warn( estr("distance of ") + dist + ' ' +
			    gduAsString(loci.getGDU()) + " at locus " + row.name +
			    " exceeds threshold; starts new chromosome." );
		    row.distance.makeNewChromosome();
		    }
		}
	    } // End scope


	    // --- Column 4 (optional): chromosome label: ---
	    const Token tok = lexToken();
	    if ( ! tok.isToken() ) // EOL or EOF reached
		row.chromNum = -1; // The no-chromosome-label code
	    else
		{
		if ( tok.isType( T_INTEGER ) )
		    row.chromNum = tok.getIntVal();
		else if ( tok.isType( T_STRING ) && (tok.getStrVal().size() == 1) &&
			    ((tok.getStrVal()[0] == 'x') || (tok.getStrVal()[0] == 'X')) )
		    row.chromNum = SimpleLocus::SECRET_X_CHROM_NUM;
		else
		    throwError( estr("invalid chromosome label ") + tok.asString() );

		if ( lexToken().isToken() ) // Not EOL or EOF
		    throwError( "spurious garbage at end of line (too many fields)" );
		}

	    // --- A few more validation checks ---
	    if ( loci.empty() )		// First locus in file:
		{
		if ( row.isCompositeWithPrevious() )
		    throwError( "first locus in file is composite with previous?!?!" );
		++loci.nComposite;

		if ( ! row.startsNewChromosome() )
		    throwError( "first locus in file doesn't start a new chromosome?!?!" );
		++loci.nChromosomes;

		explicit_chrom_labels = row.hasChromLabel();
		}
	    else			// Not the first locus in file:
		{

		// Check for duplicate locus-name:
		for ( SimpleLocusArray::Iter iter = loci.begin() ; iter != loci.end() ; ++iter )
		    if ( iter->name.equalsCaseInsens( row.name ) )
			throwError( estr("Locus name ") + row.name +
				" occurs multiple times in the locus-file" );

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

		if ( row.startsNewChromosome() )
		    ++loci.nChromosomes;
		}

	    if ( row.hasChromLabel() != explicit_chrom_labels )
		throwError( "mixed explicit/default-numbered chromosome labels in the same locus-file" );

	    // Assign default chromosome numbering if necessary:
	    if ( ! row.hasChromLabel() )
		{
		if ( row.startsNewChromosome() )
		    ++chromLabelCtr;
		row.chromNum = chromLabelCtr;
		}

	    // --- Add the newly-parsed row to the container of simple loci. ---
	    loci.push_back( row );

	    }

	}

    // Catch and rethrow DataValidError's just so that we can isolate all other exceptions below
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


    if ( loci.size() == 0 )
	throwError( "No loci in locus-file" );


    #if STATUS_TO_COUT
	std::cout << loci.size() << " simple-loci; "
	    << loci.getNComposite() << " composite loci; "
	    << loci.getNChromosomes() << " chromosomes.";
	if ( getNWarnings() != 0 )
	    std::cout << "  " << getNWarnings() << " warnings.";
	std::cout << '\n';
    #endif

    }



} // ---- end namespace genepi
