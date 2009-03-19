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
/// \file AlleleFreqParser.cc
/// Implementation of the AlleleFreqParser class.
//=============================================================================

#include "AlleleFreqParser.h"



#define STATUS_TO_COUT		1

#if STATUS_TO_COUT
    #include <iostream>
#endif



namespace genepi { // ----



typedef Genotype::AlleleType AlleleType;



//-----------------------------------------------------------------------------
// Constructor:
//-----------------------------------------------------------------------------

AlleleFreqParser::AlleleFreqParser( const char * fileName, const SimpleLocusArray & _loci,
				    PopArray & _outPop, AlleleProbVect & _outProb ) :
	GFileLexer ( fileName ) ,
	loci	   ( _loci    ) ,
	outProb	   ( _outProb ) ,
	populations( _outPop  )
    {
    outProb.reserve( loci.size() );

    // We don't allocate the AlleleProbTable's here because we don't yet know
    // the number of populations (until we parse the header line).
    }



//-----------------------------------------------------------------------------
//  checkLocusComplete() [protected]
/// Helper for parse()
//-----------------------------------------------------------------------------

void AlleleFreqParser::checkLocusComplete( Genotype::AlleleType lastAllele,
					   const SimpleLocus *  lastLocus ) const
    {
    if ( lastAllele != lastLocus->getNumAlleles() )
	warn( estr("locus ") + lastLocus->getName() +
		": contiguous locus block incomplete: ended at " +
		lastAllele + ", but expected up to " + lastLocus->getNumAlleles() );
    }



//-----------------------------------------------------------------------------
// parse()
//-----------------------------------------------------------------------------

void AlleleFreqParser::parse()
    {

    #if STATUS_TO_COUT
	std::cout << "Loading " << getFileName() << "... ";
    #endif


    try
	{

	//------------------------------------------------------------------
	// Skip blank lines, comments, etc:
	//------------------------------------------------------------------
	if ( ! skipToToken() )
	    throwError( "empty allele-frequency file" );



	//------------------------------------------------------------------
	// The first (header) line defines the populations, and may tell us if
	// the file is "allele-labeled":
	//------------------------------------------------------------------

	// Throw away the first column's header, which is for the locus:
	lexString();

	// Look for an optional second column header, which if present will
	// determine the file-format to be allele-labeled:
	const Token & tok = lexToken();
	const bool allele_labeled =
	    tok.isType(T_STRING) && equalsCaseInsens(tok.getStrVal(),"allele");
	if ( ! allele_labeled )
	    pushback( tok );


	// The rest are population-names:
	lexLineStrings( populations );


	//------------------------------------------------------------------
	// We allocate the AlleleProbArray's here (didn't allocate in the
	// constructor because we didn't know the number of populations until we
	// parsed the header line).  Too bad that STL forces us to make a copy
	// rather than construct-in-place at the end of the storage.
	//------------------------------------------------------------------

	const PopIdx K = populations.size();

	for ( SLocIdxType idx = 0 ; idx < loci.size() ; ++idx )
	    outProb.push_back( AlleleProbTable( loci[idx], K ) );


	//------------------------------------------------------------------
	// Read the frequency lines, skipping blank lines and comments.  Since
	// adjacent lines are guaranteed (in non-allele-labeled files) or likely
	// (in allele-labeled files) to be for the same locus, we cache some
	// locus-related data from one line to the next.
	//------------------------------------------------------------------

	// We cache the locus from the previous line and associated data:
	string		    lastLocusName;
	SLocIdxType	    lastLocusIdx= 0; // Initialize to suppress compiler warning
	const SimpleLocus * lastLocus	= 0; // Initialize to suppress compiler warning
	AlleleProbTable *   lastAFTab	= 0; // Initialize to suppress compiler warning
	AlleleType	    lastAllele	= 0; // Initialize to suppress compiler warning

	while ( skipToToken() )
	    {

	    // Column 1: locus name
	    const string & locusName = lexString();


	    AlleleType expectedAllele;
	    if ( locusName == lastLocusName )
		{
		expectedAllele = lastAllele + 1;
		// To-do: format a better error message for allele_labeled files
		if ( expectedAllele > lastLocus->getNumAlleles() )
		    throwError( estr("Too many lines for locus") + locusName );
		}
	    else
		{
		const SLocIdxType lIdx = loci.findIndexOf( locusName );

		if ( lastLocusName.empty() ) // If this is the first line in the file
		    {
		    if ( lIdx != 0 )
			warn( string("first line (") + locusName +
				") is not the first locus in the locus file" );
		    }
		else
		    {
		    checkLocusComplete( lastAllele, lastLocus );

		    if ( lIdx != (lastLocusIdx + 1) )
			warn( estr("missing locus, or loci order differs from locus file: "
				    "jumps from ") + lastLocusName + " to " + locusName );
		    }

		lastLocusName = locusName;
		lastLocusIdx  = lIdx;
		lastAFTab     = &( outProb.at(lIdx) );
		lastLocus     = &( loci[ lIdx ] );

		expectedAllele = 1;
		}


	    if ( allele_labeled )
		{
		const AlleleType thisAllele = lexInteger();

		if ( (thisAllele == 0) ||
			(static_cast<unsigned>(thisAllele) > lastLocus->getNumAlleles()) )
		    throwError( estr("Labeled-allele value ") + thisAllele
				+ " out-of-range for locus " + locusName );

		if ( thisAllele != expectedAllele )
		    warn( estr("labeled alleles out of order: expected ") + expectedAllele
				+ "; found " + thisAllele );
		}


	    lastAllele = expectedAllele;


	    // Maybe should validate here that freqs are non-negative:
	    for ( PopIdx pNum = 0 ; pNum < K ; ++pNum )
		lastAFTab->at( expectedAllele, pNum ) = lexFloat( "allele-population frequency" );

	    if ( lexToken().isToken() )
		throwError( "number of allele-frequencies on line exceeds number of populations" );

	    }


	checkLocusComplete( lastAllele, lastLocus );

	if ( lastLocusIdx != (loci.size() - 1) )
	    warn( string("last locus in allele-freq file (") + lastLocusName +
		    ") is not last in locus file (" + loci.last().getName() + ')' );

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
	std::cout << populations.size() << " populations";
	if ( getNWarnings() != 0 )
	    std::cout << ", " << getNWarnings() << " warnings";
	std::cout << "... normalizing...\n";
    #endif


    for ( SLocIdxType lIdx = 0 ; lIdx < loci.size() ; ++lIdx )
	outProb[lIdx].normalizeProbs();

    }



} // ---- end namespace genepi
