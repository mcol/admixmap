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
/// \file GenotypeParser.cc
/// Implementation of the GenotypeParser class
//=============================================================================

#pragma implementation
#include "GenotypeParser.h"


#define TREAT_GFILES_AS_PEDFILES 1 ///< Should we build Pedigree objects for Genotype files?
#define STATUS_TO_COUT		 1


#if STATUS_TO_COUT
    #include <iostream>
#endif
#include <set>	// Used for pedigree-connected-graph traversal

#include "estr.h"



namespace genepi { // ----



//=============================================================================
// Parser class:
//=============================================================================


//----------------------------------------------------------------------------
// sexFromInt() [protected]
//----------------------------------------------------------------------------

Organism::SexType GenotypeParser::sexFromInt( long s )
    {
    if ( (s < 0) || (s > 2) )
	throwError( estr("invalid sex code: ") + s );

    return static_cast<Organism::SexType>( s );
    }



//-----------------------------------------------------------------------------
// traverse() [static helper function]
//
//-----------------------------------------------------------------------------

/// Recursive descent on the call-stack to traverse the pedigree graph in both
/// the parents and children directions, accumulating the complete connected
/// subgraph, which should contain the entire pedigree
static void traverse( const Organism & org, std::set<const Organism *> & connected )
    {
    // This needs more work: we should definitely not have cycles if we just
    // traverse the parent or child digraphs alone, and this should be enforced,
    // but will presumably require another traversal or pass.
    // Check for cycles:
    if ( connected.find(&org) != connected.end() ) // STL should have: connected.contains(&org)
	{
	// Foo...
	return; // ***** RETURN HERE *****
	}


    // Put 'self' into the connected set:
    connected.insert( &org );


    // Traverse the children (depth-first), inserting into the connected set:
    const Organism::ChConstIter limit = org.childrenEnd();
    for ( Organism::ChConstIter iter = org.childrenBegin() ; iter != limit ; ++iter )
	traverse( **iter, connected );


    // Traverse the parents (depth-first), inserting into the connected set:
    if ( org.hasFather() )
	traverse( *org.getFather(), connected );
    if ( org.hasMother() )
	traverse( *org.getMother(), connected );
    }



//-----------------------------------------------------------------------------
//
// buildAndCheckPedGraphs() [private]
//
/// In the case of pedigree-files, after parsing the pedfile to build the rows
/// array of Individuals, we make a second, third, and fourth pass through the
/// array, building the std::multimap of pedigree-to-individuals, and
/// dereferencing the father and mother IDs into pointers to Individual objects,
/// respectively.  We wait until this time to keep a pointer in the pedigree
/// std::multimap, since we need a persistent pointer, and objects in the
/// std::vector may be relocated when elements are added.
//
//-----------------------------------------------------------------------------

void GenotypeParser::buildAndCheckPedGraphs()
    {

    //------------------------------------------------------------------
    // Second pass: build the pmap pedigreeID-to-individual mapping, while
    // checking for duplicate family/individual IDs.
    //------------------------------------------------------------------

    for ( RowIter iter = rows.begin(); iter != rows.end(); ++iter )
	{
	Organism & org = *iter;

	// Check for the uniqueness of the organism-ID within the family:
	const Organism * const dup = findByIdIfExists( org.famId, org.orgId );
	if ( dup != 0 )
	    org.throwError( estr("ID \"") + org.idDesc() +
			    "\" is duplicated on line #" + dup->getLineNum() );

	insert_pmap( org );
	}
    gp_assert_eq( pmap.size(), rows.size() );



    //------------------------------------------------------------------
    // From this point forward (and not before), the following methods are valid:
    //	    findOrgsInPed() [several overloaded versions]
    //	    findById()
    //	    findByIdIfExists()
    //	    exists()
    //------------------------------------------------------------------


    //------------------------------------------------------------------
    // Third pass: resolve parent pointers, checking for various validation
    // errors; build each parent's children-pointer list (i.e. the reverse
    // direction of the digraph from the child-to-parent reference presented in
    // the input data).
    //------------------------------------------------------------------

    for ( RowIter iter = rows.begin(); iter != rows.end(); ++iter )
	{

	Organism & row = *iter;

	if ( row.getFatherId() == MISSING_PARENT_ID )
	    row.father = 0;
	else
	    {
	    const Organism * const father = findByIdIfExists( row.getFamId(), row.getFatherId() );

	    if ( father == 0 )
		row.throwError( estr("father (") + row.getFatherId() +
			    ") of " + row.idDesc() +
			    " does not exist in input pedigree-file" );

	    row.father = const_cast<Organism*>( father );
	    row.father->children.push_front( &row );

	    if ( ! father->sexKnown() )
		{
		std::cerr << row.inLineDesc() << ": WARNING: father at " << father->inLineDesc()
			<< " sex was previously unknown: marking as male";
		row.father->sex = Organism::SEX_MALE;
		}
	    else if ( father->isFemale() )
		row.throwError( estr("father of ") + row.idDesc() +
			    " (" + father->getOrgId() +
			    ") at " + father->inLineDesc() + " is female" );
	    }

	if ( row.getMotherId() == MISSING_PARENT_ID )
	    row.mother = 0;
	else
	    {
	    const Organism * const mother = findByIdIfExists( row.getFamId(), row.getMotherId() );

	    if ( mother == 0 )
		row.throwError( estr("mother (") + row.getMotherId() +
			    ") of " + row.idDesc() +
			    " does not exist in input pedigree-file" );

	    row.mother = const_cast<Organism*>( mother );
	    row.mother->children.push_front( &row );

	    if ( ! mother->sexKnown() )
		{
		std::cerr << row.inLineDesc() << ": WARNING: mother at " << mother->inLineDesc()
		    << " sex was previously unknown: marking as female";
		row.mother->sex = Organism::SEX_FEMALE;
		}
	    else if ( mother->isMale() )
		row.throwError( row.idDesc() + "'s mother " + mother->getOrgId()
				+ " at " + mother->inLineDesc() + " is male" );
	    }
	}


    //------------------------------------------------------------------
    // From this point forward (and not before), we have connected the
    // pointers which form the digraph for each family; thus the following
    // methods are valid to be called:
    //	    Organism::getFather()
    //	    Organism::getMother()
    //	    Organism::getChildren()
    //	    Organism::isFounder()
    //	    etc.
    //------------------------------------------------------------------


    //------------------------------------------------------------------
    // Fourth pass: every pedigree should consist of a digraph which:
    //	    (a) is connected, and
    //	    (b) contains no cycles
    // The "special" case of a single organism with no other family
    // members also meets these conditions.
    //
    // By iterating over the organisms using the pmap rather than rows
    // container, we will visit them all in such an order that all of the
    // members of a pedigree (family) will be adjacent to each other.
    //
    // When we reach the first organism in a new family, we do a depth-first
    // search of the pedigree digraph (both parents and children), collecting
    // each visited node in a set, which will build the connected sub-graph.
    // Any subsequent organisms in that pedigree which are not in the set of
    // connected nodes are orphans or members of a disjoint connected sub-graph.
    //------------------------------------------------------------------

    std::set<const Organism*> connected;
    const ConstPedIter limit( pmap.end() );
    for ( ConstPedIter iter( pmap.begin() ); iter != limit ; ++iter )
	{
	const Organism & org = **iter;

	// Determine if this organism is in the same family as the last
	// one, or if we're starting a new family:

	if ( connected.empty() || ((*connected.begin())->getFamId() != org.getFamId()) )
	    {
	    // Start new family:
	    connected.clear();

	    // Depth-first traverse the children and parent graphs, adding
	    // to the connected set.  For the moment, we will use the
	    // call-stack to perform the recursion:
	    traverse( org, connected );
	    }
	else
	    {
	    // Another member of existing family:
	    if ( connected.find(&org) == connected.end() ) // STL should have: ! connected.contains(&org)
		org.throwError( estr("family has disjoint connected sub-graphs: ") + org.idDesc() );
	    }
	}

    }



//-----------------------------------------------------------------------------
// Constructor:
//-----------------------------------------------------------------------------

GenotypeParser::GenotypeParser( const char * fileName, const SimpleLocusArray & sLoci ) :
	GFileLexer ( fileName ) ,
	simpleLoci ( sLoci    )
    {
    Token t;


    #if STATUS_TO_COUT
	std::cout << "Loading " << fileName << "... ";
    #endif


    const size_t nLoci = getNSimpleLoci();


    // Set up a try-block during the parse so that we can catch exceptions that
    // aren't DataValidError's, and wrap them in a DataValidError exception so
    // that filename/line-number are remembered and reported back to the user.
    try
	{

	//------------------------------------------------------------------
	// Skip blank lines, comments, etc:
	//------------------------------------------------------------------
	if ( ! skipToToken() )
	    throwError( "empty genotype file" );


	//------------------------------------------------------------------
	// Read the header line:
	//------------------------------------------------------------------

	lexLineStrings( headers );


	//------------------------------------------------------------------
	// Check the number of columns and categorize the file format:
	//------------------------------------------------------------------

	const size_t nCols = headers.size();
	if ( nCols == (nLoci + 1) )
	    {
	    hasSexFlag = false;
	    isPedFileFlag = false;
	    gtypeOffset = 1;
	    #if STATUS_TO_COUT
		std::cout << "format: no-sex genotype file";
	    #endif
	    }
	else if ( nCols == (nLoci + 2) )
	    {
	    hasSexFlag = true;
	    isPedFileFlag = false;
	    gtypeOffset = 2;
	    #if STATUS_TO_COUT
		std::cout << "format: has-sex genotype file";
	    #endif
	    }
	else if ( nCols == (nLoci + 6) )
	    {
	    hasSexFlag = false;
	    isPedFileFlag = true;
	    gtypeOffset = 6;
	    #if STATUS_TO_COUT
		std::cout << "format: pedigree file";
	    #endif
	    }
	else
	    throwError( estr("incorrect number of columns (") + nCols +
			    ") for number of loci (" + nLoci + ')' );



	//------------------------------------------------------------------
	// Prepare a "buffer" Organism object, to read each row into before
	// inserting in the vector.
	// TO-DO: currently copies the loci data; should insert the new row into the
	// vector with 0 loci, allocate in place, then read the loci from the file
	// directly into their final home.
	//------------------------------------------------------------------

	Organism row( *this );
	if ( ! isPedFileFlag )
	    {
	    // Since we're using a combined parser/container for old-style
	    // "genotype" as well as "pedigree" file formats, we store the
	    // genotype-file's "organism ID" in the "family ID" field, and use
	    // a fixed "organism ID" in that case:
	    row.orgId = NON_PEDFILE_ORG_ID;

	    // In non-pedfiles, every organism is a founder:
	    row.father = 0;
	    row.mother = 0;
	    row.depth = 0;

	    row.outcome = 0;

	    if ( ! hasSexColumn() )
		row.sex = Organism::SEX_UNKNOWN;
	    }


	//------------------------------------------------------------------
	// Read the data lines, skipping blank lines and comments:
	//------------------------------------------------------------------

	while ( skipToToken() )
	    {

	    row.lineNum = curLineNum();

	    if ( isPedFile() )
		{
		// Maybe just pass in the headers from the header-line as field names:
		row.famId = lexString();
		row.orgId = lexInteger( "pedigree-member-ID" );

		// For the first pass (reading in the file), we store the father/mother
		// IDs in the place of the pointers (punning); a second pass
		// dereferences them.
		row.setFatherId( lexInteger( "father-ID" ) );
		row.setMotherId( lexInteger( "mother-ID" ) );

		// Per http://linkage.rockefeller.edu/soft/linkage/sec2.7.html,
		// parent-ID of 0 indicates founder:
		if ( (row.getFatherId() == 0) && (row.getMotherId() == 0) )
		    row.depth = 0;
		else
		    row.depth = Organism::UNKNOWN_DEPTH;

		row.sex = sexFromInt( lexInteger( "sex" ) );

		row.outcome = lexInteger();
		}
	    else
		{
		row.famId = lexString();


		// The organism-ID (here called "family ID") in a
		// genotype-file is essentially unused; yet we will check here
		// for its uniqueness within the file and print a warning if
		// duplicates are found:
		const Organism * const dup = findByIdIfExists( row.famId, row.orgId );
		if ( dup != 0 )
		    std::cerr << getFileName() << ':' << curLineNum() << ": warning: organism-ID "
			    << row.famId << " is duplicated on line " << dup->getLineNum() << std::endl;


		if ( hasSexColumn() )
		    row.sex = sexFromInt( lexInteger( "sex" ) );
		else
		    row.sex = Organism::SEX_UNKNOWN;
		}


	    // Perhaps we should format a nicer error message here, in case
	    // there are insufficient fields on the line:
	    for ( size_t i = 0 ; i < nLoci ; ++i )
		row.gtypes[i] = lexGType();


	    // That should be the last data on the line:
	    if ( ! lexToken().isType( T_EOL ) )
		throwError( "spurious garbage at end of line (i.e. too many fields)" );


	    // Add the newly-parsed row to the container of organisms.  If we
	    // are storing them in a std::vector, we cannot yet get a persistent
	    // pointer, as these objects may be relocated as the vector expands
	    // its storage.
	    rows.push_back( row );
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
	std::cout << "... " << rows.size() << " organisms.\n";
    #endif

    #if ! TREAT_GFILES_AS_PEDFILES
      if ( isPedFile() )
    #endif
	buildAndCheckPedGraphs();

    }



//-----------------------------------------------------------------------------
// Destructor:
//-----------------------------------------------------------------------------

GenotypeParser::~GenotypeParser()
    {
    }



//=============================================================================
// Access to contents (i.e. container class):
//=============================================================================

const Organism & GenotypeParser::
	findById( const FamIdType & famId, const OrgIdType & orgId ) const
    {
    ConstPedIter end;
    for ( ConstPedIter iter = findOrgsInPed(famId,end) ; iter != end ; ++iter )
	{
	const Organism & org = **iter;
	gp_assert_eq( org.getFamId(), famId );
	if ( org.getOrgId() == orgId )
	    return org;
	}

    throw std::runtime_error( estr("{fam:") + famId + ";ind:" + orgId + "} not found" );
    }


const Organism * GenotypeParser::
		findByIdIfExists( const FamIdType & famId, const OrgIdType & orgId ) const
    {
    ConstPedIter end;
    for ( ConstPedIter iter = findOrgsInPed(famId,end) ; iter != end ; ++iter )
	{
	const Organism & org = **iter;
	gp_assert_eq( org.getFamId(), famId );
	if ( org.getOrgId() == orgId )
	    return &org;
	}

    return 0;
    }



} // ---- end namespace genepi
