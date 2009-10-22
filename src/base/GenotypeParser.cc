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
/// \file GenotypeParser.cc
/// Implementation of the GenotypeParser class
//=============================================================================

#pragma implementation
#include "GenotypeParser.h"


#define TREAT_GFILES_AS_PEDFILES 1 ///< Should we build Pedigree objects for non-ped genotype files?
				   ///< But see also Pedigree::generatePedigrees()
#define STATUS_TO_COUT		 1


#if STATUS_TO_COUT
    #include <iostream>
#endif
#include <set>	// Used for pedigree-connected-graph traversal

#include "bclib/estr.h"



//=============================================================================
// Compile-time options for data validation:
//=============================================================================

/// Should we consider male with diploid genotype-data on X chromosome fatal error?
#define MALE_DIPLOID_X_FATAL		0

/// Should male with diploid genotype-data on X chromosome issue warning if the
/// two gametes' alleles are the same?
#define MALE_DIPLOID_X_WARN_HOMO	0

/// Male diploid genotypes on X chromosome should be regarded as missing or haploid?
#define MALE_DIPLOID_X_HOMO_MISSING	0
#define MALE_DIPLOID_X_HETERO_MISSING	0

/// Should an individual with no genotypes be a fatal error?
#define NO_GENOTYPES_FATAL		0


/// Not sure why we would, but we could count the number of loci on the X chromosome
#define COUNT_X_LOCI			0


// If WAIT_DONT_THROW is turned on, when we encounter an error, just print a
// diagnostic to cerr, then throw one error after all processing is finished.
// If not, throw immediately.
#define WAIT_DONT_THROW			1

//=============================================================================



namespace genepi { // ----



//=============================================================================
// Parser class:
//=============================================================================


//----------------------------------------------------------------------------
// lexSex() [protected]
//----------------------------------------------------------------------------

Organism::SexType GenotypeParser::lexSex()
    {
    Organism::SexType rv;

    const Token & tok = lexToken();

    if ( tok.isType( T_INTEGER ) )
	{
	const long s = tok.getIntVal();
	if ( (s < 0) || (s > 2) )
	    throwError( estr("invalid sex code: ") + s );
	rv = static_cast<Organism::SexType>( s );
	}
    else if ( tok.isType(T_STRING) && equalsCaseInsens(tok.getStrVal(),"na") )
	rv = Organism::SEX_UNKNOWN;
    else
	throwError( estr("expected integer or \"NA\" in sex column; found ") + tok.getTypeName() );

    return rv;
    }



//-----------------------------------------------------------------------------
//
// traverse() [static helper function]
//
/// Recursive descent on the call-stack to traverse the pedigree graph in both
/// the parents and children directions, accumulating the complete connected
/// subgraph, which should contain the entire pedigree
//
//-----------------------------------------------------------------------------

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


    // Recursively traverse the children (depth-first), inserting each into the
    // connected set if not already there:
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
			    "\" is duplicated on line numbers " + org.getLineNum()
			    + " and " + dup->getLineNum() );

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
		std::cerr << row.inLineDesc() << ' ' << row.idDesc()
		    << ": WARNING: father at "
		    << father->inLineDesc() << ' ' << father->idDesc()
		    << " sex was previously unknown: marking as male.\n";
		row.father->sex = Organism::SEX_MALE;
		}
	    else if ( father->isFemale() )
		row.throwError( row.idDesc() + "'s father (org:" + father->getOrgId()
				+ ") at " + father->inLineDesc() + " is female." );
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
		std::cerr << row.inLineDesc() << ' ' << row.idDesc()
		    << ": WARNING: mother at "
		    << mother->inLineDesc() << ' ' << mother->idDesc()
		    << " sex was previously unknown: marking as female.\n";
		row.mother->sex = Organism::SEX_FEMALE;
		}
	    else if ( mother->isMale() )
		row.throwError( row.idDesc() + "'s mother (org:" + mother->getOrgId()
				+ ") at " + mother->inLineDesc() + " is male." );
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
// Check the input data for some basic validation.
//-----------------------------------------------------------------------------

void GenotypeParser::validateOrg( const Organism & org ) const
    {

    const SimpleLocusArray & sLoci = getSLoci();

    // Some counts:
    int nMissing    = 0;
    int nNotMissing = 0;

    for ( SLocIdxType locIdx = sLoci.size() ; locIdx-- != 0 ; )
	{

	// Male X-chromosome must be haploid; everything else must be diploid.

	const SimpleLocus & locus  = sLoci[ locIdx ];
	const bool	    isX	   = locus.isXChrom();
	const Genotype &    gType  = org.getGType( locIdx );
	const bool	    isMale = org.isMale();

	if ( gType.isMissing2() )
	    ++nMissing;
	else
	    {
	    ++nNotMissing;

	    if ( gType.isHaploid() )
		{
		if ( ! (isMale && isX) )
		    warn( estr("organism ") + org.idDesc() + " at locus " + locus.getName() +
			    ": haploid genotype data on non-male organism or non-X chromosome" );
		}
	    else
		{
		// Diploid data on male's X chromosome: if homozygous, just
		// convert to haploid; otherwise issue a warning.
		if ( isMale && isX )
		    {
		    if ( gType.isHeterozygous() )
			warn( estr("organism ") + org.idDesc() + " at locus " + locus.getName() +
			    ": diploid heterozygous genotype data on male organism's X chromosome" );
		    else
			gType.forceHaploid();
		    }
		}
	    }

	}


    // The problem with generating this warning here for pedfiles is that these
    // are often valid organisms in the context of pedigrees, typically because
    // they are the parent of an organism that _does_ have observed data; but we
    // don't know at this stage weather an organism has children since the
    // pedigree graphs have not yet been connected.
    if ( (! isPedFile()) && (nNotMissing == 0) )
	{
	#if NO_GENOTYPES_FATAL
	    org.throwError( "no observed genotypes" );
	#else
	    warn( estr("organism ") + org.idDesc() + ": has no observed genotypes." );
	#endif
	}

    }



//-----------------------------------------------------------------------------
// Constructor:
//-----------------------------------------------------------------------------

GenotypeParser::GenotypeParser( const char * fileName, const SimpleLocusArray & sLoci ) :
	GFileLexer ( fileName ) ,
	simpleLoci ( sLoci    )
    {
    Token		t;
    std::vector<string> headers;


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


	bool pedHasOutcome = true;


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

	    // Warning if have data for X chromosomes and no sex supplied:
	    for ( SimpleLocusArray::ConstIter it = simpleLoci.begin(); it != simpleLoci.end(); ++it )
		if ( it->isXChrom() )
		    {
		    warn( "no-sex file-format yet locus-file contains loci for X chromosome." );
		    break;
		    }
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
	else if ( (nCols == (nLoci + 5)) || (nCols == (nLoci + 6)) )
	    {
	    hasSexFlag = false;
	    isPedFileFlag = true;
	    pedHasOutcome = (nCols == (nLoci + 6));
	    gtypeOffset = nCols - nLoci;
	    #if STATUS_TO_COUT
		std::cout << "format: pedigree file";
		if ( ! pedHasOutcome )
		    std::cout << " (no outcomes)";
	    #endif
	    }
	else
	    throwError( estr("incorrect number of columns (") + nCols +
			    ") for number of loci (" + nLoci + ')' );



	//------------------------------------------------------------------
	// Verify that the genotype columns' headers correspond to the
	// locus-names as read in from the locus file:
	//------------------------------------------------------------------

	for ( size_t idx = 0 ; idx < nLoci ; ++idx )
	    if ( ! equalsCaseInsens( headers[ idx + gtypeOffset ], simpleLoci[ idx ].getName() ) )
		warn( estr("genotype column #") + (idx+gtypeOffset) + "'s header ("
			+ headers[ idx + gtypeOffset ] + ") does not match the "
			"corresponding locus name (" + simpleLoci[ idx ].getName() +
			')' );



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

	    // This relies on skipToToken() leaving lastTokenLine as the last
	    // _parsed_ token's line, even though it's been pushbacked.  This
	    // would be clearer if we had a getLastToken() so that skipToToken()
	    // wouldn't need to pushback.
	    row.lineNum = getLastTokenLine();

	    if ( isPedFile() )
		{
		// Maybe just pass in the headers from the header-line as field names:
		row.famId = lexString(); // "pedigree-ID"
		row.orgId = lexInteger( "organism-ID" );

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

		row.sex = lexSex();

		if ( pedHasOutcome )
		    row.outcome = lexInteger( "affected-status" );
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
		    warn( estr("organism-ID ") + row.famId + " is duplicated on line #" + dup->getLineNum() );

		if ( hasSexColumn() )
		    row.sex = lexSex();
		else
		    row.sex = Organism::SEX_UNKNOWN;
		}


	    row.gtypedFlg = false;

	    for ( size_t i = 0 ; i < nLoci ; ++i )
		{
		const SimpleLocus & sLoc = simpleLoci[ i ];
		const size_t nAlleles = sLoc.getNumAlleles();

		Genotype & gtype = row.gtypes[i];

		// Perhaps we should format a nicer error message here should
		// there be insufficient fields on the line:
		gtype = lexGType( sLoc.getName().c_str() );

		// Alleles are numbered from 1, so we use strict inequality comparison:
		if ( gtype.getVal1(0) > nAlleles )
		    throwError( estr("genotype-data for organism ") + row.idDesc() +
			    " in locus " + sLoc.getName() + " first allele value is " +
			    gtype.getVal1() + " which exceeds maximum value of "
			    + nAlleles );

		if ( gtype.isDiploid() && (gtype.getVal2(0) > nAlleles) )
		    throwError( estr("genotype-data for organism ") + row.idDesc() +
			    " in locus " + sLoc.getName() + " second allele value is " +
			    gtype.getVal2() + " which exceeds maximum value of "
			    + nAlleles );

		if ( ! gtype.isMissing2() )
		    row.gtypedFlg = true;
		}


	    if ( hasSexColumn() && (row.sex == Organism::SEX_UNKNOWN) && row.isGenotyped() )
		warn( estr("sex is missing for ") + row.idDesc()
			+ " yet locus-file contains loci for X chromosome." );


	    // That should be the last data on the line:
	    if ( ! lexToken().isType( T_EOL ) )
		throwError( "spurious garbage at end of line (i.e. too many fields)" );


	    validateOrg( row );


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


    if ( rows.size() == 0 )
	throwError( "No organisms in genotype/pedigree file" );


    #if STATUS_TO_COUT
	std::cout << "... " << rows.size() << " organisms.";
	if ( getNWarnings() != 0 )
	    std::cout << "  " << getNWarnings() << " warnings.";
	std::cout << '\n';
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



//-----------------------------------------------------------------------------
// Check the input data for some basic validation.
//-----------------------------------------------------------------------------

void GenotypeParser::validate() const
    {

    const SimpleLocusArray & sLoci = getSLoci();

    #define KEEP_LIST_OF_PLOID_ERRS 0
    #if KEEP_LIST_OF_PLOID_ERRS
	cvector< const SimpleLocusArray * > diploidShouldBeHaploid;
	cvector< const SimpleLocusArray * > haploidShouldBeDiploid;
    #endif

    for ( RowsType::const_iterator it = rows.begin() ; it != rows.end() ; ++it )
	{

	const Organism & org = *it;

	// Some counts:
	int nMissing	= 0;
	int nNotMissing = 0;

	#if KEEP_LIST_OF_PLOID_ERRS
	    diploidShouldBeHaploid.clear();
	    haploidShouldBeDiploid.clear();
	#endif

	for ( SLocIdxType locIdx = sLoci.size() ; locIdx-- != 0 ; )
	    {

	    // Male X-chromosome must be haploid; everything else must be diploid.

	    const SimpleLocus & locus	= sLoci[ locIdx ];
	    const bool		isX	= locus.isXChrom();
	    const Genotype &	gType	= org.getGType( locIdx );
	    const bool		isMale	= org.isMale();

	    if ( gType.isMissing2() )
		++nMissing;
	    else
		{
		++nNotMissing;

		if ( gType.isHaploid() )
		    {
		    if ( ! (isMale && isX) )
			#if KEEP_LIST_OF_PLOID_ERRS
			    haploidShouldBeDiploid.push_back( &locus );
			#else
			    warn( org.inLineDesc() + ": " + org.idDesc() +
				    ": haploid genotype data on non-male organism or non-X chromosome" );
			#endif
		    }
		else
		    {
		    if ( isMale && isX )
			#if KEEP_LIST_OF_PLOID_ERRS
			    diploidShouldBeHaploid.push_back( &locus );
			#else
			    warn( org.inLineDesc() + ": " + org.idDesc() +
				    ": diploid genotype data on male organism's X chromosome" );
			#endif
		    }
		}

	    }


	if ( nNotMissing == 0 )
	    {
	    #if NO_GENOTYPES_FATAL
		org.throwError( "no observed genotypes" );
	    #else
		warn( org.inLineDesc() + ": " + org.idDesc()
			+ ": has no observed genotypes." );
	    #endif
	    }

	} // End organism loop

    }



} // ---- end namespace genepi
