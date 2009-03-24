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
/// \file pedigree-tester.cc
//=============================================================================


#include <iostream>
#include <iomanip>
#include <string>
#include <cctype>	// toupper()
#include <cstring>	// strcasecmp()
#include <cstdlib>	// exit()


#include "AlleleArray.h"
#include "AlleleFreqParser.h"
#include "GenotypeParser.h"
#include "HiddenStateSpace.h"
#include "InheritanceVector.h"
#include "Pedigree.h"
#include "SimpleLocusArray.h"
#include "SimpleLocusParser.h"


using namespace std;
using namespace genepi;



static bool print_afreqs = false;



//-----------------------------------------------------------------------------
// Print a Haplotype to an ostream
//-----------------------------------------------------------------------------

ostream & operator<<( ostream & os, const Haplotype & hap )
    {
    return os << hap.paternal() << ',' << hap.maternal();
    }



//-----------------------------------------------------------------------------
// runTest()
//-----------------------------------------------------------------------------

static void runTest( const char * locusFileName, const char * pedFileName, const char * afFileName )
    {
    SimpleLocusArray loci;
    SimpleLocusParser::parse( locusFileName, loci );

    GenotypeParser pedFileArray( pedFileName, loci );

    if ( pedFileArray.getNumOrganisms() == 0 )
	throw std::runtime_error( string("No organisms in ") + pedFileName );

    vector<Pedigree> peds;
    Pedigree::generatePedigrees( pedFileArray, peds );

    size_t minNMembers	= INT_MAX;
    size_t maxNMembers	= 0;
    size_t totNMembers	= 0;
    size_t minNFounders = INT_MAX;
    size_t maxNFounders = 0;
    size_t totNFounders = 0;
    size_t minNSibs	= INT_MAX;
    size_t maxNSibs	= 0;
    size_t totNSibs	= 0;

    for	 ( vector<Pedigree>::const_iterator iter = peds.begin(); iter != peds.end(); ++iter )
	{
	const Pedigree & ped = *iter;

	if ( ped.getNMembers() < minNMembers )
	    minNMembers = ped.getNMembers();
	if ( ped.getNMembers() > maxNMembers )
	    maxNMembers = ped.getNMembers();
	totNMembers += ped.getNMembers();

	if ( ped.getNFounders() < minNFounders )
	    minNFounders = ped.getNFounders();
	if ( ped.getNFounders() > maxNFounders )
	    maxNFounders = ped.getNFounders();
	totNFounders += ped.getNFounders();

	if ( ped.getNNonFndrs() < minNSibs )
	    minNSibs = ped.getNNonFndrs();
	if ( ped.getNNonFndrs() > maxNSibs )
	    maxNSibs = ped.getNNonFndrs();
	totNSibs += ped.getNNonFndrs();
	}

    gp_assert( totNMembers == pedFileArray.getNumOrganisms() );

    const double avgNMembers  = double(totNMembers)  / peds.size();
    const double avgNFounders = double(totNFounders) / peds.size();
    const double avgNSibs     = double(totNSibs)     / peds.size();

    cout.setf(ios::fixed);
    cout << setprecision(2) << '\n'
	<< peds.size() << " pedigrees:\n"
	"  smallest number of members:  " << setw(3) << minNMembers	<< "\n"
	"  largest number of members:   " << setw(3) << maxNMembers	<< "\n"
	"  average number of members:   " << setw(6)
						     << avgNMembers	<< "\n"
	"  smallest number of founders: " << setw(3) << minNFounders	<< "\n"
	"  largest number of founders:  " << setw(3) << maxNFounders	<< "\n"
	"  average number of founders:  " << setw(6)
						     << avgNFounders	<< "\n"
	"  smallest number of sibs:     " << setw(3) << minNSibs	<< "\n"
	"  largest number of sibs:      " << setw(3) << maxNSibs	<< "\n"
	"  average number of sibs:      " << setw(6)
						     << avgNSibs	<< "\n"
	"\n";



    //--------------------------------------------------------------
    // Read in the allele-frequency file:
    //--------------------------------------------------------------

    AlleleProbVect afVect;
    vector<string> populations;

    AlleleFreqParser::parse( afFileName, loci, populations, afVect );

    if ( print_afreqs )
	for ( SLocIdxType sLocIdx = 0 ; sLocIdx < loci.size() ; ++sLocIdx )
	    afVect[sLocIdx].print( std::cout << '\n', populations );


    //--------------------------------------------------------------
    // Generate and print the emission probabilities founder-haplotype-states/IVs:
    //--------------------------------------------------------------

    for	 ( vector<Pedigree>::const_iterator iter = peds.begin(); iter != peds.end(); ++iter )
	{

	const Pedigree & ped = *iter;

	cout << "\n\nGenerating space of hidden states for pedigree \""
		<< ped.getId() << "\" (" << ped.getNMembers() << " members, "
		<< ped.getNFounders() << " founders, "
		<< ped.getNNonFndrs() << " non-founders)...\n\n"

		"  For this pedigree, AVs are listed in founder organism order:";
	for ( Pedigree::Iterator it = ped.getFirstFounder() ; it != ped.getEndFounder() ; ++it )
	    cout << ' ' << (*it)->getOrgId();
	cout << "\n"
		"  IVs are listed in non-founder organism order:";
	for ( Pedigree::Iterator it = ped.getFirstNonFndr() ; it != ped.getEndNonFndr() ; ++it )
	    cout << ' ' << (*it)->getOrgId();
	cout << '\n';


	ped.genPossibleStatesInternal( populations.size(), afVect );


	for ( SLocIdxType sLocIdx = 0 ; sLocIdx < loci.size() ; ++sLocIdx )
            {
	    const HiddenStateSpace & space = ped.getStateProbs( sLocIdx );

            cout << "\nSimple-locus \"" << loci[ sLocIdx ].getName() << "\":\n";

	    size_t n_states = 0;

	    for ( HiddenStateSpace::Iterator it(space); it; ++it )
		{
		++n_states;
                const HiddenStateSpace::Iterator::State & st = *it;
                cout << "  " << st.av << ' ' << st.iv << ' ' << st.emProb << '\n';
		}

            cout << "\n  (" << n_states << " states).\n\n";
            }

	}

    }



//-----------------------------------------------------------------------------
// usage()
//-----------------------------------------------------------------------------

static void usage( const char * prog ) __attribute__((noreturn));
static void usage( const char * prog )
    {
    cerr << "\nUsage: " << prog << " [options] <locus-filename> <pedigree-filename> <allele-freq-filename>\n"
	"\tOptions:\n"
	"\t--print-afreqs makes a print-out of the allele frequencies before generating probabilities.\n"
	"\t--char-iv prints inheritance vectors in a non-standard format.\n"
	"\t--debug-emission prints debugging information about emission-probability calculations to stdout.\n"
	"\t--debug-recursion prints debugging information about pedigree-graph recursion to stdout.\n"
	"\n";
    exit(1);
    }



//-----------------------------------------------------------------------------
// main()
//-----------------------------------------------------------------------------

int main( int argc, const char * const argv [] )
    {

    const char * const * argPtr = argv + 1;
    --argc;


    while ( (argc != 0) && (**argPtr == '-') )
	{

	const char * const option = *argPtr;


	// --char-iv option: prints inheritance vectors in a non-standard format
	if ( strcasecmp(option,"--char-iv") == 0 )
	    {
	    --argc;
	    ++argPtr;
	    setIVOutputStyle( IV_ALPHA );
	    continue;
	    }


	// --print-afreqs option: prints the allele frequencies
	if ( strcasecmp(option,"--print-afreqs") == 0 )
	    {
	    --argc;
	    ++argPtr;
	    print_afreqs = true;
	    continue;
	    }


	if ( strcasecmp(option,"--debug-recursion") == 0 )
	    {
	    --argc;
	    ++argPtr;
	    Pedigree::dbgRecursion( true );
	    continue;
	    }

	if ( strcasecmp(option,"--debug-emission") == 0 )
	    {
	    --argc;
	    ++argPtr;
	    Pedigree::dbgEmission( true );
	    continue;
	    }


	cerr << argv[0] << ": unknown option " << option << '\n';
	usage( argv[0] );

	}


    if ( argc != 3 )
	{
	cerr << argv[0] << ": wrong number of arguments; require 3 file-names.\n";
	usage( argv[0] );
	}


    try
	{
	runTest( argPtr[0], argPtr[1], argPtr[2] );
	}
    catch ( exception & err )
	{
	cerr << '\n' << argv[0] << ": error: " << err.what() << '\n';
	}

    return 0;
    }
