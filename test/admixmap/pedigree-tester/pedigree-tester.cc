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


#include "GenotypeParser.h"
#include "SimpleLocusParser.h"
#include "SimpleLocusArray.h"
#include "Pedigree.h"


using namespace std;
using namespace genepi;


static bool iv_binary = true; // Print IVs in traditional binary format

// Generate all hidden states ahead of time in a vector versus generating and
// printing on-the-fly:
static bool st_in_vector = false;



//-----------------------------------------------------------------------------
// Print an InheritanceVector to an ostream
//-----------------------------------------------------------------------------

inline static char si2char( const InheritanceVector::SegInd & si )
    {
    return (si == InheritanceVector::SI_PATERNAL) ? 'p' : 'm';
    }

inline static char si2digit( const InheritanceVector::SegInd & si )
    {
    return (si == InheritanceVector::SI_PATERNAL) ? '1' : '0';
    }

inline static ostream & operator<<( ostream & os, const InheritanceVector::Bits & b )
    {
    if ( iv_binary )
	return os << si2digit(b.paternal()) << ',' << si2digit(b.maternal());
    else
	return os << char(toupper(si2char(b.paternal()))) << si2char(b.maternal());
    }

ostream & operator<<( ostream & os, const InheritanceVector & iv )
    {
    os << "IV(";

    if ( iv.getNSibs() == 0 )
	os << "-no-sibs-)";
    else
	{
	const size_t limit = iv.getNMembers() - 1;
	for ( size_t sib = iv.getNFounders() ; sib < limit ; ++sib )
	    os << iv.getMember(sib) << ';';
	os << iv.getMember(limit) << ')';
	}

    return os;
    }



//-----------------------------------------------------------------------------
// Print a Haplotype to an ostream
//-----------------------------------------------------------------------------

ostream & operator<<( ostream & os, const Haplotype & hap )
    {
    return os << hap.paternal() << ',' << hap.maternal();
    }



//-----------------------------------------------------------------------------
// Print a "State" to an ostream
//-----------------------------------------------------------------------------

ostream & operator<<( ostream & os, const State & st )
    {
    os << st.getIV() << " {";

    const size_t limit = st.getIV().getNFounders() - 1;
    for ( size_t fIdx = 0 ; fIdx < limit ; ++fIdx )
	os << st.getFounderHap(fIdx) << ';';

    return os << st.getFounderHap(limit) << '}';
    }



//-----------------------------------------------------------------------------
// prState()
//
// Print a "state" to stdout: this is a Pedigree::StateReceiver, used only when
// --in-vector is *not* specified.
//-----------------------------------------------------------------------------

static void prState( const Pedigree &	       ped	       ,
		     size_t		       sLocIdx	       ,
		     const InheritanceVector & iv	       ,
		     const Haplotype *	       founderHapState )
    {
    cout << ped.getSLoci()[sLocIdx].getName() << ' ' << State(iv,founderHapState) << '\n';
    }



//-----------------------------------------------------------------------------
// runTest()
//-----------------------------------------------------------------------------

static void runTest( const char * locusFileName, const char * pedFileName )
    {
    SimpleLocusArray sLocArray;
    SimpleLocusParser::parse( locusFileName, sLocArray );

    GenotypeParser pedFileArray( pedFileName, sLocArray );

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

	if ( ped.getNSibs() < minNSibs )
	    minNSibs = ped.getNSibs();
	if ( ped.getNSibs() > maxNSibs )
	    maxNSibs = ped.getNSibs();
	totNSibs += ped.getNSibs();
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
	;



    //--------------------------------------------------------------
    // Generate and print the founder-haplotype-states/IVs:
    //--------------------------------------------------------------

    for	 ( vector<Pedigree>::const_iterator iter = peds.begin(); iter != peds.end(); ++iter )
	{

	const Pedigree & ped = *iter;

	cout << "\n\nGenerating space of hidden states for pedigree \""
		<< ped.getId() << "\" (" << ped.getNMembers() << " members, "
		<< ped.getNFounders() << " founders, "
		<< ped.getNSibs() << " sibs)...\n\n";


	// Two different ways to accomplish the same thing: generate the full
	// list in a std::vector then print; or print each state on-the-fly as
	// it is generated.
	if ( st_in_vector )
	    {
	    ped.genPossibleStatesInt();
	    for ( SLocIdxType sLocIdx = 0 ; sLocIdx < sLocArray.size() ; ++sLocIdx )
		{
		const std::vector<State> & states = ped.getConsistentStates( sLocIdx );

		cout << "  Simple-locus \"" << sLocArray[ sLocIdx ].getName()
			<< "\" " << states.size() << " states:\n";

		std::vector<State>::const_iterator stIt;
		for ( stIt = states.begin(); stIt != states.end(); ++stIt )
		    cout << "    " << *stIt << '\n';

		cout << '\n';
		}
	    }
	else
	    ped.genPossibleStates( prState );

	}

    }



//-----------------------------------------------------------------------------
// main()
//-----------------------------------------------------------------------------

int main( int argc, const char * const argv [] )
    {

    const char * const * argPtr = argv + 1;
    --argc;


    // Check for --char-iv option, which prints inheritance vectors in a
    // non-standard format:
    if ( (argc != 0) && (strcasecmp(*argPtr,"--char-av") == 0) )
	{
	--argc;
	++argPtr;
	iv_binary = false;
	}


    // Check for --in-vector option, which accumulates the states in a vector
    // first, then prints them out.
    if ( (argc != 0) && (strcasecmp(*argPtr,"--in-vector") == 0) )
	{
	--argc;
	++argPtr;
	st_in_vector = true;
	}


    if ( argc != 2 )
	{
	cerr << "\nUsage: " << argv[0] << " [--char-iv] [--in-vector] <locus-filename> <pedigree-filename>\n"
	    "\t--char-iv prints inheritance vectors in a non-standard format.\n"
	    "\t--in-vector accumulates all of the states in a vector first, then prints them out.\n\n";
	return 1;
	}


    try
	{
	runTest( argPtr[0], argPtr[1] );
	}
    catch ( exception & err )
	{
	cerr << '\n' << argv[0] << ": error: " << err.what() << '\n';
	}

    return 0;
    }
