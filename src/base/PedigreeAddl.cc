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
/// \file PedigreeAddl.cc
/// Some additional functions for Pedigrees, outside of the class.
//=============================================================================


#include "PedigreeAddl.h"


#include <iostream>
#include <iomanip>


using namespace std;


namespace genepi { // ----



//-----------------------------------------------------------------------------
// print_aggregate_summary()
//-----------------------------------------------------------------------------

ostream & print_aggregate_summary( ostream & os, const vector<Pedigree> & peds )
    {
    if ( peds.size() == 0 )
	return os;

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

    //gp_assert( totNMembers == pedFileArray.getNumOrganisms() );

    const double avgNMembers  = double(totNMembers)  / peds.size();
    const double avgNFounders = double(totNFounders) / peds.size();
    const double avgNSibs     = double(totNSibs)     / peds.size();

    os.setf(ios::fixed);
    os << setprecision(2) << '\n'
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


    return os;
    }




} // ---- end namespace genepi
