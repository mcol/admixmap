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
/// \file AlleleArray.cc
/// Implementation of AlleleArray, AlleleFreqTable, AlleleProbTable, AlleleProbVect
//=============================================================================

#include "AlleleArray.h"


#if AA_OSTREAM
    #include <iostream>
    #include <iomanip>
#endif



namespace genepi { // ----



//-----------------------------------------------------------------------------
// normalizeProbs()
//-----------------------------------------------------------------------------

void AlleleProbTable::normalizeProbs()
    {
    const Genotype::AlleleType lim = getLoc().getNumAlleles();

    for ( PopIdx pIdx = 0 ; pIdx < getK() ; ++pIdx )
	{
	double total = 0.0;

	for ( Genotype::AlleleType al = 1 ; al <= getLoc().getNumAlleles() ; ++al )
	    total += at( al, pIdx );

	if ( total == 0.0 )
	    {
	    // If there's no data at all for a given population, use even
	    // probabilities for all possible alleles:
	    const double evenProb = 1.0 / lim;
	    for ( Genotype::AlleleType al = 1 ; al <= getLoc().getNumAlleles() ; ++al )
		at( al, pIdx ) = evenProb;
	    }
	else
	    {
	    for ( Genotype::AlleleType al = 1 ; al <= getLoc().getNumAlleles() ; ++al )
		at( al, pIdx ) /= total;
	    }
	}
    }



//-----------------------------------------------------------------------------
// Print an AlleleProbTable to an ostream
//-----------------------------------------------------------------------------

#if AA_OSTREAM

    void AlleleProbTable::print( std::ostream & os, const std::vector<std::string> & pops ) const
	{
	os << "Allele probability table for locus " << getLoc().getName() << ":\n";

	os.setf(std::ios::fixed);
	os.precision(8);

	os << "Allele";
	for ( PopIdx pIdx = 0 ; pIdx < getK() ; ++pIdx )
	    os << std::setw(18) << pops[pIdx];
	os << '\n';

	for ( Genotype::AlleleType al = 1 ; al <= getLoc().getNumAlleles() ; ++al )
	    {
	    os << std::setw(4) << al << "   ";

	    for ( PopIdx pIdx = 0 ; pIdx < getK() ; ++pIdx )
		os << std::setw(18) << at( al, pIdx );
	    os << '\n';
	    }
	}

#endif // AA_OSTREAM



} // ---- end namespace genepi
