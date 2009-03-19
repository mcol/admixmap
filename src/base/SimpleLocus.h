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
/// \file SimpleLocus.h
//=============================================================================

#ifndef __base_SimpleLocus_h
#define __base_SimpleLocus_h



#include "estr.h"

#define START_CHROMOSOME_INDICATOR 0



namespace genepi { // ----



/** \addtogroup base
 * @{ */



class SimpleLocus
    {
    friend class SimpleLocusParser;

    private:
	estr	name	   ;
	size_t	numAlleles ;
	double	distance   ;
	int	chromNum   ;
	bool	inComposite;

	static const double UNLINKED_DIST  = 1.31484755040568800000e+16;
	static const double NEW_CHROM_DIST = 1.41484755040568800000e+16;

    public:
	const std::string & getName	  () const { return name	; }
	size_t		    getNumAlleles () const { return numAlleles	; }
	double		    getDistance	  () const { return distance	; }
	bool		    hasChrom	  () const { return (chromNum>=0); }
	int		    getChromNum	  () const { return chromNum	; }

	bool isInComposite() const { return inComposite; }
	void makeInComposite() { inComposite = true; }

	// "Unlinked To Previous" was formerly known as "isMissing"
	bool isUnlinkedToPrevious() const { return (distance == UNLINKED_DIST); }
	bool isLinkedToPrevious  () const { return (distance != UNLINKED_DIST); }
	void makeUnlinkedToPrevious() { distance = UNLINKED_DIST; }

	bool isCompositeWithPrevious() const { return isLinkedToPrevious() && (distance == 0.0); }
	void makeCompositeWithPrevious() { distance = 0.0; }

	#if START_CHROMOSOME_INDICATOR
	    bool startsNewChromosome() const { return (distance==NEW_CHROM_DIST); }
	    void makeStartsNewChromosome() { distance = NEW_CHROM_DIST; }
	#endif
    };



} // ---- end namespace genepi



/** @} */



#endif // ! __base_SimpleLocus_h
