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
/// Definition of the genepi::SimpleLocus class.
//=============================================================================

#ifndef __base_SimpleLocus_h
#define __base_SimpleLocus_h



#include "bclib/estr.h"
#include "GeneticDistance.h"



#define START_CHROMOSOME_INDICATOR 1



namespace genepi { // ----


/** \addtogroup base
 * @{ */



/// "Boolean" datatype to indicate whether a given locus is on the X chromosome
/// or some other chromosome.

enum IsXChromType
    {
    CHR_IS_X	 ,
    CHR_IS_NOT_X
    };



class SimpleLocus
    {
    friend class SimpleLocusParser;

    public:
	typedef GeneticDistance GDist; ///< Convenience, shorter

    private:
	estr	name	   ;
	size_t	numAlleles ;
	GDist	distance   ; ///< Distance to the preceding locus
	int	chromNum   ;
	bool	inComposite;

	static const int SECRET_X_CHROM_NUM = -2; // Yuck

    public:
	const std::string & getName	  () const { return name	; }
	size_t		    getNumAlleles () const { return numAlleles	; }
	const GDist &	    getDistance	  () const { return distance	; } ///< Distance to preceding locus
	int		    getChromNum	  () const { return chromNum	; }
	IsXChromType	    isXChrom	  () const { return (chromNum==SECRET_X_CHROM_NUM) ? CHR_IS_X : CHR_IS_NOT_X; }
	std::string	    getChromLabel () const { return (chromNum==SECRET_X_CHROM_NUM) ? "X" : estr(chromNum); }
	bool		    hasChromLabel () const { return (chromNum==SECRET_X_CHROM_NUM) || (chromNum>=0); }

	bool isInComposite() const { return inComposite; }
	void makeInComposite() { inComposite = true; }

	// "Unlinked To Previous" was formerly known as "isMissing"
	bool isUnlinkedToPrevious() const { return distance.isUnlinked(); }
	bool isLinkedToPrevious  () const { return ! distance.isUnlinked(); }
	void makeUnlinkedToPrevious() { distance.makeUnlinked(); }

	bool isCompositeWithPrevious() const { return distance.isComposite(); }
	void makeCompositeWithPrevious() { distance.makeComposite(); }

	#if START_CHROMOSOME_INDICATOR
	    bool startsNewChromosome() const { return distance.isNewChromosome(); }
	    void makeStartsNewChromosome() { distance.makeNewChromosome(); }
	#endif

	std::string getDesc() const { return name + '(' + getChromLabel() + ')'; }
    };



} // ---- end namespace genepi


/** @} */



#endif // ! __base_SimpleLocus_h
