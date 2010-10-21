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
/// \file GeneticDistance.h
/// Definition of the genepi::GeneticDistance class.
//=============================================================================

#ifndef __base_GeneticDistance_h
#define __base_GeneticDistance_h


#include "GeneticDistanceUnit.h"
#include <bclib/exceptions.h>



namespace genepi { // ----



/** \addtogroup base
 * @{ */



//-----------------------------------------------------------------------------
//
// GeneticDistance
/// Class capable of representing distance on the genome in either base-pair
/// (basepair, kilobase, megabase) or Morgan (morgan, centimorgan) units and,
/// within each distance-type, capable of automatically converting between
/// units.
///
/// @sa GeneticDistanceUnit GeneticDistanceUnit.h
///
/// Currently is allowing conversions between base-pairs and Morgans (see
/// MB_PER_MGN): should this be forbidden?
//
//-----------------------------------------------------------------------------

class GeneticDistance
    {
    private:
	double distance ; ///< Morgans or megabases
	bool   isBP	; ///< Is @a distance in megabases (true) or Morgans (false)

	static const double UNLINKED_DIST  = 1.31484755040568800000e+16;
	static const double NEW_CHROM_DIST = 1.41484755040568800000e+16;
	static const double COMPOSITE_DIST = 0.0;

    public:

	/// Perhaps more appropriate would be to forbid conversions between
	/// base-pairs and Morgans, but for the moment we choose an arbitrary
	/// conversion factor
	static const double MB_PER_MGN = 100.0;

	GeneticDistance() {}
	GeneticDistance( GeneticDistanceUnit unit, double val ) { set(unit,val); }

	template < GeneticDistanceUnit > double value() const;

	void set( GeneticDistanceUnit unit, double val );

	double inMorgans     () const { gp_assert( ! isUnlinked	    () );
					gp_assert( ! isNewChromosome() );
					return isBP ? (distance / MB_PER_MGN) : distance; }
	double inCentimorgans() const { return (inMorgans() * 100.0); }
	double inBasepairs   () const { return (inMegabases() * 1e6); }
	double inKilobases   () const { return (inMegabases() * 1e3); }
	double inMegabases   () const { gp_assert( distance != UNLINKED_DIST  );
					gp_assert( distance != NEW_CHROM_DIST );
					return isBP ? distance : (distance * MB_PER_MGN); }

	bool isUnlinked	    () const { return (distance == UNLINKED_DIST ); }
	bool isNewChromosome() const { return (distance == NEW_CHROM_DIST); }
	bool isComposite    () const { return (distance == COMPOSITE_DIST); }

	void makeUnlinked     () { distance = UNLINKED_DIST; }
	void makeNewChromosome() { distance = NEW_CHROM_DIST; }
	void makeComposite    () { distance = COMPOSITE_DIST; }

	bool exceedsThreshold() const
	    { return isBP ? (inMegabases() >= 10.0) : (inMorgans() >= 100.0); }

	#define RELOP(X) bool operator X( const GeneticDistance & rhs ) const \
	    { return isBP ? \
			(inMegabases() X rhs.inMegabases()) : \
			(inMorgans() X rhs.inMorgans()); }
	RELOP(==)
	RELOP(!=)
	RELOP(>=)
	RELOP(<=)
	RELOP(>)
	RELOP(<)
	#undef RELOP
    };

// Specializations of value<> template (there is no generalized version):
template<> inline double GeneticDistance::value<basepairs   >() const { return inBasepairs   (); }
template<> inline double GeneticDistance::value<kilobases   >() const { return inKilobases   (); }
template<> inline double GeneticDistance::value<megabases   >() const { return inMegabases   (); }
template<> inline double GeneticDistance::value<centimorgans>() const { return inCentimorgans(); }
template<> inline double GeneticDistance::value<Morgans	    >() const { return inMorgans     (); }



/** @} */



} // ---- end namespace genepi



#endif // ! __base_GeneticDistance_h
