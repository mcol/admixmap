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
/// \file Genotype.h
/// Defines Genotype, GenotypeArray, and Haplotype
//=============================================================================


#ifndef __base_Genotype_h
#define __base_Genotype_h



#include <climits>  // USHRT_MAX
#include <vector>

#include "estr.h"
#include "exceptions.h"



namespace genepi { // ----



/** \addtogroup base
 * @{ */



//---------------------------------------------------------------------
/// Class to hold a "genotype", two unphased alleles.  This should
/// perhaps be moved out of this class into its own header.  Since the
/// pair of integers are unordered, we will assure that val2 is missing
/// before val1 is.
//---------------------------------------------------------------------

class Genotype
    {
    friend class GFileLexer;

    public:
	typedef unsigned short AlleleType;

	static const AlleleType HAPLOID_VAL = USHRT_MAX ;
	static const AlleleType MISSING_VAL = 0		;

    private:
	AlleleType val1;
	AlleleType val2;

    public:
	AlleleType getVal1() const { gp_assert(val1!=MISSING_VAL); return val1; }
	AlleleType getVal2() const { gp_assert(val2!=MISSING_VAL); return val2; }

	AlleleType getVal1( AlleleType missingVal ) const { return (val1==MISSING_VAL) ? missingVal : val1; }
	AlleleType getVal2( AlleleType missingVal ) const { return (val2==MISSING_VAL) ? missingVal : val2; }

	/// Did the input value contain a delimiter?
	bool isDiploid () const { return (val2!=HAPLOID_VAL); }

	/// Was the input value a single integer?
	bool isHaploid () const { return (val2==HAPLOID_VAL); }

	/// Is one of the two values missing or 0?
	bool isMissing1() const { return (val2==MISSING_VAL); }

	/// Are both of the two values missing or 0?
	bool isMissing2() const { return (val1==MISSING_VAL); }

	/// Is there exactly one observed value (haploid or one missing)?
	bool hasOneVal() const { return (val1!=MISSING_VAL) &&
				((val2==MISSING_VAL) || (val2==HAPLOID_VAL)); }

	/// Is there exactly two observed values?
	bool hasTwoVals() const { return (val1!=MISSING_VAL) &&
				((val2!=MISSING_VAL) || (val2!=HAPLOID_VAL)); }

	/// Is at least one the two values neither missing nor 0?
	bool isObserved() const { return (val1!=MISSING_VAL); }

	/// Is al present as either of the observed alleles?
	bool contains( AlleleType al ) const { return ((val1 == al) || (val2 == al)); }

	/// Is the presence of al consistent with the observed alleles?
	bool consistent( AlleleType al ) const
	    { return (val1 == MISSING_VAL) || (val1 == al) || (val2 == al); }

	/// Are the potential and observed genotypes consistent?
	// ---- How to deal with haploid? ----
	bool consistent( const Genotype & hap ) const
	    { return ( // Perhaps this can be rewritten as a shorter boolean expression?
			(val1 == MISSING_VAL) ||
			((val1 == hap.val1) && (val2 == hap.val2)) ||
			((val2 == hap.val1) && (val1 == hap.val2)) ||
			((val1 == hap.val1) && (val2 == MISSING_VAL))
		     ); }

	/// Create a human-readable representation of the value
	estr desc() const;

	/// This somewhat ugly hack was created especially for the
	/// male-has-diploid-X-chromosome backwards-compatability hack
	/// in AdmixmapGenotypeConverter.cc
	void forceHaploid() const
	    {
	    const_cast<Genotype*>(this)->val2 = HAPLOID_VAL;
	    }

	void setVals( AlleleType a1, AlleleType a2 )
	    {
	    val1 = a1;
	    val2 = a2;
	    }

	/// Normally only needs to be called at parse time.
	/// Canonicalize unordered pair so they can be easily compared.
	void condense()
	    {
	    if ( val1 == MISSING_VAL )
		val1 = val2;
	    }

	// NOTE *3*: if uses a constructor can't be used in a union
	// (Token) nor memcpy()d (GenotypeParser::Individual):
	#if 0
	    Genotype() {}
	    Genotype( AlleleType a1, AlleleType a2 ) : val1(a1), val2(a2) {}
	#endif

    } __attribute__ ((packed));




//-----------------------------------------------------------------------------
/// Array of Genotypes
//-----------------------------------------------------------------------------

typedef std::vector<Genotype> GenotypeArray;



//-----------------------------------------------------------------------------
/// A Haplotype is a phased genotype, so we reuse the same class but consider
/// the order of the alleles to be significant.
///
/// The first element was received from the father, the second from the mother.
//-----------------------------------------------------------------------------

class Haplotype : public Genotype
    {
    public:
	Haplotype() {}
	Haplotype( const AlleleType & a1, const AlleleType & a2 )
	    {
	    setVals( a1, a2 );
	    }

	AlleleType getPaternalAllele() const { return getVal1(); }
	AlleleType getMaternalAllele() const { return getVal2(); }
	AlleleType paternal() const { return getVal1(); } ///< Shorthand
	AlleleType maternal() const { return getVal2(); } ///< Shorthand
    };



} // ---- end namespace genepi



/** @} */



#endif // ! __base_Genotype_h
