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
/// \file AlleleArray.h
/// Interface of AlleleArray class.
//=============================================================================

#ifndef __base_AlleleArray_h
#define __base_AlleleArray_h



#include "SimpleLocusArray.h"
#include "Genotype.h"		// AlleleType
#include "Organism.h"		// PopIdx
#include "exceptions.h"

#include <vector>



/// Should the array locally cache the number of populations (K)?
#define AA_K_REF	1

/// Should we compile ostream support for Allele probability tables (typically
/// used for debugging)?
#define AA_OSTREAM	1



#if AA_OSTREAM
    #include <iosfwd>
#endif



namespace genepi { // ----



/** \addtogroup base
 * @{ */



//-----------------------------------------------------------------------------
//
/// A double-indexed array on allele value (max number of alleles for that
/// simple-locus), and ancestry-index (i.e. population-index).
/// Used to store allele frequencies and probabilities.
///
/// @warning
/// <SPAN STYLE="font-weight: bold; color: red;">IMPORTANT!</SPAN>: see
/// <A HREF="Pedigree_8cc.html#note-1"><B>NOTE *1*</B> in Pedigree.cc</A> regarding
/// copy constructors and assignment operators.
//
//-----------------------------------------------------------------------------

template < typename T /*,size_t K*/ > class AlleleArray
    {
    private:
	const SimpleLocus &	loc;
	#if AA_K_REF
	    PopIdx		K;
	#endif
	T *			data; // not "*const" for NOTE *1* in Pedigree.cc

    public:
	AlleleArray( const SimpleLocus & _loc, size_t _K ) :
	    loc	 ( _loc					) ,
	    #if AA_K_REF
		K( _K					) ,
	    #endif
	    data ( new T [ _loc.getNumAlleles() * _K ]	) {}

	AlleleArray( const AlleleArray & rhs );

	~AlleleArray() { delete[] data; }

	AlleleArray & operator=( const AlleleArray & rhs );

    public:
	T & at( Genotype::AlleleType al, PopIdx pop )
	    {
	    gp_assert_le( al, loc.getNumAlleles() );
	    #if AA_K_REF
		gp_assert_lt( pop, K );
	    #endif

	    --al; // AlleleType is "1-based", indexes are "0-based"

	    #if AA_K_REF
		return data[ (al * K) + pop ];
	    #else
		return data[ (pop * loc.getNumAlleles()) + al ];
	    #endif

	    //return data[ (al * K) + pop ]; // Fastest but must know K at compile-time
	    }

	/// Const version
	const T & at( Genotype::AlleleType al, PopIdx pop ) const
	    { return const_cast<AlleleArray*>(this)->at(al,pop); }

	const SimpleLocus & getLoc() const { return loc; }

	#if AA_K_REF
	    PopIdx getK() const { return K; }
	#endif
    };



//-----------------------------------------------------------------------------
/// Table of allele frequencies (normalized probabilities) at _one_ locus
//-----------------------------------------------------------------------------

class AlleleProbTable : public AlleleArray<double>
    {
    public:
	AlleleProbTable( const SimpleLocus & _loc, size_t _K ) :
		AlleleArray<double>(_loc,_K)
	    {
	    }

	/// Normalize the probability distribution for each population across
	/// alleles.  If any population contains no frequency for any allele,
	/// all values will be initialized to (1/N).
	void normalizeProbs();


	/// Useful for debugging: output an AlleleProbTable to an ostream
	#if AA_OSTREAM
	    void print( std::ostream & os, const std::vector<std::string> & pops ) const;
	#endif // AA_OSTREAM
    };



//-----------------------------------------------------------------------------
/// Table of allele frequencies (normalized probabilities) at each locus.
/// Vector of AlleleProbTable s, indexed on Simple-Locus index.
///
/// To-do: reimplement at private subclass, with range-checking enabled for
/// operator[]().  In the meantime, use at() instead.
//-----------------------------------------------------------------------------

typedef std::vector<AlleleProbTable> AlleleProbVect;



//-----------------------------------------------------------------------------
// Inlined template method implementations:
//-----------------------------------------------------------------------------

template < typename T > AlleleArray<T>::AlleleArray( const AlleleArray & rhs ) :
	loc	    ( rhs.loc	) ,
	#if AA_K_REF
	    K	    ( rhs.K	) ,
	#endif
	data	    ( rhs.data	)
    {
    // !!!WARNING!!! -- see NOTE *1* in Pedigree.cc
    const_cast<AlleleArray&>(rhs).data = 0;
    }


template < typename T > AlleleArray<T> & AlleleArray<T>::operator=( const AlleleArray & rhs )
    {
    gp_assert( &loc == &rhs.loc );
    #if AA_K_REF
	K = rhs.K;
    #endif
    data = rhs.data;

    // !!!WARNING!!! -- see NOTE *1* in Pedigree.cc
    const_cast<AlleleArray&>(rhs).data = 0;

    return *this;
    }



/** @} */



} // ---- end namespace genepi



#endif // ! __base_AlleleArray_h
