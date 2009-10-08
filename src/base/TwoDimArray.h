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
/// \file TwoDimArray.h
/// TwoDimArray template.
//=============================================================================

#ifndef __base_TwoDimArray_h
#define __base_TwoDimArray_h


#include <cstring>	// memcpy(), size_t

#include "bclib/exceptions.h"



namespace genepi { // ----

/** \addtogroup base
 * @{ */




//-----------------------------------------------------------------------------
//
/// Two-dimensional array of fixed dimensions which are known at the time of construction.
///
/// A two-dimensional array which is stored in a "flattened" one-dimensional
/// array, allocated on the heap.  The sizes of the dimensions are fixed and
/// constant, i.e. each of the "rows" must have the same number of elements, and
/// the sizes of each dimension cannot be changed after the time of
/// construction.  Neither of these constraints are present in the
/// <CODE>std::vector<std::vector<T> ></CODE> approach, but the "flattened"
/// approach may be advantageous in performance-critical applications,
/// especially for objects that are allocated and freed often.
///
/// The data-element type is template parameter @a T, as are the data-types of
/// the two indexes (@a I1 and @a I2), although these default to size_t.  We
/// require that both dimensions must be of non-zero size.
///
/// Indexed access is via get(), at(), and set().  Range-checking is implemented
/// for all access, but can be compiled out via <CODE>-DAGGRESSIVE_RANGE_CHECK=0</CODE>.
/// There is no operator[](), it's not straightforward to overload for multiple indexes.
///
/// There is also an Iterator, which allows iteration over all elements of the
/// array with no guarantee as to order; and some support for applying functors
/// to the elements via applyToAll().  The iterator knows the end of its range,
/// i.e. is not the "STL" style iterator which requires a test for range-end but
/// rather a "Java" style iterator with a hasNext() method.
//
//-----------------------------------------------------------------------------

template < typename T, typename I1 = size_t, typename I2 = size_t > class TwoDimArray
    {

    private:
	const I1 i1size ;
	const I2 i2size ;
	T *	 storage;


    protected:

	size_t idxOf( const I1 & i1, const I2 & i2 ) const
	    {
	    gp_assert_lt( i1, i1size );
	    gp_assert_lt( i2, i2size );
	    return (i2 * i1size) + i1;
	    }

	size_t n_els() const { return (i1size * i2size); }



    public:

	//------------------------------------------------------------------
	// Standard constructor
	//------------------------------------------------------------------
	TwoDimArray( I1 _i1size, I2 _i2size ) :
		i1size ( _i1size ) ,
		i2size ( _i2size )
	    {
	    storage = new T [ n_els() ];
	    }


	//------------------------------------------------------------------
	// Copy constructor
	//------------------------------------------------------------------
	TwoDimArray( const TwoDimArray & rhs ) :
		i1size ( rhs.i1size ) ,
		i2size ( rhs.i2size )
	    {
	    gp_assert_gt( i1size, 0 );
	    gp_assert_gt( i2size, 0 );

	    storage = new T [ n_els() ];
	    memcpy( storage, rhs.storage, n_els() * sizeof(*storage) );
	    }


	//------------------------------------------------------------------
	// Destructor
	//------------------------------------------------------------------
	~TwoDimArray()
	    {
	    delete[] storage;
	    }


	//------------------------------------------------------------------
	// Assignment operator
	//------------------------------------------------------------------
	TwoDimArray & operator=( const TwoDimArray & rhs )
	    {
	    i1size = rhs.i1size;
	    i2size = rhs.i2size;
	    delete[] storage;
	    storage = new T [ n_els() ];
	    memcpy( storage, rhs.storage, n_els() * sizeof(*storage) );
	    }



	//------------------------------------------------------------------
	// Element access
	//------------------------------------------------------------------

	T &	  get( I1 i1, I2 i2 )	    { return storage[ idxOf( i1, i2 ) ]; } ///< Element access (as non-const reference).
	const T & get( I1 i1, I2 i2 ) const { return storage[ idxOf( i1, i2 ) ]; } ///< Element access (const version).
	T &	  at ( I1 i1, I2 i2 )	    { return storage[ idxOf( i1, i2 ) ]; } ///< Element access.  Synonym for get().

	/// Set an element.  This can also be accomplished via <CODE>at(i1,i2) = val;</CODE>
	void set( const I1 & i1, const I2 & i2, const T & val ) { at(i1,i2) = val; }


	I1 getI1size() const { return i1size; } ///< Get the first dimension size ("number of rows").
	I1 getI2size() const { return i2size; } ///< Get the second dimension size ("number of columns").



	//------------------------------------------------------------------
	//
	/// This is perhaps a little dodgy, but we'll offer iterative access to
	/// the whole contents; but no particular ordering of the elements is
	/// guaranteed.  Since the dimensions are required to be of non-zero
	/// size, there is always at least one element in the array, so
	/// hasNext() does not need to be checked initially; rather, after
	/// construction, operator*() can be immediately accessed.
	///
	/// For example,
	/// <PRE>
	/// typedef TwoDimArray<size_t,size_t,double> DArray;
	/// void print( const DArray & a )
	///	{
	///	for ( DArray::ConstIterator it( a ) ; it.hasNext() ; it.advance() )
	///	    cout << *it << endl;
	///	}
	/// </PRE>
	//
	//------------------------------------------------------------------

	class Iterator
	    {
	    private:
		TwoDimArray & array;
		size_t	      idx;
	    public:
		Iterator( TwoDimArray & _array ) : array( _array ), idx( _array.n_els() - 1 ) {}
		bool hasNext() const { return (idx != 0); }
		void advance() { gp_assert(hasNext()); --idx; }
		Iterator & operator++() { advance(); return *this; }
		T & operator*() { return array.storage[ idx ]; }
	    };

	class ConstIterator
	    {
	    private:
		const TwoDimArray & array;
		size_t		    idx;
	    public:
		ConstIterator( const TwoDimArray & _array ) : array( _array ), idx( _array.n_els() - 1 ) {}
		bool hasNext() const { return (idx != 0); }
		void advance() { gp_assert(hasNext()); --idx; }
		ConstIterator & operator++() { advance(); return *this; }
		const T & operator*() { return array.storage[ idx ]; }
	    };



	//------------------------------------------------------------------
	/// The ability to apply a functor to every element.
	//------------------------------------------------------------------
	template< typename Functor > void applyToAll( Functor ftor )
	    {
	    for ( ConstIterator it( *this ) ; it.hasNext() ; it.advance() )
		ftor( *it );
	    };



	//------------------------------------------------------------------
	// A few simple aggregate functions:
	//------------------------------------------------------------------

	/// Sum of all elements.
	T sum() const
	    {
	    T rv = 0;
	    for ( ConstIterator it( *this ) ; it.hasNext() ; it.advance() )
		rv += *it;
	    return rv;
	    }

	/// Assign value @a x to all elements.
	void assign( const T & x )
	    {
	    for ( Iterator it( *this ) ; it.hasNext() ; it.advance() )
		*it = x;
	    }
    };



/** @} */

} // ---- end namespace genepi



#endif // ! __base_TwoDimArray_h
