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
/// \file cvector.h
/// Definition of the cvector template.
//=============================================================================

#ifndef __bclib_cvector_h
#define __bclib_cvector_h


#include <vector>

#include <stdexcept>	// std::out_of_range

#include "config.h"	// AGGRESSIVE_RANGE_CHECK

#include "bclib/exceptions.h"



namespace genepi { // ----

/** \addtogroup bclib
 * @{ */



template<typename T> class divideBy
    {
    private:
	const T denominator;
    public:
	divideBy( const T & _denominator ) : denominator( _denominator ) {}
	template < typename U > void operator()( U & lhs ) const { lhs /= denominator; }
    };



template<typename T> class multiplyBy
    {
    private:
	const T factor;
    public:
	multiplyBy( const T & _factor ) : factor( _factor ) {}
	template < typename U > void operator()( U & lhs ) const { lhs *= factor; }
    };



//-----------------------------------------------------------------------------
// cvector template
//
/// A drop-in replacement for std::vector, but with (conditionally-compiled) range-checking.
///
/// While the std::vector template provides range-checking for the @a at()
/// method, it lacks any for @a operator[].  This template aims to be a
/// method-for-method identical replacement, but with range-checking implemented
/// for both at() and operator[], which can then be "compiled out" conditionally
/// based on the preprocessor macro AGGRESSIVE_RANGE_CHECK.
//
//-----------------------------------------------------------------------------

template < typename T, typename Alloc = std::allocator<T> > class cvector
    {
    private:
	std::vector<T,Alloc> v;


    protected:
	#if AGGRESSIVE_RANGE_CHECK
	    void range_check( size_t el ) const
		{
		if ( el >= size() )
		    throw std::out_of_range( estr("cvector-index ") + el +
			" out of range (" + size() + ')' );
		}
	#else
	    void range_check( size_t /*el*/ ) const {}
	#endif

    public:

	// Constructors:
	cvector() {}
	explicit cvector( size_t init_size ) : v( init_size ) {}
	explicit cvector( size_t init_size, const T & init_val ) : v( init_size, init_val ) {}
	cvector( const cvector & rhs ) : v( rhs.v ) {}
	explicit cvector( const std::vector<T,Alloc> & rhs ) : v( rhs ) {}


	// Assignment:
	cvector & operator=( const cvector & rhs ) { v = rhs.v; return *this; }
	cvector & operator=( const std::vector<T,Alloc> & rhs ) { v = rhs; return *this; }

	// Access:
	size_t size() const { return v.size(); }

	const T & operator[]( size_t el ) const
	    {
	    range_check( el );
	    return v[el];
	    }

	T & operator[]( size_t el )
	    {
	    range_check( el );
	    return v[el];
	    }

	T &	  at( size_t el )	{ return operator[]( el ); }
	const T & at( size_t el ) const { return operator[]( el ); }

	T &	  at_unsafe( size_t el )       { return v[el]; }
	const T & at_unsafe( size_t el ) const { return v[el]; }

	void push_back ( const T & el ) { v.push_back( el ); }
	void pop_back  ( const T & el ) { v.pop_back ( el ); }


	// Iterators:
	typedef typename std::vector<T,Alloc>::iterator	      iterator;
	typedef typename std::vector<T,Alloc>::const_iterator const_iterator;

	iterator	begin()	      { return v.begin(); }
	const_iterator	begin() const { return v.begin(); }
	iterator	end  ()	      { return v.end(); }
	const_iterator	end  () const { return v.end(); }

	T &		back ()	      { range_check(0); return v.back (); }
	const T &	back () const { range_check(0); return v.back (); }
	T &		front()	      { range_check(0); return v.front(); }
	const T &	front() const { range_check(0); return v.front(); }
	T &		last ()	      { return back(); } ///< synonym for back()
	const T &	last () const { return back(); } ///< synonym for back()

	void reserve( size_t expected_size )	  { v.reserve( expected_size ); }
	void assign( size_t size, const T & el )  { v.assign( size, el ); }
	void resize( size_t size )		  { v.resize( size ); }
	void resize( size_t size, const T & el )  { v.resize( size, el ); }

	bool empty  () const { return v.empty(); }
	bool isEmpty() const { return empty(); } ///< Not in STL: synonym for empty()

	void clear() { v.clear(); }

	iterator erase( iterator position )		{ v.erase( position ); }
	iterator erase( iterator first, iterator last ) { v.erase( first, last ); }

	void erase( size_t position ) { v.erase( v.begin() + position ); } ///< Not in STL

	// --- Functor generic algorithms: ---

	template < typename Functor > void for_each( Functor func )
	    {
	    for ( iterator it = begin() ; it != end() ; ++it )
		func( *it );
	    }

	template< typename U > cvector & operator/=( const U & rhs ) { for_each( divideBy  <U>( rhs ) ); return *this; }
	template< typename U > cvector & operator*=( const U & rhs ) { for_each( multiplyBy<U>( rhs ) ); return *this; }

	template< typename U > cvector operator*( const U & rhs )
	    {
	    cvector rv( size() );
	    for ( size_t idx = size() ; idx-- != 0 ; )
		rv[idx] = at(idx) * rhs;
	    return rv;
	    }


	/// Direct access to the STL std::vector -- avoid if possible.
	const std::vector<T,Alloc> & getVector_unsafe() const { return v; }
	std::vector<T,Alloc> &	     getVector_unsafe()	      { return v; }

	/// Direct access to the data as a "C array" -- avoid if possible.
	const T * data_unsafe() const { return v.data(); }
	T *	  data_unsafe()	      { return v.data(); }

    };



//-----------------------------------------------------------------------------
// Some special care need be taken for STL's "bool" specialization to avoid
// compiler errors.
//-----------------------------------------------------------------------------
#if 0
    template<> inline const bool & cvector<bool>::operator[]( size_t el ) const
	{
	range_check( el );
	return v[el];
	}

    template bool & cvector<bool>::operator[]( size_t el )
	{
	range_check( el );
	return v[el];
	}
#endif



/** @} */

} // ---- end namespace genepi



#endif // ! __bclib_cvector_h
