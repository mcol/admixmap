/*
 *  Class to represent an array of probabilities.
 *  Copyright (c) 2006-2007 David O'Donnell and Paul McKeigue
 *  Portions Copyright (C) 2009 David Favro
 *
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 *
 */


// *-*-C++-*-*
#ifndef PVECTOR_H_
#define PVECTOR_H_

#include "bclib/Array.hh"

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */


#define PVECTOR_THRESHOLD 1e-20


/// Extension of STL Vector to handle vectors of probabilities
template<class T>
class pvector : public Array<T>{

    typedef Array<T> SUPER;

public:


  //-----------------------------------------------------------------------------
  // Constructors:
  //-----------------------------------------------------------------------------

  pvector() :
	SUPER() ,
	threshold( PVECTOR_THRESHOLD )
    {
    }

  pvector( size_t n, const T & def_val ) :
	SUPER( n, def_val ) ,
	threshold( PVECTOR_THRESHOLD )
    {
    }

  pvector( size_t n ) :
	SUPER( n ) ,
	threshold( PVECTOR_THRESHOLD )
    {
    }


  void normalize();
  bool verify();
  bool is_normalized();
  void snapToZero();
  void snapToZero( const T & t_threshold )
    {
    threshold = t_threshold;
    snapToZero();
    }

  //void print(std::ostream&, const T precision, const char* sep)const;



  //-------------------------------------------------------------------------
  // Due to C++ templates' mind-numbing complexity, these must be "re-declared"
  // here, rather than simply inheriting:
  //-------------------------------------------------------------------------
  size_t size() const { return SUPER::size(); }
  typedef typename SUPER::iterator iterator;
  typedef typename SUPER::const_iterator const_iterator;
  const_iterator begin() const { return SUPER::begin(); }
  const_iterator end  () const { return SUPER::end  (); }
  iterator	 begin()       { return SUPER::begin(); }
  iterator	 end  ()       { return SUPER::end  (); }
  //-------------------------------------------------------------------------

  static bool greater_than_0( T x ) { return x > 0.0; }
  static bool not_equal_to_0( T x ) { return x != 0.0; }

  template< typename DestIter			    > void inv_softmax( DestIter dest ) const;
  template< typename DestIter, typename QualFunctor > void inv_softmax( DestIter dest, QualFunctor qualifies, const T & defVal = 0 ) const;
  ///< Only peforms the tranformation on those elements which return true from
  ///< @a qualifies; all other elements of the @a dest array are set to @a
  ///< defVal.

  /// DDF: problem: how to make these specialized versions get chosen in preference to above
  ///	template-methods for arguments derived from cvector, such as pvector?
  void inv_softmax( pvector & dest ) const
	{
	if ( dest.size() < size() )
	    dest.resize( size() );
	inv_softmax( dest.begin() );
	}

  /// Only peforms the tranformation on those elements which return true from @a
  /// qualifies; all other elements of the @a dest array are set to @a defVal.
  template< typename QualFunctor > void inv_softmax( pvector & dest, QualFunctor qualifies, const T & defVal = 0 ) const
	{
	if ( dest.size() < size() )
	    dest.resize( size() );
	inv_softmax( dest.begin(), qualifies, defVal );
	}


  /// Convenience: inv_softmax_gt0(x) == inv_softmax(x,greater_than_0)
  /// DDF: perhaps this should be _ne0 -- I don't understand its true function.
  template< typename DestType > void inv_softmax_gt0( DestType dest, const T & defVal = 0 ) const
	{
	inv_softmax( dest, not_equal_to_0, defVal );
	}

  /// We must define this rather than use the above template-method because it
  /// effectively makes a *copy* of dest rather than passing in by reference,
  /// which then throws away the results when finished executing.
  void inv_softmax_gt0( pvector & dest, const T & defVal = 0 ) const
	{
	inv_softmax( dest, not_equal_to_0, defVal );
	}


  template< typename DestIter, typename QualFunctor > void softmax( DestIter dest, QualFunctor qualifies, const T & defVal = 0 ) const;
  template< typename QualFunctor > void softmax( pvector & dest, QualFunctor qualifies, const T & defVal = 0 ) const
	{
	if ( dest.size() < size() )
	    dest.resize( size() );
	softmax( dest.begin(), qualifies, defVal );
	}

  template< typename DestType > void softmax_gt0( DestType dest, const T & defVal = 0 ) const
	{
	softmax( dest, not_equal_to_0, defVal );
	}
  /// We must define this rather than use the above template-method because it
  /// effectively makes a *copy* of dest rather than passing in by reference,
  /// which then throws away the results when finished executing.
  void softmax_gt0( pvector & dest, const T & defVal = 0 ) const
	{
	softmax( dest, not_equal_to_0, defVal );
	}



  //-----------------------------------------------------------------------------
  // Element-wise operators
  //-----------------------------------------------------------------------------

  /// Simple element-wise division operator.
  template< typename U > pvector & operator/=( const U & rhs )
	    {
	    for ( iterator it = begin() ; it != end() ; ++it )
		*it /= rhs;
	    return *this;
	    }

  /// Simple element-wise multiplication operator.
  template< typename U > pvector & operator*=( const U & rhs )
	    {
	    for ( iterator it = begin() ; it != end() ; ++it )
		*it *= rhs;
	    return *this;
	    }


private:
  T sum; ///< *NOT* maintained, this is undefined except in very limited scope
  T threshold;

};


/** @} */

END_BCLIB_NAMESPACE

#endif /*PVECTOR_H_*/
