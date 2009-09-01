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
  pvector<T>() : SUPER(){
    threshold = PVECTOR_THRESHOLD;
  }
  pvector<T>(unsigned n) : SUPER(n){};
  void normalize();
  bool verify();
  bool is_normalized();
  void snapToZero();
  void snapToZero(const T t_threshold);
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

  template< typename DestIter			    > void inv_softmax( DestIter dest ) const;
  template< typename DestIter, typename QualFunctor > void inv_softmax( DestIter dest, QualFunctor qualifies ) const;

  void inv_softmax( pvector & dest ) const
	{
	inv_softmax( dest.begin() );
	}
  template< typename QualFunctor > void inv_softmax( pvector & dest, QualFunctor qualifies ) const
	{
	inv_softmax( dest.begin(), qualifies );
	}


  // Convenience: inv_softmax_gt0(x) == inv_softmax(x,greater_than_0)
  template< typename DestIter > void inv_softmax_gt0( DestIter dest ) const
	{
	inv_softmax( dest, greater_than_0 );
	}
  void inv_softmax_gt0( pvector & dest ) const
	{
	inv_softmax( dest.begin(), greater_than_0 );
	}


  template< typename DestIter, typename QualFunctor > void softmax( DestIter dest, QualFunctor qualifies ) const;
  template< typename QualFunctor > void softmax( pvector & dest, QualFunctor qualifies ) const
	{ softmax( dest.begin(), qualifies ); }

  template< typename DestIter > void softmax_gt0( DestIter dest ) const
	{
	softmax( dest, greater_than_0 );
	}
  void softmax_gt0( pvector & dest ) const
	{
	softmax( dest.begin(), greater_than_0 );
	}


  template< typename U > pvector & operator/=( const U & rhs )
	    {
	    for ( iterator it = begin() ; it != end() ; ++it )
		*it /= rhs;
	    return *this;
	    }
  template< typename U > pvector & operator*=( const U & rhs )
	    {
	    for ( iterator it = begin() ; it != end() ; ++it )
		*it *= rhs;
	    return *this;
	    }


private:
  T sum;
  T threshold;


};


/** @} */

END_BCLIB_NAMESPACE

#endif /*PVECTOR_H_*/
