#ifndef PMBZARRAY_H
#define PMBZARRAY_H

#include "bclib/bclib.h"
#include "bclib/cvector.h"
#include <iostream>
#include <iterator>
#include <stdexcept>


BEGIN_BCLIB_NAMESPACE


//-----------------------------------------------------------------------------
/// Base class for arrays, same as std::vector but with extra functions, tuned
/// towards arithmetic operations.
//-----------------------------------------------------------------------------

template <class T>
class Array : public genepi::cvector<T>{

  typedef genepi::cvector<T> SUPER;

public:
  Array() : SUPER() {}

  Array(unsigned n) : SUPER (n){}
  virtual ~Array(){};



  //-------------------------------------------------------------------------
  // Due to C++ templates' mind-numbing complexity, these must be "re-declared"
  // here, rather than simply inheriting from std::vector.
  //-------------------------------------------------------------------------
  size_t size() const { return SUPER::size(); }
  typedef typename SUPER::iterator iterator;
  typedef typename SUPER::const_iterator const_iterator;
  const_iterator begin() const { return SUPER::begin(); }
  const_iterator end  () const { return SUPER::end  (); }
  iterator	 begin()       { return SUPER::begin(); }
  iterator	 end  ()       { return SUPER::end  (); }
  //-------------------------------------------------------------------------



  ///write to output stream, with optional separator
  virtual void print(std::ostream& os, const char* sep= " ")const{
    copy(this->begin(), this->end(), std::ostream_iterator<T>(os, sep));
    //os << std::endl;
  }

  ///scalar multiplication
  virtual void operator*=(const T& t){
    typename SUPER::iterator i = this->begin();
    for(;  i != this->end(); ++ i)
      *i *= t;
  }

  ///scalar multiplication
  Array operator*(const T& t)const{
    Array A ;
    typename SUPER::const_iterator i = this->begin();
    for(;  i != this->end(); ++ i)
      A.push_back( (*i) * t);
    return A;
  }


  /// Element-wise addition (RHS is any iterator that can produce an element
  /// type that makes arithmetic sense for "+=").
  ///
  /// The RHS iterator must have at least (and presumably should have exactly)
  /// the number of elements in the LHS array; no range checking is provided,
  /// however, beyond whatever the \a rhs iterator itself performs.

  template< typename T_Iter > Array & operator+=( T_Iter rhs )
    {
    for ( iterator it = begin() ; it != end() ; )
	*it++ += *rhs++;
    return *this;
    }

  /// Element-wise addition (RHS is an same-type array; it must have exactly the
  /// same number of elements as the LHS array).
  Array & operator+=( const Array & rhs )
    {
    if ( size() != rhs.size() )
	throw std::out_of_range( "+= array sizes don't match" );
    return operator+=(rhs.begin());
    }


  /// Apply a function to each element of the right-hand-side, then add to the corresponding element of the left-hand-side.
  template < typename T_Func > void addAfterFunction( const Array & rhs, T_Func func )
    {
    if ( size() != rhs.size() )
	throw std::out_of_range( "+= array sizes don't match" );

    for ( iterator lhs_it = begin(), rhs_it = rhs.begin() ; lhs_it != end() ; )
	*lhs_it++ += func(*rhs++);
    }

};


///class to wrap a 2d array as a 1d one, a kind of poor-man's 2d Blitz array. Useful only for rectangular arrays
template <class T>
class array2D : public Array<T>{
  typedef Array<T> SUPER;
public:
  ///default c'tor
  array2D(){
    mdim2 = 1;
  }
  ///constructor with initialisation
  array2D(unsigned dim1, unsigned dim2, T init = 0){
    assign(dim1, dim2, init);
  }
  ///assign dimensions and optional initial value
  void assign(unsigned dim1, unsigned dim2, T init = 0){
    mdim2 = dim2;
    SUPER::assign(dim1*dim2, init);
  }
  ///for write access
  T& operator()(unsigned i, unsigned j){
    return (*this)[i*mdim2 +j];
  }
  ///for read-only access
  const T& operator()(unsigned i, unsigned j)const{
    return (*this)[i*mdim2 +j];
  }

  ///write array to an output stream with one row per line
  void print(std::ostream& os)const{
    for(unsigned i = 0; i < this->size() / mdim2; ++i){
      copy(this->begin()+i*mdim2, this->begin()+ (i+1)*mdim2, std::ostream_iterator<T>(os, " "));
      os << std::endl;
    }
    os << std::endl;
  }

private:
  unsigned mdim2;
};

///class to wrap a 3d array as a 1d one, a kind of poor-man's 3d Blitz array. Useful only for rectangular arrays
template <class T>
class array3D : public Array<T>{
  typedef Array<T> SUPER;
public:
  ///default c'tor
  array3D(){
    mdim2 = 1;
    mdim3 = 1;
  }
  ///constructor with initialisation
  array3D(unsigned dim1, unsigned dim2, unsigned dim3, T init = 0){
    assign(dim1, dim2, dim3, init);
  }
  ///assign dimensions and optional initial value
  void assign(unsigned dim1, unsigned dim2, unsigned dim3, T init = 0){
    mdim2 = dim2;
    mdim3 = dim3;
    SUPER::assign(dim1*dim2*dim3, init);
  }

  ///for write access
  T& operator()(unsigned i, unsigned j, unsigned k){
    return (*this)[i*mdim2*mdim3 +j*mdim3 +k];
  }
  ///for read-only access
  const T& operator()(unsigned i, unsigned j, unsigned k)const{
    return (*this)[i*mdim2*mdim3 +j*mdim3 +k];
  }

  ///write array to an output stream with one 'row'(dim2) per line
  void print(std::ostream& os)const{
    for(unsigned i = 0; i < this->size() / (mdim2*mdim3); ++i){
      os << "[" << i+1 << ",,]" << std::endl;
      for(unsigned j = 0; j < mdim2; ++j){
	copy(this->begin()+i*mdim2*mdim3 + j*mdim3, this->begin()+ +i*mdim2*mdim3 + (j+1)*mdim3, std::ostream_iterator<T>(os, " "));
	os << std::endl;
      }
      os << std::endl;
    }
    os << std::endl;
  }

private:
  unsigned mdim2, mdim3;
};

END_BCLIB_NAMESPACE

#endif
