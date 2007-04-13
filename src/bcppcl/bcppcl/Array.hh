#ifndef PMBZARRAY_H
#define PMBZARRAY_H

#include "bcppcl/bcppcl.h"
#include <vector>
#include <iostream>
#include <iterator>


BEGIN_BCPPCL_NAMESPACE

//base class for arrays, same as std::vector but with extra functions
template <class T>
class Array : public std::vector<T>{
public:
  Array() : std::vector<T>() {}

  Array(unsigned n) : std::vector<T> (n){}
  virtual ~Array(){};

  ///write to output stream, with optional separator
  virtual void print(std::ostream& os, const char* sep= " ")const{
    copy(this->begin(), this->end(), std::ostream_iterator<T>(os, sep));
    //os << std::endl;
  }

  ///scalar multiplication
  virtual void operator*=(const T& t){
    typename std::vector<T>::iterator i = this->begin();
    for(;  i != this->end(); ++ i)
      *i *= t;
  }
  ///scalar multiplication
  Array operator*(const T& t)const{
    Array A ;
    typename std::vector<T>::const_iterator i = this->begin();
    for(;  i != this->end(); ++ i)
      A.push_back( (*i) * t);
    return A;
  }



};


///class to wrap a 2d array as a 1d one, a kind of poor-man's 2d Blitz array. Useful only for rectangular arrays
template <class T>
class array2D : public Array<T>{
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
    ((std::vector<T>*)this)->assign(dim1*dim2, init);
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
    ((std::vector<T>*)this)->assign(dim1*dim2*dim3, init);
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

END_BCPPCL_NAMESPACE

#endif
