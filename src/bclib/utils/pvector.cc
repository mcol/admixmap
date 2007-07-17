#include "bclib/pvector.h"
#include <cmath>//for fabs
#include <algorithm>
#include <numeric>
//#include <string>//for exceptions
#include "IsNegative.h" //for verify
#include <iterator>

#define PVECTOR_PRECISION 1e-10
#define PVECTOR_SUM 1.0

BEGIN_BCLIB_NAMESPACE

template <class T>
bool pvector<T>::is_normalized(){
  // TODO: Accumulate only when necessary, i.e.
  // when vector is changed.
  sum = accumulate(this->begin(), this->end(), 0.0);
  return (fabs(sum - PVECTOR_SUM) < PVECTOR_PRECISION);
}

template <class T>
void pvector<T>::normalize(){
  if (!this->is_normalized()) { // Updates sum
    (*this) *= 1.0 / sum; 
  }
}

template <class T>
bool pvector<T>::verify(){
  if (!this->is_normalized()) {
    throw ("pvector<T>::verify(): doesn't sum to one.");
  }
  if((unsigned)count_if(this->begin(), this->end(), IsNegative<T>) > 1) {
    throw ("pvector<T>::verify(): Some elements are less than zero.");
  }
  return true;
}

template <class T>
void pvector<T>::snapToZero(){
 for(typename pvector<T>::iterator i = this->begin(); i != this->end(); ++i){
  if(*i < threshold) *i = 0;
 } 
 snapToZero(threshold);
 this->normalize();
}

template <class T>
void pvector<T>::snapToZero(const T t_threshold){
  threshold = t_threshold;
  snapToZero();
}

// Providing explicit template instantiation to avoid linker errors.
template bool pvector<double>::is_normalized();
template void pvector<double>::normalize();
template bool pvector<double>::verify();
template void pvector<double>::snapToZero();
template void pvector<double>::snapToZero(double);

template bool pvector<float>::is_normalized();
template void pvector<float>::normalize();
template bool pvector<float>::verify();
template void pvector<float>::snapToZero();
template void pvector<float>::snapToZero(float);

END_BCLIB_NAMESPACE
