#include "pvector.h"

//pvector::pvector()
//{
//}
//
//pvector::~pvector()
//{
//}


template <class T>
bool pvector<T>::is_normalized()
{
  // TODO: Accumulate only when necessary, i.e.
  // vector was changed.
  sum = accumulate(this->begin(), this->end(), 0.0);
  return (fabs(sum - PVECTOR_SUM) < PVECTOR_PRECISION);
}

template <class T>
void pvector<T>::normalize()
{
  if (!this->is_normalized()) { // Updates sum
    divisor.setDivisor(sum);
    transform(this->begin(), this->end(), this->begin(), divisor);
  }
}

template <class T>
void pvector<T>::verify()
{
  if (!this->is_normalized()) {
    throw string("pvector<T>::verify(): doesn't sum up to one.");
  }
  ZeroOrMore<T> z = ZeroOrMore<T>();
  if((unsigned)count_if(this->begin(), this->end(), z) != this->size()) {
    throw string("pvector<T>::verify(): Some elements are less than zero.");
  }
}

template <class T>
void pvector<T>::snapToZero()
{
  transform(this->begin(), this->end(), this->begin(), snapper);
  this->normalize();
}

// Providing explicit template instantiation to avoid linker errors.
// It needs to be _after_ function definitions.
template bool pvector<double>::is_normalized();
template void pvector<double>::normalize();
template void pvector<double>::verify();
template void pvector<double>::snapToZero();
