#include "DivideBy.h"

template <class T>
DivideBy<T>::DivideBy()
{
}

template <class T>
DivideBy<T>::DivideBy(const T& d)
: divisor(d)
{
  reciprocal = 1 / divisor;
}

template <class T>
void DivideBy<T>::setDivisor(const T& d)
{
  divisor = d;
  reciprocal = 1 / divisor;
}

template <class T>
DivideBy<T>::~DivideBy()
{
}

// Explicit template instantiation
template DivideBy<double>::DivideBy(double const&);
template DivideBy<double>::DivideBy();
template const double DivideBy<double>::operator()(const double &) const;
template void DivideBy<double>::setDivisor(double const&);
