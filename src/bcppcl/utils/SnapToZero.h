// *-*-C++-*-*
#ifndef SNAPTOZERO_H_
#define SNAPTOZERO_H_

#include "../bcppcl.h"

BEGIN_BCPPCL_NAMESPACE

template <class T>
void SnapToZero(T& x, const T& threshold = PVECTOR_THRESHOLD){
    if (x < threshold)  x = 0;
}

template <class T>
bool IsNegative(T x){
  return (x < 0);
}

END_BCPPCL_NAMESPACE

#endif /*SNAPTOZERO_H_*/
