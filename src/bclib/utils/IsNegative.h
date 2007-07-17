// *-*-C++-*-*
#ifndef ISNEGATIVE_H_
#define ISNEGATIVE_H_

#include "bclib/bclib.h"

BEGIN_BCLIB_NAMESPACE


template <class T>
bool IsNegative(T x){
  return (x < 0);
}

END_BCLIB_NAMESPACE

#endif
