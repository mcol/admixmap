// *-*-C++-*-*
#ifndef ISNEGATIVE_H_
#define ISNEGATIVE_H_

#include "bcppcl/bcppcl.h"

BEGIN_BCPPCL_NAMESPACE


template <class T>
bool IsNegative(T x){
  return (x < 0);
}

END_BCPPCL_NAMESPACE

#endif
