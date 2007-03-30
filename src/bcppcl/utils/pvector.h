// *-*-C++-*-*
#ifndef PVECTOR_H_
#define PVECTOR_H_

#include "Array.hh"

#define PVECTOR_THRESHOLD 1e-20

BEGIN_BCPPCL_NAMESPACE

/// Extension of STL Vector to handle vectors of probabilities
template<class T>
class pvector : public Array<T>{

public:
  pvector<T>() : Array<T>(){}
  pvector<T>(unsigned n) : Array<T>(n){};
  void normalize();
  bool verify();
  bool is_normalized();
  void snapToZero();
  void print(std::ostream&, const T precision, const char* sep)const;

private:
  T sum;

};

END_BCPPCL_NAMESPACE
#endif /*PVECTOR_H_*/
