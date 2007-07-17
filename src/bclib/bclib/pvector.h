// *-*-C++-*-*
#ifndef PVECTOR_H_
#define PVECTOR_H_

#include "bclib/Array.hh"

BEGIN_BCLIB_NAMESPACE

#define PVECTOR_THRESHOLD 1e-20

/// Extension of STL Vector to handle vectors of probabilities
template<class T>
class pvector : public Array<T>{

public:
  pvector<T>() : Array<T>(){
    threshold = PVECTOR_THRESHOLD;
  }
  pvector<T>(unsigned n) : Array<T>(n){};
  void normalize();
  bool verify();
  bool is_normalized();
  void snapToZero();
  void snapToZero(const T t_threshold);
  //void print(std::ostream&, const T precision, const char* sep)const;

private:
  T sum;
  T threshold;


};

END_BCLIB_NAMESPACE
#endif /*PVECTOR_H_*/
