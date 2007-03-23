#ifndef PVECTOR_H_
#define PVECTOR_H_

#include <cmath>
#include <vector>
#include <algorithm>
#include <string>
#include <numeric>

#include "DivideBy.h"
#include "ZeroOrMore.h"
#include "SnapToZero.h"

using std::string;
using std::vector;

#define PVECTOR_PRECISION 1e-10
#define PVECTOR_SUM 1.0
#define PVECTOR_SNAP 1e-20

/// Extension of STL Vector to handle vectors of probabilities
template<class T>
class pvector : public vector<T>
{
private:
  T sum;
  SnapToZero<T> snapper;
  DivideBy<T> divisor;
public:
  pvector<T>() {
    snapper.setThreshold(PVECTOR_SNAP);
  }
  void normalize();
  void verify();
  bool is_normalized();
  void snapToZero();
};

#endif /*PVECTOR_H_*/
