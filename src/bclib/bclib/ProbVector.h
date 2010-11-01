// *-*-C++-*-*
#ifndef PVECTOR_H_
#define PVECTOR_H_

#include "bclib/bclib.h"
#include <ostream>
#include <vector>

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */


/// Extension of STL Vector to handle vectors of probabilities
class ProbVector : public std::vector<double>{
public:
  ProbVector(){};
  ProbVector(unsigned n) : std::vector<double> (n){};
  void normalize();
  bool verify();
  bool is_normalized();
  void snapToZero();
  void print(std::ostream& os);

private:
  double sum;
};

/** @} */

END_BCLIB_NAMESPACE

#endif
