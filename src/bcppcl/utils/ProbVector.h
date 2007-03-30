// *-*-C++-*-*
#ifndef PVECTOR_H_
#define PVECTOR_H_


#include <vector>

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

#endif
