#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

// need an STL vector to hold x, h, gradient, z1, z2 for each x 
// then just loop over vector to evaluate and normalize cumulative probs

double ::AreaUnderTangentCurve(const double x1, const double h, const double g, 
			       const double z1, const double z2, 
			       const bool lower, const bool upper) {
  // returns unnormalized area under exponential curve formed by tangent to log density  
  // arguments: x1 abscissa, h height of log density at x1, g gradient at x1, 
  // z1 and z2 lower and upper bounds on x axis, 
  // lower=true if lower bound, upper=true if upper bound
  const double absg = 0.00000001;
  double unscaled = 0.0;
  // should check that z1 <= x1 if lower, and x1 <= z2 if upper
  // to prevent computational underflow, could return unscaled with h = log scale factor
  
  if(lower & upper & fabs(g) > absg) {
    unscaled = (exp(g*(z2 - x1)) - exp(g*(z1 - x1))) / g;
  } else if(lower & g < 0) {
    unscaled = exp(g*(z2 - x1)) / g; 
  } else if(upper & g > 0) {
    unscaled = exp(g*(z1 - x1)) / g;
  } else if(lower & upper & fabs(g) <= absg) { // gradient close to 0
    // evaluate difference between exponentials using first two terms of series 
    unscaled = (z2 - z1 + 0.5*(z2 + z1 - 2.0*x1)*(z2 - z1)*g);
  } else {
    cout << "Error in arguments passed to AreaUnderTangentCurve \n";
    exit(1);
  }
  return exp(h) * unscaled;
} 


int main() {
  // test algorithm
  double x1 = 3.0;
  double h =  -10.0;
  double z1 = -5.0;
  double z2 = 15.0;
  double gradient = 1.5;
  bool lower = false;
  bool upper=true;

  cout << AreaUnderTangentCurve(x1, h, gradient, z1, z2, lower, upper) << "\n";
 
}
