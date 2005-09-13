#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

void ::RecursionProbs2(const double ff, const double f[2], 
			 double *stateArrivalProbs, double *oldProbs, double *newProbs) {
  double Sum = 0.0, scaleFactor = 1.0;
  double row0Prob;
  double col0Prob;
  double Expectation0;
  double Expectation1;
  double cov;
  double Product;
  double RealF1;
  double ImaginaryF1;

  // version for K = 2
  // scale array oldProbs so that elements sum to 1, and accumulate row and col sums
  for( int j = 0; j <  4; ++j ) {
    Sum += oldProbs[j];
  }
  scaleFactor = 1.0 / Sum;
  row0Prob = ( oldProbs[0] + oldProbs[2] ) * scaleFactor;
  col0Prob = ( oldProbs[1] + oldProbs[3] ) * scaleFactor;

  // calculate expectations of indicator variables for each ancestry state on each gamete
  Expectation0 = f[0]*row0Prob + stateArrivalProbs[0];
  Expectation1 = f[1]*col0Prob + stateArrivalProbs[1];
  Product = Expectation0 * Expectation1;
  
  // calculate covariance of ancestry states as ff * deviation from product of row and col probs
  cov = ff * ( oldProbs[0]*scaleFactor - row0Prob * col0Prob );

  // calculate expectation of product as covariance plus product of expectations
  // evaluates newProbs[j0*K + j1] = cov + Expectation0 * Expectation1;
  // uses discrete Fourier transform  
  //   P0 = Expectation0 * Expectation1 + cov;
  //   P1 = (1 - Expectation0) * Expectation1 - cov;
  //   P2 = Expectation0 * (1 - Expectation1) - cov;
  //   P3 = (1 - Expectation0) * (1 - Expectation1) + cov;
  // F0 = 0.5*(P0 + P1 + P2 + P3) = 0.5
  // F1 = 0.5*((P0 - P2) + i*(P1 - P3)) = Real + i*Imaginary
  // F2 = 0.5*(P0 - P1 + P2 - P3)
  // F3 = 0.5*((P0 - P2) - i*(P1 - P3)) 
  RealF1      = Product - 0.5 * Expectation0 + cov;
  ImaginaryF1 = Expectation1 - Product - 0.5 * (1 - Expectation0) - cov;
  
  newProbs[0] = 0.5*Expectation0 + RealF1;
  newProbs[1] = 0.5*(1.0 - Expectation0) + ImaginaryF1;
  newProbs[2] = 0.5*Expectation0 - RealF1;
  newProbs[3] = 0.5*(1.0 - Expectation0) - ImaginaryF1;
  //undo scaling 
  for(int j0 = 0; j0 < 2; ++j0) {
    for(int j1 =0; j1 < 2; ++j1) {
      newProbs[j0*K + j1] *= Sum;
    }
  }
}

int main() {
  // test algorithm
  const int K = 2;
  const int States = 4;
  double Theta[K][2] = { {0.25, 0.2}, {0.75, 0.8} };
  double f[2] = {0.3, 0.3};
  double ff = f[0]*f[1];

  double oldProbs[States] = { 1, 9, 3, 7 };
  double newProbs[States] = { 0, 0, 0, 0 };
  double stateArrivalProbs[States] = { 0, 0, 0, 0 };

  // calculate state arrival probs on each gamete
  // store and re-use for backwards recursion 
  for(int i = 0; i < K; ++i) {
    stateArrivalProbs[i*K] = (1 - f[0]) * Theta[i][0];
    stateArrivalProbs[i*K + 1] = (1 - f[1]) * Theta[i][1];
  }

  RecursionProbs2(ff, f, stateArrivalProbs, oldProbs, newProbs);
 
  cout << "\nNewProbs\n";
  for( int j0 = 0; j0 <  K; ++j0 ) { 
    for( int j1 =0; j1 < K; ++j1 ) { 
      cout << newProbs[j0*K + j1] << " ";
    }
    cout << "\n";
  }

}
