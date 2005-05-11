#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

double **alloc2D_d(int m,int n)
{
  double **M;
  M = new double*[m];
  for(int i=0; i<m; ++i) M[i]=new double[n];
  return M;
}


// argument oldProbs is square array of size K, K
// for forward recursions, pass alpha_t and multiply newProbs by observation probs 
// for backward recursions, pass array of products lambda_t+1[jj] * beta_t+1[jj] 
// updates oldProbs (scaled to sum to 1), newProbs, scaleFactor
void ::RecursionProbs(const int K, const double p, const double *f, 
		      const double stateArrivalProbs[][2], double oldProbs[][2], 
		      double & scaleFactor, double newProbs[][2]) {
  static double Sum = 0;
  double *rowProb = new double[K];
  double *colProb = new double[K];
  double *rowSum = new double[K];
  double *colSum = new double[K];
  double **cov;
  cov = alloc2D_d(K, K);

  for( int i = 0; i < K; i++) {
    rowProb[i] = 0;
    colProb[i] = 0;
    rowSum[i] = 0;
    colSum[i] = 0;
  }

  // scale array oldProbs so that elements sum to 1, and accumulate row and col sums
  for( int j0 = 0; j0 <  K; ++j0 ) { // leave out last row
    for( int j1 =0; j1 < K; ++j1 ) { // leave out last col
      Sum += oldProbs[j0][j1];
    }
  }
  scaleFactor = 1 / Sum;
  for( int j0 = 0; j0 <  K; ++j0 ) { 
    for( int j1 =0; j1 < K; ++j1 ) { 
      oldProbs[j0][j1] = oldProbs[j0][j1] * scaleFactor;
      rowProb[j0] += oldProbs[j0][j1];
      colProb[j1] += oldProbs[j0][j1];
    }
  }

  // calculate covariance of ancestry states as deviation from product of row and col prob
  for(int j0 = 0; j0 <  K-1; ++j0) { // leave out last row
    for(int j1 =0; j1 < K-1; ++j1) { // leave out last col
      cov[j0][j1] = p * ( oldProbs[j0][j1] - rowProb[j0] * colProb[j1] );
      rowSum[j0] += cov[j0][j1];
      colSum[j1] += cov[j0][j1];
    }
  }
  // calculate last row except for last col, by subtracting colSum from 0
  for( int j1 = 0; j1 < K-1; ++j1 ) {
    cov[K-1][j1] =  -p *  colSum[j1];
    rowSum[K-1] += cov[K-1][j1];
  }
  // calculate last col by subtracting rowSum from 0
  for( int j0 = 0; j0 < K; ++j0 ) {
    cov[j0][K-1] =  -p * rowSum[j0];
  }

  // calculate expectation of product as covariance plus product of expectations
  // can speed up with Fourier transform 
  for(int j0 = 0; j0 < K; ++j0) {
    for(int j1 =0; j1 < K; ++j1) {
      newProbs[j0][j1] = cov[j0][j1] + 
      ( f[0]*rowProb[j0] + stateArrivalProbs[j0][0] ) * ( f[1]*colProb[j1] + stateArrivalProbs[j1][1] );
    }
  }

}

int main() {
  // test algorithm
  const int K = 2;
  double Theta[K][2] = { {0.2, 0.4}, {0.8, 0.6} };
  double f[2] = {02, 0.2};
  double p = f[0]*f[1];
  double scaleFactor = 0;

  double oldProbs[K][K] = { {1, 9}, {3, 7} };
  double newProbs[K][K] = { {0, 0}, {0, 0} };
  double stateArrivalProbs[K][K] = { {0, 0}, {0, 0} };

  // calculate state arrival probs on each gamete
  // store and re-use for backwards recursion 
  for(int i = 0; i < K; ++i) {
    stateArrivalProbs[i][0] = (1 - f[0]) * Theta[i][0];
    stateArrivalProbs[i][1] = (1 - f[1]) * Theta[i][1];
  }

  RecursionProbs(K, p, f, stateArrivalProbs, oldProbs, scaleFactor, newProbs);
 
  cout << "\nscalefactor " << scaleFactor << "\n";
  for( int j0 = 0; j0 <  K; ++j0 ) { 
    for( int j1 =0; j1 < K; ++j1 ) { 
      cout << newProbs[j0][j1] << " ";
    }
    cout << "\n";
  }

}
