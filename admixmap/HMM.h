// *-*-C++-*-*
//RevMCMC.h

/* Class to implement Hidden Markov Models (see MacDonald and Zucchini)
   Instantiated in class Chromosome
*/


#ifndef HMM_H
#define HMM_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <values.h>

#define FILENAMELENGTH 100

class Vector_i;
class Matrix_d;
class MatrixArray_d;

class HMM
{
public:
  HMM();
  HMM( int, int );
  ~HMM();
  void SetDimensions( int, int );
  // functions Sample, CheckArguments and UpdateFrwrdBackwdProbabilities have been edited to take
  // pointer to  transition matrices as argument so that HMM doesn't have to store a copy of this array
  Vector_i Sample(MatrixArray_d&);  /* samples hidden states */
  void UpdateParameters( Matrix_d&, MatrixArray_d&, MatrixArray_d&, bool );
  Vector_d GetStateProbs( int );
  double getLikelihood();
  
private:
  void CheckArguments(MatrixArray_d&);
  
  // void UpdateFwrdBckwdProbabilities();
  void UpdateFwrdBckwdProbabilities(MatrixArray_d&);
  
  int States; //number of states of Markov chain, m in book
  //There are k*k (=D in Chromosome)states since k populations and 2 chromosomes
  
  int Transitions; //length of chain
  // = # composite Loci, (=L in Chromosome)
  
  bool _CalculateBeta;
  double factor;
  
  Matrix_d StationaryDist;
  // as transition probs are calculated and stored in Chromosome, no need to 
  // store another copy in HMM object
  // MatrixArray_d TransitionProbs;
  
  //(L-1) * D * D array
  MatrixArray_d Likelihood;
  MatrixArray_d alpha;//forward and
  MatrixArray_d beta;//backward probabilities
  
  //alpha and beta are arrays of length L corresponding to {\alpha_t} and {\beta_t}
  //each beta(t) is a m * 1 matrix
  //each alpha(t) is an m * 1 matrix 
};

#endif /* ! HMM_H */
