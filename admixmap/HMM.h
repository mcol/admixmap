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
  Vector_i Sample();
  void UpdateParameters( Matrix_d&, MatrixArray_d&, MatrixArray_d&, bool );
   Vector_d GetStateProbs( int );
   double getLikelihood();
   
private:
  void CheckArguments();
  void UpdateFwrdBckwdProbabilities();

  int States; //number of states of Markov chain, m in book
  //There are k*k (=D in Chromosome)states since k populations and 2 chromosomes

  int Transitions; //length of chain
  // = # composite Loci, (=L in Chromosome)

   bool _CalculateBeta;
   double factor;

  Matrix_d StationaryDist;
  MatrixArray_d TransitionProbs;
  //(L-1) * D * D array
  MatrixArray_d Likelihood;
  MatrixArray_d alpha;//forward and
  MatrixArray_d beta;//backward proabilities

  //alpha and beta are arrays of length L corresponding to {\alpha_t} and {\beta_t}
  //each beta(t) is a m * 1 matrix
  //each alpha is a ?? 
};

#endif /* ! HMM_H */
