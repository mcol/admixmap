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
  void Sample(int *);  /* samples hidden states */
  void UpdateFwrdBckwdProbabilities( double *StationaryDist, double **Likelihood,  bool CalculateBeta);
  void SetTProb(int, int, int, double);
  void GetStateProbs(double *,  int );
  double getLikelihood();
  
private:
  
  int States; //number of states of Markov chain, m in book
  //There are k*k (=D in Chromosome)states since k populations and 2 chromosomes
  
  int Transitions; //length of chain
  // = # composite Loci, (=L in Chromosome)
  double ***TransitionProbs;//L * D * D array of transition probabilities 
 
  double factor;//used for adjustment of alpha to avoid underflow
  
  double ** alpha;//forward and
  double **beta;//backward probabilities
  //alpha and beta are L x m arrays 
  double *betarow;//used to calculate backward probs;
};

#endif /* ! HMM_H */
