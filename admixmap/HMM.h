// *-*-C++-*-*

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

class HMM
{
public:
  HMM();
  HMM( int inTransitions, int pops, bool isdiploid );
  ~HMM();
  void SetDimensions( int inTransitions, int pops, bool isdiploid );
  //void Sample(int *);  /* samples hidden states */
  void NewSample(int *C, Matrix_d &Admixture, double *f[], int StartLocus,int Mcol,bool isdiploid);
  //void UpdateFwrdBckwdProbabilities( double *StationaryDist, double **Likelihood,  bool CalculateBeta);
  //void SetTProb(int, int, int, double);
  //void GetStateProbs(double *,  int );
  void NewGetStateProbs( double * probs, int t);
  //double getLikelihood();
  double NewgetLikelihood();
  void UpdateProbsDiploid(double *f[],int StartLocus, Matrix_d& Admixture, double ***lambda, int Mcol, bool CalculateBeta);
  void UpdateProbsHaploid(double *f[], int StartLocus,Matrix_d& Admixture, double ***lambda, bool CalculateBeta);
  
private:
  int K;
  int States; //number of states of Markov chain, m in book
  //There are k*k (=D in Chromosome)states since k populations and 2 chromosomes
  
  int Transitions; //length of chain
  // = # composite Loci, (=L in Chromosome)
  //double ***TransitionProbs;//L * D * D array of transition probabilities 

  double sumfactor; 
  double factor;//used for adjustment of alpha to avoid underflow
  
  //double **alpha;//forward and
  //double **beta;//backward probabilities
  //alpha and beta are L x m arrays 

  //forward and backward probabilities
  //L x K x K arrays
  double ***newalpha, ***newbeta;
  double *p,*q,*r,*s, **Q ,**R, ***S, **U, **V, **W;
};

#endif /* ! HMM_H */
