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
  /* samples hidden states */
  void Sample(int *C, Matrix_d &Admixture, double *f[], bool isdiploid);
  void GetStateProbs( double * probs, int t);
  void NewGetStateProbs( double * probs, int t, double ***StateArrivalProbs, double *f[]);

  double getLikelihood();
  void UpdateProbsDiploid(double ***StateArrivalProbs,double *f[], Matrix_d& Admixture, double ***lambda, int Mcol, bool CalculateBeta);
  void UpdateForwardProbsDiploid(double ***StateArrivalProbs, double *f[], Matrix_d& Admixture, double ***lambda, int Mcol);
  void UpdateRevAlphas(double ***StateArrivalProbs, double *f[], Matrix_d& Admixture, double ***lambda, int Mcol);

  void UpdateBackwardProbsDiploid(double ***StateArrivalProbs,double *f[], Matrix_d& Admixture, double ***lambda, int Mcol);

  void UpdateProbsHaploid(double *f[], Matrix_d& Admixture, double ***lambda, bool CalculateBeta);

  void RecursionProbs(const double ff, const double f[2], 
		      double **stateArrivalProbs, double **oldProbs, 
		      double **newProbs, bool sumFactor);  
private:
  int K;
  int States; //number of states of Markov chain, m in book
  //There are k*k (=D in Chromosome)states since k populations and 2 chromosomes
  
  int Transitions; //length of chain
  // = # composite Loci, (=L in Chromosome)

  double sumfactor; 
  double factor;//used for adjustment of alpha to avoid underflow
  
  //forward and backward probabilities
  //L x K x K arrays
  double ***alpha, ***beta,**LambdaBeta;
  double ***ralpha;
  double *p,**Q ,**R, ***S;

  void PrintLikelihoods();
};

#endif /* ! HMM_H */
