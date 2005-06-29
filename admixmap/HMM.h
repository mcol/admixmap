// *-*-C++-*-*

/* Class to implement Hidden Markov Models (see MacDonald and Zucchini)
   Instantiated in class Chromosome
*/


#ifndef HMM_H
#define HMM_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
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

  void SetStateArrivalProbs(double *f[], double *Theta, int Mcol);
  /* samples hidden states */
  void Sample(int *SStates, double *Admixture, double *f[], bool isdiploid);
  void GetStateProbs( double * probs, int t);
  void Get3WayStateProbs( int j, double AncestryProbs[][3]);

  double getLikelihood();

  void UpdateForwardProbsDiploid(double *f[], double *lambda);

  void UpdateBackwardProbsDiploid(double *f[], double *lambda);

  void UpdateProbsHaploid(double *f[], double *Admixture, double *lambda, bool CalculateBeta);

  void RecursionProbs(const double ff, const double f[2], double *stateArrivalProbs,
		      double *oldProbs, double *newProbs); 

  void SampleJumpIndicators(int *LocusAncestry, double *f[], const unsigned int gametes, 
			    double *Distances, const int startLocus,  
			    int *sumxi, double *Sumrho0, int *SumLocusAncestry, int *SumLocusAncestry_X, bool isX, 
			    unsigned int SumN[], unsigned int SumN_X[], bool RhoIndicator) ;
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
  double *alpha, *beta,*LambdaBeta;
  double *p;
  double *StateArrivalProbs;
  double *ThetaThetaPrime;

};

#endif /* ! HMM_H */
