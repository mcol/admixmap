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

class HMM
{
public:
  HMM();
  HMM( int inTransitions, int pops);
  ~HMM();
  void SetDimensions( int inTransitions, int pops);

  void SetStateArrivalProbs(const double* const f[], const double* const Theta, int Mcol);
  /* samples hidden states */
  void Sample(int *SStates, const double* const Admixture, const double* const f[], bool isdiploid)const;
  //void GetStateProbs( double * probs, int t)const;
  std::vector<std::vector<double> > Get3WayStateProbs( int t)const;

  double getLogLikelihood()const;

  void UpdateForwardProbsDiploid(const double* const f[], const double* const lambda, double coolness);

  void UpdateBackwardProbsDiploid(const double* const f[], const double* const lambda);

  void UpdateForwardProbsHaploid(const double* const f[], const double* const Admixture, const double* const lambda);
  void UpdateBackwardProbsHaploid(const double* const f[], const double* const Admixture, const double* const lambda);

  void RecursionProbs(const double ff, const double f[2], const double* const stateArrivalProbs,
		      double* oldProbs, double *newProbs); 

  void SampleJumpIndicators(const int* const LocusAncestry, const double* const f[], const unsigned int gametes, 
			    int *SumLocusAncestry, int *SumLocusAncestry_X, bool isX, 
			    unsigned int SumN[], unsigned int SumN_X[], bool isGlobalRho)const;
private:
  int K;
  int States; //number of states of Markov chain, m in book
  //There are k*k (=D in Chromosome)states since k populations and 2 chromosomes
  
  int Transitions; //length of chain
  // = # composite Loci, (=L in Chromosome)

  double sumfactor; 
  
  //forward and backward probabilities
  //L x K x K arrays
  double *alpha, *beta,*LambdaBeta;
  double *p;
  double *StateArrivalProbs;
  double *ThetaThetaPrime;

  //arrays for RecursionProbs
  double *rowProb;
  double *colProb;
  double *Expectation0;
  double *Expectation1;
  double *rowSum;
  double *colSum;
  double **cov;

void RecursionProbs2(const double ff, const double f[2], 
		     const double* const stateArrivalProbs, const double* const oldProbs, double *newProbs) ;

};

#endif /* ! HMM_H */
