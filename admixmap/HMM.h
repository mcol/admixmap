// *-*-C++-*-*

/* Class to implement hidden Markov model for haploid or diploid Poisson arrivals
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
#include <vector>

class HMM
{
public:
  HMM();
  HMM( int inTransitions, int pops);
  ~HMM();
  void SetDimensions( int inTransitions, int pops);
  void SetStateArrivalProbs(const double* const f, const double* const Theta, int Mcol);

  void Sample(int *SStates, const double* const Admixture, const double* const f, bool isdiploid)const;
  std::vector<std::vector<double> > Get3WayStateProbs( const bool isDiploid, int t )const;
  double getLogLikelihood(bool isDiploid)const;
  void SampleJumpIndicators(const int* const LocusAncestry, const double* const f, const unsigned int gametes, 
			    int *SumLocusAncestry, std::vector<unsigned> &SumNumArrivals, bool SampleArrivals, unsigned startlocus)const;


  // these four methods should be private
  void UpdateForwardProbsDiploid(const double* const f, const double* lambda, bool* const missing);
  void UpdateBackwardProbsDiploid(const double* const f, const double* const lambda);
  void UpdateForwardProbsHaploid(const double* const f, const double* const Admixture, const double* const lambda);
  void UpdateBackwardProbsHaploid(const double* const f, const double* const Admixture, const double* const lambda);

private:
  int K;
  int DStates; //number of diploid states 
  //There are K*K (=D in Chromosome) states since K populations and 2 chromosomes
  
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

  void RecursionProbs(const double ff, const double f[2], const double* const stateArrivalProbs,
		      double* oldProbs, double *newProbs); 
  void RecursionProbs2(const double ff, const double f[2], const double* const stateArrivalProbs, 
		       const double* const oldProbs, double *newProbs);
};

#endif /* ! HMM_H */
