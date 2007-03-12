// *-*-C++-*-*

#ifndef HMM_H
#define HMM_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string>
#include <cstdlib>
#include <values.h>
#include <vector>
#include "GPI.h"
#include "MixturePropsWrapper.hh"

/** Class to implement hidden Markov model for haploid or diploid Poisson arrivals.
   Instantiated in class Chromosome
*/

class HMM
{
public:
  HMM();
  HMM( int inTransitions, int pops, const double* const f);
  ~HMM();
  void SetDimensions( int inTransitions, int pops, const double* const fin);
  void SetGenotypeProbs(const GenotypeProbIterator& lambdain, const bool* const missing);
  void SetTheta(const MixturePropsWrapper& Theta, const int Mcol, const bool isdiploid);
  void SetStateArrivalProbs(const int Mcol, bool isdiploid);

  void Sample(int *SStates, bool isdiploid);
  const std::vector<double> GetHiddenStateProbs( const bool isDiploid, int t );
  double getLogLikelihood(bool isDiploid);
  void SampleJumpIndicators(const int* const LocusAncestry, const unsigned int gametes, int *SumLocusAncestry, 
			    std::vector<unsigned> &SumNumArrivals, bool SampleArrivals, unsigned startlocus)const;
  void SampleJumpIndicators(const int* const HiddenStates,  const unsigned int gametes, 
			    int *SumHiddenStates)const ;
private:
  int K;
  int DStates; //number of diploid states 
  //There are K*K (=D in Chromosome) states since K populations and 2 chromosomes
  
  int Transitions; //length of chain
  // = # composite Loci, (=L in Chromosome)

  double sumfactor; 
  
  //forward and backward probabilities
  //L x K x K arrays
  double *alpha, *beta, *LambdaBeta;
  double *p;
  double *StateArrivalProbs;
  double *ThetaThetaPrime;

  const double* f;
  MixturePropsWrapper theta;
  GenotypeProbIterator LambdaGPI;
  const bool* missingGenotypes;
  bool alphaIsBad, betaIsBad;

  //arrays for RecursionProbs
  double *rowProb;
  double *colProb;
  double *Expectation0;
  double *Expectation1;
  double *cov;

  void UpdateForwardProbsDiploid();
  void UpdateForwardProbsHaploid();
  void UpdateBackwardProbsDiploid();
  void UpdateBackwardProbsHaploid();
  void RecursionProbs(const double ff, const double f[2], const double* const stateArrivalProbs,
		      const double* const oldProbs, double *newProbs); 
  void RecursionProbs2(const double ff, const double f[2], const double* const stateArrivalProbs, 
		       const double* const oldProbs, double *newProbs);
};

#endif /* ! HMM_H */
