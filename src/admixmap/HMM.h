// *-*-C++-*-*
/*
 *   HMM.h 
 *   Class to implement hidden Markov model for haploid or diploid Poisson arrivals.
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef HMM_H
#define HMM_H 1

#include "GPI.h"
#include "utils/pvector.h"
#include "MixturePropsWrapper.hh"

///
class HMM
{
public:
  HMM();
  HMM( int inTransitions, int pops, const double* const f);
  ~HMM();
  void SetDimensions( int inTransitions, int pops, const double* const fin);
  void SetGenotypeProbs(const GenotypeProbIterator& lambdain, const bool* const missing);
  void SetTheta(const MixturePropsWrapper& Theta, const MixturePropsWrapper& ThetaSq, 
		const MixturePropsWrapper& ThetaSqInv);
  void SetStateArrivalProbs(int Mcol, bool isdiploid);

  void Sample(int *SStates, bool isdiploid);
  const pvector<double>& GetHiddenStateProbs( bool isDiploid, int t );
  double getLogLikelihood(bool isDiploid);
  void SampleJumpIndicators(const int* const LocusAncestry, unsigned int gametes, int *SumLocusAncestry, 
			    std::vector<unsigned> &SumNumArrivals, bool SampleArrivals, unsigned startlocus)const;
  void SampleJumpIndicators(const int* const HiddenStates, unsigned int gametes, 
			    int *SumHiddenStates)const;
private:
  int K;
  int DStates; ///<number of diploid states 
  //There are K*K (=D in Chromosome) states since K populations and 2 chromosomes
  
  int Transitions; ///<length of chain
  // = # composite Loci, (=L in Chromosome)

  double sumfactor; 
  
  //forward and backward probabilities
  //L x K x K arrays
  double *alpha, *beta, *LambdaBeta;
  double *p;
  double *StateArrivalProbs;

  const double* f;
  MixturePropsWrapper theta;
  MixturePropsWrapper ThetaThetaPrime;
  MixturePropsWrapper ThetaThetaInv;
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
  void RecursionProbs(double ff, const double f[2], const double* const stateArrivalProbs,
		      const double* const oldProbs, double *newProbs); 
  void RecursionProbs2(double ff, const double f[2], const double* const stateArrivalProbs, 
		       const double* const oldProbs, double *newProbs);
  
  // Storing results so vectors are not being created every time
  // the GetHiddenStateProbs() function is called
  pvector<double> hiddenStateProbs;
};

#endif /* ! HMM_H */
