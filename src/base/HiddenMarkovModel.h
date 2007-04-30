// *-*-C++-*-*
/*
 *   HiddenMarkovModel.h 
 *   Class to implement hidden Markov model for haploid or diploid Poisson arrivals.
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef HIDDENMARKOVMODEL_H
#define HIDDENMARKOVMODEL_H 1

#include <vector>
#include "bcppcl/pvector.h"

class GenotypeProbIterator;

/** Class to implement hidden Markov model for haploid or diploid Poisson arrivals.
    Resides in class Chromosome
*/
class HiddenMarkovModel{
public:
  HiddenMarkovModel();
  HiddenMarkovModel( int inTransitions, int pops, const double* const f);
  virtual ~HiddenMarkovModel();

  void SetGenotypeProbs(const double* const lambdain, const bool* const missing);
  //temporary until interface is finished
  //final version will use bcppcl::arrays instead of GPI
  virtual void SetGenotypeProbs(const GenotypeProbIterator& , const bool* const ){};

  virtual void SetStateArrivalProbs(const double* const Theta, const int Mcol, const bool isdiploid);

  void SampleHiddenStates(int *SStates, bool isdiploid);
  const bcppcl::pvector<double>& GetHiddenStateProbs( const bool isDiploid, int t );
  double getLogLikelihood(bool isDiploid);
  void SampleJumpIndicators(const int* const LocusAncestry, const unsigned int gametes, int *SumLocusAncestry, 
			    std::vector<unsigned> &SumNumArrivals, bool SampleArrivals, unsigned startlocus)const;

  //required for HapMixHMM class
  virtual void SampleJumpIndicators(const int* const , unsigned int, int *)const{};

protected:
  int K;
  int DStates; //number of diploid states  = K*K
  
  int Transitions; //length of chain
  // = # composite Loci, (=L in Chromosome)

  double sumfactor;///< for accumulating log-likelihood
  
  //forward and backward probabilities
  //L x K x K arrays
  double *alpha, *beta, *LambdaBeta;
  double *p;
  double* StateArrivalProbs[2];
  double* ThetaThetaPrime;
  double* ThetaThetaInv;

  const double* f;
  const double* theta;
  const double* Lambda;
  const bool* missingGenotypes;
  bool alphaIsBad, betaIsBad;

  // Storing results so vectors are not being created every time
  // the GetHiddenStateProbs() function is called
  bcppcl::pvector<double> hiddenStateProbs;

  virtual void UpdateForwardProbsDiploid();
  virtual void UpdateForwardProbsHaploid();
  virtual void UpdateBackwardProbsDiploid();
  virtual void UpdateBackwardProbsHaploid();

  void SetArraysForRecursionProbs(unsigned K);
  void SetNullValues();
  void RecursionProbs(const double ff, const double f[2], const double* const stateArrivalProbs0,const double* const stateArrivalProbs1,
		      const double* const oldProbs, double *newProbs); 
  void RecursionProbs2(const double ff, const double f[2], const double* const stateArrivalProbs0,const double* const stateArrivalProbs1, 
		       const double* const oldProbs, double *newProbs);

private:
  //arrays for RecursionProbs
  double *rowProb;
  double *colProb;
  double *Expectation0;
  double *Expectation1;
  double *cov;

  virtual void SetDimensions( int inTransitions, int pops, const double* const fin, bool diploid);

  // UNIMPLEMENTED
  // to avoid use
  HiddenMarkovModel(const HiddenMarkovModel&);
  HiddenMarkovModel& operator=(const HiddenMarkovModel&);

};

#endif /* ! HIDDENMARKOVMODEL_H */
