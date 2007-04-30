/*
 *   HAPMIXMAP
 *   HapMixHMM.hh 
 *   Extension of ADMIXMAP's HMM class with locus-specific Hidden-State probs
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef HAPMIXHMM_H
#define HAPMIXHMM_H 1

#include "HiddenMarkovModel.h"
#include "GPI.h"

///Extension of ADMIXMAP's HMM class with locus-specific Hidden-State probs
class HapMixHMM :public HiddenMarkovModel{

public:
  HapMixHMM( int inTransitions, int pops, const double* const f);
  ~HapMixHMM();

  void SetDimensions( int inTransitions, int pops, const double* const fin, bool diploid);

  void SetGenotypeProbs(const double* const lambdain, const bool* const missing);

  void SetGenotypeProbs(const GenotypeProbIterator& lambdain, const bool* const missing);

  void SetStateArrivalProbs(const double* const Theta, const int Mcol, const bool isdiploid);

  void SampleJumpIndicators(const int* const HiddenStates, unsigned int gametes,
			    int *SumHiddenStates)const;
  double getLogLikelihood(bool isDiploid){return HiddenMarkovModel::getLogLikelihood(isDiploid);};
private:
  GenotypeProbIterator LambdaGPI;


  //overloaded update functions
  void UpdateForwardProbsDiploid();
  void UpdateForwardProbsHaploid();
  void UpdateBackwardProbsDiploid();
  void UpdateBackwardProbsHaploid();

  // UNIMPLEMENTED
  // to avoid use
  HapMixHMM();
  HapMixHMM(const HapMixHMM&);
  HapMixHMM& operator=(const HapMixHMM&);
};


#endif
