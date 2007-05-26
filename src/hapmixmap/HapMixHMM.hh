/*
 *   HAPMIXMAP
 *   HapMixHMM.hh 
 *   Extension of ADMIXMAP's HMM class with locus-specific Hidden-State probs
 *   Copyright (c) 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
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
  /**
     initialising constructor.
     \param inTransitions number of transitions (loci)
     \param numhiddenstates number of hidden states 
     \param f pointer to locus correlation, in Chromosome
  */
  HapMixHMM( int inTransitions, int numhiddenstates, const double* const f);
  ///Destructor
  ~HapMixHMM();

  /**
     sets HMM dimensions.
     initialising constructor.
     \param inTransitions number of transitions (loci)
     \param numhiddenstates number of hidden states 
     \param f pointer to locus correlation, in Chromosome
     \param diploid indicator for diploidy
  */
  void SetDimensions( int inTransitions, int pops, const double* const fin, bool diploid);

  ///sets genotype probs with genotype probs pointer
  void SetGenotypeProbs(const double* const lambdain, const bool* const missing);

  ///sets genotype probs with GPI
  void SetGenotypeProbs(const GenotypeProbIterator& lambdain, const bool* const missing);

  /**
     set state arrival probs.
     \param Theta mixture proportions
     \param Mcol Maternal gamete column (0 if assortative mating, 1 if random mating)
     \param isdiploid indicator for diploidy
     \see HiddenMarkovModel::SetStateArrivalProbs.
  */
  void SetStateArrivalProbs(const double* const Theta, const int Mcol, const bool isdiploid);

  /**
     Samples jump indicators.
     \param HiddenStates array of sampled hidden states
     \param gametes number of gametes (1 or 2)
     \param SumHiddenStates sums of hidden states where jump indicators are 1
  */
  void SampleJumpIndicators(const int* const HiddenStates, unsigned int gametes,
			    int *SumHiddenStates)const;

  ///returns loglikelihood
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
