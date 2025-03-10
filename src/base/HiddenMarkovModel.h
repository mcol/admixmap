//=============================================================================
//
// Copyright (C) 2002-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file HiddenMarkovModel.h
/// Definition of the HiddenMarkovModel class.
//=============================================================================

#ifndef HIDDENMARKOVMODEL_H
#define HIDDENMARKOVMODEL_H 1



/** \addtogroup base
 * @{ */


#include <vector>
#include "bclib/pvector.h"


class GenotypeProbIterator;


/** Class to implement hidden Markov model for haploid or diploid Poisson arrivals.
    Resides in class Chromosome
*/
class HiddenMarkovModel {
public:
  /**
     initialising constructor.
     \param inTransitions number of transitions (loci on chromosome)
     \param pops number of hidden states
     \param f pointer to locus correlation (held in Chromosome
  */
  HiddenMarkovModel( int _transitions, size_t _K, size_t _nStates, const double * _f);
  ///Destructor
  virtual ~HiddenMarkovModel();

  /**
     sets pointer to genotype probs.
     \param lambdain the genotype probs
     \param missing indicators for missing genotypes
  */
  void SetGenotypeProbs(const double* const lambdain, const bool* const missing);

  ///set genotype probs 
  //temporary until interface is finished
  //final version will use bcppcl::arrays instead of GPI
  virtual void SetGenotypeProbs(const GenotypeProbIterator& , const bool* const ){};

  /**
     set state arrival probs.
     \param Theta mixture proportions
     \param Mcol Maternal gamete column (0 if assortative mating, 1 if random mating)
     \param isdiploid indicator for diploidy
  */
  virtual void SetStateArrivalProbs(const double* const Theta, const int Mcol, const bool isdiploid);

  /**
     Samples Hidden States.
     \param SStates an int array to store the sampled states
     \param isDiploid indicator for diploidy
  */
  void SampleHiddenStates(int *SStates, bool isdiploid);

  ///returns conditional probs of hidden states
  const bclib::pvector<double>& GetHiddenStateProbs( const bool isDiploid, int t );

  ///returns loglikelihood
  double getLogLikelihood(bool isDiploid);

  /**
     Samples jump indicators.
     \param LocusAncestry array of hidden states
     \param gametes number of gametes (1 or 2)
     \param SumLocusAncestry
     \param SumNumArrivals
     \param SampleArrivals indicates whether to sample arrivals or not
     \param startlocus index of first locus on chromosome
  */
  void SampleJumpIndicators(const int* const LocusAncestry, const unsigned int gametes, int *SumLocusAncestry, 
			    std::vector<unsigned> &SumNumArrivals, bool SampleArrivals, unsigned startlocus)const;

  //required for HapMixHMM class
  virtual void SampleJumpIndicators(const int* const , unsigned int, int *)const{};

  /// Ensure that the forward/backward probabilities are up-to-date
  void UpdateForwardBackwardProbs(bool isDiploid) {
    UpdateForwardProbs(isDiploid);
    UpdateBackwardProbs(isDiploid);
  }

protected:

  /// Number of populations
  const int K;
  const int nStates; ///< Number of hidden states (for diploid individual, K*K)
  const int Transitions; //length of chain
  // = # composite Loci, (=L in Chromosome)

  double sumfactor; ///< for accumulating log-likelihood
  
  //forward and backward probabilities
  //L x K x K arrays
  double *alpha, *beta, *LambdaBeta;
  double *p;
  double* StateArrivalProbs[2];

  /// Stationary distribution
  double* pi;

  /// Inverse of the stationary distribution
  double* piInv;

  /// Locus correlations
  const double* f;

  /// Admixture probabilities
  const double* theta;

  /// Emission probabilities
  const double* Lambda;
  const bool* missingGenotypes;
  bool alphaIsBad, betaIsBad;

  // Storing results so vectors are not being created every time
  // the GetHiddenStateProbs() function is called
  bclib::pvector<double> hiddenStateProbs;

  void UpdateForwardProbs(bool);
  void UpdateBackwardProbs(bool);
  
  virtual void UpdateForwardProbsDiploid();
  virtual void UpdateForwardProbsHaploid();
  virtual void UpdateBackwardProbsDiploid();
  virtual void UpdateBackwardProbsHaploid();

  void SetArraysForRecursionProbs();
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

  /// Allocate the arrays and set the pointer of locus correlations f
  virtual void SetDimensions(const double* const fin);

  // UNIMPLEMENTED
  // to avoid use
  ///default constructor
  HiddenMarkovModel();
  HiddenMarkovModel(const HiddenMarkovModel&);
  HiddenMarkovModel& operator=(const HiddenMarkovModel&);

};


/** @} */


#endif /* ! HIDDENMARKOVMODEL_H */
