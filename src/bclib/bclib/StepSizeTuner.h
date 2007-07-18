// *-*-C++-*-*
/** 
 *   StepSizeTuner.h
 *   Copyright (c) 2005, 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef STEPSIZETUNER_H
#define STEPSIZETUNER_H 1

#include "bclib/bclib.h"
#include "bclib/rand.h"
#include <math.h>

BEGIN_BCLIB_NAMESPACE

/**
 *   This class is used to tune the step size in a Metropolis update, eg proposal sd in a Random Walk , 
 *   in order to reach a specified acceptance rate.
*/
class StepSizeTuner
{
public:
  /** Constructor.
     supply: 
            step0 = initial step size
            min, max = min and max values for step size
            target = target acceptance rate
  */
  StepSizeTuner(double sigma0, double min, double max, double target);
  /// default constructor, initialises k only
  StepSizeTuner();
  ~StepSizeTuner();
  
  void resetStepSizeApproximator(int newk);
  void SetParameters(double sigma0, double min, double max, double target);
  ///returns log step size
  double GetSigma()const;
  /// returns step size after updating from current acceptance prob
  double UpdateStepSize(double AcceptanceProb);
  double getStepSize()const; 
  double getExpectedAcceptanceRate()const;

private:
  double sigma0; // Initial value of stepsize.
  double step; // step size: must be positive
  double sigma; // log step size 
  double min; // Minimum value of stepsize
  double max; // Maximum value of stepsize
  double target; // Target acceptance probability
  double SumAcceptanceProb; // cumulative sum of acceptance probs
  int k; // Number of times stepsize is to be updated
  int count; // Number of iterations since last update
  int NumberAccepted; // Number of accepted proposals in since last update
};

END_BCLIB_NAMESPACE
#endif /* ! StepSizeTuner_H */
