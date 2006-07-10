/* 
 *   StepSizeTuner.cc
 *   This class is used to tune the step size in a Metropolis update, eg proposal sd in a Random Walk , 
 *   in order to reach a specified acceptance rate.
 *   Copyright (c) 2005, 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "StepSizeTuner.h"

StepSizeTuner::StepSizeTuner()
{
  k = 1;   
}


StepSizeTuner::StepSizeTuner(double insigma0, double inmin, double inmax, double intarget)
{
  SetParameters(insigma0, inmin, inmax, intarget);
}

StepSizeTuner::~StepSizeTuner()
{
}

void StepSizeTuner::resetApproximator(int newk) {
  if(k > newk) k = newk;   
}

void StepSizeTuner::SetParameters(double step0, double inmin, double inmax, double intarget)
{
  // step, inmin, inmax must be positive real numbers
  step = step0;
  sigma0 = log(step);
  sigma = sigma0;
  min = log(inmin);
  max = log(inmax);
  target = intarget;
  k = 1;
  count = 0;
  SumAcceptanceProb = 0.0;
}

double StepSizeTuner::GetSigma()const
{
  return sigma;
}

double StepSizeTuner::UpdateStepSize(double AcceptanceProb)
{
  if(AcceptanceProb > 1.0) AcceptanceProb = 1.0;
  sigma = sigma + ( AcceptanceProb - target ) / k; 
  /// initial adjustment to step size will be of the order of exp(0.5) fold 
  if( sigma > max )
    sigma = max;
  else if( sigma < min )
    sigma = min;
  step = exp(sigma);
  SumAcceptanceProb += AcceptanceProb; // accumulate sum of acceptance probs
  count++;
  k++;
  return step; /// returns step as a positive real number
}

double StepSizeTuner::getStepSize()const
{
  return step; 
}

double StepSizeTuner::getExpectedAcceptanceRate()const
{
  if(count > 0)
    return SumAcceptanceProb / count;
  else return 0.0;
}


