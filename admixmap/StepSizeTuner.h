// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   StepSizeTuner.h (formerly AdaptiveRandomWalkMH, TuneRW.h)
 *   This class is used to tune the step size in a Metropolis update, eg proposal sd in a Random Walk , 
 *   in order to reach a specified acceptance rate.
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#ifndef STEPSIZETUNER_H
#define STEPSIZETUNER_H 1

#include "rand.h"
#include <math.h>

class StepSizeTuner
{
public:
  /*supply: w = number of samplings between tunings
            step0 = initial step size
            min, max = min and max values for step size
            target = target acceptance rate
  */
  StepSizeTuner(int w, double sigma0, double min, double max, double target);
  StepSizeTuner();//default constructor, initialises k only
  ~StepSizeTuner();
  
  void SetParameters(int w, double sigma0, double min, double max, double target);
    
  double GetSigma();//returns log step size
  double UpdateStepSize(double AcceptanceProb);// returns step size after updating from current acceptance prob
  double getStepSize(); 
  double getExpectedAcceptanceRate();

private:
  double sigma0; // Initial value of stepsize.
  double sigma; // log step size 
  double step; // step size: must be positive
  double min; // Minimum value of stepsize
  double max; // Maximum value of stepsize
  double target; // Target acceptance probability
  double SumAcceptanceProb; // cumulative sum of acceptance probs
  int k; // Number of times stepsize is to be updated
  int w; // Frequency of updates - 10 is good
  int count; // Number of iterations since last update
  int NumberAccepted; // Number of accepted proposals in since last update
};

#endif /* ! StepSizeTuner_H */
