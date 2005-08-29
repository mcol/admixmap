/** 
 *   ADMIXMAP
 *   StepSizeTuner.cc (formerly AdaptiveRandomWalkMH.cc, TuneRW.cc)
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
#include "StepSizeTuner.h"


StepSizeTuner::StepSizeTuner()
{
  k = 1;   
}

StepSizeTuner::StepSizeTuner(int inw, double insigma0, double inmin, double inmax, double intarget)
{
  SetParameters(inw, insigma0, inmin, inmax, intarget);
}

StepSizeTuner::~StepSizeTuner()
{
}

void StepSizeTuner::SetParameters(int inw, double step0, double inmin, double inmax, double intarget)
{
  w = inw;
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

double StepSizeTuner::GetSigma()
{
  return sigma;
}

double StepSizeTuner::UpdateStepSize(double AcceptanceProb)
{
  sigma = sigma + ( AcceptanceProb - target ) / k; 
  // initial adjustment to step size will be of the order of exp(0.5) fold 
  if( sigma > max )
    sigma = max;
  else if( sigma < min )
    sigma = min;
  step = exp(sigma);
  SumAcceptanceProb += AcceptanceProb; // accumulate sum of acceptance probs
  count++;
  k++;
  return step;
}

double StepSizeTuner::getStepSize()
{
  return step; 
}

double StepSizeTuner::getExpectedAcceptanceRate()
{
  if(count > 0)
    return SumAcceptanceProb / count;
  else return 0.0;
}


