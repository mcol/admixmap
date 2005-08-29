/** 
 *   ADMIXMAP
 *   GaussianProposalMH.cc (formerly TuneRW.cc)
 *   This class is used to implement a Metropolis Hastings update with Gaussian proposal distribution
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
#include "AdaptiveRandomWalkMH.h"


AdaptiveRandomWalkMH::AdaptiveRandomWalkMH()
{
  k = 1;   
}

AdaptiveRandomWalkMH::AdaptiveRandomWalkMH(int inw, double insigma0, double inmin, double inmax, double intarget)
{
  SetParameters(inw, insigma0, inmin, inmax, intarget);
}

AdaptiveRandomWalkMH::~AdaptiveRandomWalkMH()
{
}

void AdaptiveRandomWalkMH::SetParameters(int inw, double step0, double inmin, double inmax, double intarget)
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

double AdaptiveRandomWalkMH::GetSigma()
{
  return sigma;
}

double AdaptiveRandomWalkMH::UpdateStepSize(double AcceptanceProb)
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

double AdaptiveRandomWalkMH::UpdateSigma(int NumberAccepted)
{
  double ProportionAccepted = (double)NumberAccepted / w;
  sigma = sigma + sigma0 * ( ProportionAccepted - target ) / k;
  if( sigma > max )
    sigma = max;
  else if( sigma < min )
    sigma = min;
  k++;
  return sigma;
}

double AdaptiveRandomWalkMH::getStepSize()
{
  return step; 
}

double AdaptiveRandomWalkMH::getExpectedAcceptanceRate()
{
  return SumAcceptanceProb / count;
}

void AdaptiveRandomWalkMH::Event(bool accept)
{
  count++;
  NumberAccepted+=accept;
  if( count == w ){
    double ProportionAccepted = (double)NumberAccepted / w;
    sigma = sigma + sigma0 * ( ProportionAccepted - target ) / k;
    if( sigma > max )
      sigma = max;
    else if( sigma < min )
      sigma = min;
    k++;
    NumberAccepted = 0;
    count = 0;
  }  
}
