// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   AdaptiveRandomWalkMH.h (formerly TuneRW.h)
 *   This class is used to implement an adaptive random-walk Metropolis Hastings sampler, which tunes a proposal variance 
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
#ifndef ADAPTIVERANDOMWALKMH_H
#define ADAPTIVERANDOMWALKMH_H 1

#include "rand.h"
#include <math.h>

class AdaptiveRandomWalkMH
{
public:
  /*supply: w = number of samplings between tunings
            step0 = initial step size
            min, max = min and max values for step size
            target = target acceptance rate
  */
  AdaptiveRandomWalkMH(int w, double sigma0, double min, double max, double target);
  AdaptiveRandomWalkMH();//default constructor, initialises k only
  ~AdaptiveRandomWalkMH();
  
  void SetParameters(int w, double sigma0, double min, double max, double target);
    
  double GetSigma();//returns log step size
  double UpdateSigma(int);
  double UpdateStepSize(double AcceptanceProb);// returns step size after updating from current acceptance prob
  void Event(bool);// updates acceptance count and tunes proposal SD
  double getStepSize(); 
  double getExpectedAcceptanceRate();

private:
  double sigma0; // Initial value of parameter of random walk being tuned.
  double sigma; // log step size 
  double step; // step size: must be positive
  double min; // Minimum value of sigma
  double max; // Maximum value of sigma
  double target; // Target acceptance probability
  double SumAcceptanceProb; // cumulative sum of acceptance probs
  int k; // Number of times sigma sigma updated
  int w; // Frequency sigma is updated - 10 is good
  int count; // Number of iterations since sigma lsat updated
  int NumberAccepted; // Number of accepted proposals in random-walk since sigma last updated
};

#endif /* ! AdaptiveRandomWalkMH_H */
