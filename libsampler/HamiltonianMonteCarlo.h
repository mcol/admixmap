// *-*-C++-*-*
/** 
 *   HamiltonianMonteCarlo.h 
 *   Copyright (c) 2005, 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef HMC_H
#define HMC_H 1
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "rand.h"
#include "StepSizeTuner.h"
extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
}

/**
 *   Class to implement a Hamiltonian (or hybrid )Monte Carlo sampler.
 *   See Information Theory, Inference, and Learning Algorithms by David Mackay (1993), Neal (1993).
*/
class HamiltonianMonteCarlo{
public:
  HamiltonianMonteCarlo();
  ~HamiltonianMonteCarlo();
  /// Sample new value
  void Sample(double* const x, const void* const args);
  /// Set dimension, stepsize, number of steps and parameters for density function.
  void SetDimensions(unsigned pdim, double pepsilon, double min, double max, unsigned pTau, float target, 
		     double (*pfindE)(const double* const theta, const void* const args),
		     void (*pgradE)( const double* const theta, const void* const args, double *g));
  ///returns expected acceptance rate
  float getAcceptanceRate()const;
  ///returns current stepsize
  float getStepsize()const;

private:
  double (*findE)(const double* const theta, const void* const args); 
  void (*gradE)(const double* const theta, const void* const args, double *g);

  unsigned dim;     //<dimension
  double epsilon;   //<stepsize
  double Tau;       //<# leapfrog steps
  double* xnew, *g, *gnew, *p;
  long overall_accept_count;
  long totalsamples;//  "        "      "     "   "    "  in total
  StepSizeTuner Tuner;

  HamiltonianMonteCarlo(const HamiltonianMonteCarlo&);
  HamiltonianMonteCarlo& operator=(const HamiltonianMonteCarlo);

};


#endif /* !defined HMC_H */
