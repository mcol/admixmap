// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   HMC.h 
 *   Class to implement a Hamiltonian (or hybrid )Monte Carlo sampler
 *   (see Information Theory, Inference, and Learning Algorithms by David Mackay (1993), Neal (1993))
 *   Copyright (c) 2005 LSHTM
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

class HamiltonianMonteCarlo{
public:
  HamiltonianMonteCarlo();
  ~HamiltonianMonteCarlo();
  void Sample(double *x, const double* const* args);//call inside a loop
  void SetDimensions(unsigned pdim, double pepsilon, unsigned pTau, float target, 
		     double (*pfindE)(unsigned d, const double* const theta, const double* const* args),
		     void (*pgradE)(unsigned d, const double* const theta, const double* const *args, double *g));
  //sets dimension, stepsize, number of steps and parameters for density function
  float getAcceptanceRate();
  float getStepsize();
  void Tune();

private:
  double (*findE)(unsigned d, const double* const theta, const double* const* args); //calculate objective function
  void (*gradE)(unsigned d, const double* const theta, const double* const* args, double *g);//calculate gradient

  unsigned dim;     //dimension
  double epsilon;   //stepsize
  double Tau;       //# leapfrog steps
  double E;         //value of objective function
  double *g;        //gradient (multidim)
  long overall_accept_count;
  long totalsamples;//  "        "      "     "   "    "  in total
  StepSizeTuner Tuner;

  HamiltonianMonteCarlo(const HamiltonianMonteCarlo&);
  HamiltonianMonteCarlo& operator=(const HamiltonianMonteCarlo);

};


#endif /* !defined HMC_H */
