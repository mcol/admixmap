// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   HMCMC.h 
 *   Class to implement a Hamiltonian (or hybrid )MCMC sampler
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
#ifndef HMCMC_H
#define HMCMC_H 1
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "rand.h"
extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
}

class HMCMC{
public:
  HMCMC();
  ~HMCMC();
  void Sample(double *x, double **args);//call inside a loop
  void SetDimensions(const unsigned pdim, const double pepsilon, const unsigned pTau, 
		     double (*pfindE)(unsigned d, double *theta, double **args),
		     void (*pgradE)(unsigned d, double *theta, double **args, double *g));
  //sets dimension, stepsize, number of steps and parameters for density function
  float getAcceptanceCount();

private:
  double (*findE)(unsigned d, double *theta, double **args); //calculate objective function
  void (*gradE)(unsigned d, double *theta, double **args, double *g);//calculate gradient

  unsigned dim;     //dimension
  double epsilon;   //stepsize
  double Tau;       //# leapfrog steps
  double E;         //value of objective function
  double *g;        //gradient (multidim)
  long accept_count; //number of acceptances, incase needed for monitoring

};


#endif /* !defined HMCMC_H */
