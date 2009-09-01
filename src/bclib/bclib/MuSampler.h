// *-*-C++-*-*
/** 
 *   MuSampler.h 
 *   Class to sample the proportion parameters of a multinomial-Dirichlet distribution
 *   parameterised as \mu_1, ..., \mu_H, \eta ie a vector of proportions and a dispersion parameter
 *   Copyright (c) 2005, 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef MUSAMPLER_H
#define MUSAMPLER_H 1
#include "bclib/bclib.h"
#include "bclib/HamiltonianMonteCarlo.h"
#include "bclib/StepSizeTuner.h"

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */



typedef struct{
  const int *counts;
  const double *counts1;
  unsigned H;
  unsigned K;
  double eta;
} MuSamplerArgs;

class MuSampler{

public:
  MuSampler();
  ~MuSampler();

  void setDimensions(unsigned inK, unsigned inH, double mustep0, double mumin, double mumax, double mutarget);
  void Sample(double* const alpha, const double eta, const int* const Counts);
  void Sample1D(double* alpha, const double eta, const int* const Counts);

  float getAcceptanceRate()const;
  float getStepsize()const;

private:
  //dimensions
  unsigned K;//number of counts
  unsigned H;//dimension of alpha

  //proportions are transformed using softmax transformation and sampled with a Hamiltonian Monte Carlo sampelr
  double *params;
  HamiltonianMonteCarlo muSampler;
  MuSamplerArgs muArgs;

  static double muEnergyFunction(const double * const params, const void* const args);
  static void muGradient(const double * const params, const void* const args, double *g);
  static double fMu( double alpha, const void* const args );
  static double dfMu( double alpha, const void* const args );
  static double ddfMu( double alpha, const void* const args );
  static double logJacobian(const double* a, const double z, unsigned H);
  static double DlogJacobian(const double* const a, const double z, unsigned H, unsigned h, double delta);
};


/** @} */

END_BCLIB_NAMESPACE

#endif
