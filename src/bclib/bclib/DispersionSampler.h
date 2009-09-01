// *-*-C++-*-*
/** 
 *   DispersionSampler.h 
 *   Class to sample the dispersion parameter of a multinomial-Dirichlet distribution
 *   parameterised as \mu_1, ..., \mu_H, \eta ie a vector of proportions and a dispersion parameter
 *   Copyright (c) 2005, 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef DISPERSIONSAMPLER_H
#define DISPERSIONSAMPLER_H 1
#include "bclib/bclib.h"
#include "bclib/HamiltonianMonteCarlo.h"

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */



///struct to hold the arguments for sampler for dispersion parameters
typedef struct{
  const int** counts;
  const double** alpha;
  double priorshape;
  double priorrate;
} EtaSamplerArgs;

///Class to sample the dispersion parameter of a multinomial-Dirichlet distribution
class DispersionSampler{

public:
  DispersionSampler();
  ~DispersionSampler();

  void Initialise(double step0, double min, double max, double target);
  static void setDimensions(unsigned inL, unsigned inK, int* const inH);
  void setEtaPrior(double, double);
  void addAlphas(unsigned, const double* const);
  void addCounts(unsigned, const int* const);
  double Sample();
  double getEnergy(double);
  double getGradient(double);

  float getAcceptanceRate()const;
  float getStepsize()const;

private:
  //dimensions
  static unsigned L;
  static unsigned K;//number of counts
  static unsigned *NumStates;

  double logeta[1];
  HamiltonianMonteCarlo Sampler;
  EtaSamplerArgs Args;

  static double etaEnergyFunction(const double * const logitmu, const void* const args);
  static void etaGradient(const double * const logitmu, const void* const args, double *g);
};


/** @} */

END_BCLIB_NAMESPACE

#endif
