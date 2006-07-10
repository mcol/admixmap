// *-*-C++-*-*
/** 
 *   DirichletParamSampler.h 
 *   Class to sample parameters of a Dirichlet distribution
 *   Copyright (c) 2005, 2006 Clive Hoggart, David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef DirichletParamSampler_H
#define DirichletParamSampler_H 1

#define DIRICHLETPARAM_ARS_SAMPLER 1
#define DIRICHLETPARAM_HAMILTONIAN_SAMPLER 2
#define DIRICHLETPARAM_SAMPLERTYPE 1 // 1 = adaptive rejection sampler on pairwise elements of proportion vector mu
                      //     Hamiltonian on dispersion parameter 
                      // 2 = Hamiltonian on parameter vector alpha

#include <gsl/gsl_sf_gamma.h>
#include "rand.h"
//#include "StepSizeTuner.h"
#include "HamiltonianMonteCarlo.h"

#if DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_ARS_SAMPLER
#include "AdaptiveRejection.h"
#endif

/// struct to hold arguments for sampler for Dirichlet parameters
typedef struct{
  int dim;
  int n;
  double eps0;
  double eps1;
  const double* sumlogtheta;
} AlphaSamplerArgs;
///struct to hold arguments for samper for population admixture dispersion parameter
typedef struct{
  unsigned numpops;
  unsigned numobs;
  const double* mu;
  const double* sumlogtheta;
  double priorshape;
  double priorrate;
} PopAdmixEtaSamplerArgs;

///Class to sample the parameters of a Dirichlet distribution
class DirichletParamSampler {
public:
  DirichletParamSampler();
  DirichletParamSampler(unsigned, unsigned);
  ~DirichletParamSampler();
  
  void SetSize( unsigned, unsigned);
  void SetPriorEta( double, double );
  void SetPriorMu( const double* const);
  void Sample( const double* const, std::vector<double> *alpha, bool fixedprops);
  double getStepSize()const;
  double getExpectedAcceptanceRate()const;
    
private:
  unsigned int K;

#if DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_ARS_SAMPLER
  //StepSizeTuner TuneEta;
  HamiltonianMonteCarlo EtaSampler;
  PopAdmixEtaSamplerArgs EtaArgs;

  double eta;
  double *mu;
  double etanew;
  double *munew;
  double *muDirichletParams;
//   double EtaAlpha;//shape parameter of prior on eta
//   double EtaBeta;//rate parameter of prior on eta
  //double step, step0;
  //double LogAccProb;
  AdaptiveRejection DirParamArray;
  double AlphaParameters[5];
  // AlphaParameters is an array with 5 elements
  // element 0 is number of observations
  // element 1 is dispersion parameter
  // element 2 is mu[k]
  // element 3 is sum log p[j]
  // element 4 is sum log p[k]

  //MuSampler muSampler;

  void SampleEta(unsigned n, const double* const sumlogtheta, double *eta, const double* const mu);
  static double logf( double, const void* const );
  static double dlogf( double, const void* const );
  static double ddlogf( double, const void* const );

  static double etaEnergy( const double* const eta, const void* const vargs );
  static void etaGradient( const double* const eta, const void* const vargs, double* g );

#elif DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_HAMILTONIAN_SAMPLER
  AlphaSamplerArgs AlphaArgs;
  double *logalpha;
  HamiltonianMonteCarlo AlphaSampler;
  double initialAlphaStepsize;
  float targetAlphaAcceptRate;
  static double alphaEnergy(const double* const theta, const void* const args);
  static void alphaGradient(const double* const theta, const void* const args, double *g);
#endif

  void Initialise();
};

#endif 
/* ! DirichletParamSampler_H */
