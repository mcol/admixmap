// *-*-C++-*-*
#ifndef DirichletParamSampler_H
#define DirichletParamSampler_H 1

#include <gsl/gsl_sf_gamma.h>
#include "DARS.h"
#include "StepSizeTuner.h"
#include "MuSampler.h"
#include "rand.h"

class DirichletParamSampler
{
public:
  DirichletParamSampler();
  DirichletParamSampler(unsigned, unsigned);
  ~DirichletParamSampler();
  
  void SetSize( unsigned, unsigned );
  void SetPriorEta( double, double );
  void SetPriorMu( double* );
  void Sample( unsigned int, double*, double*, double* );
  void Sample2( unsigned int, double*, double*, double*, int* );
  double getEtaStepSize();
  double getEtaExpectedAcceptanceRate();
  double getMuStepSize();
  double getMuExpectedAcceptanceRate();
  
private:
  StepSizeTuner TuneEta;
  unsigned int d;
  double etanew;
  double *munew;
  double *gamma;
  double EtaAlpha;
  double EtaBeta;
  double step, step0;
  double LogAccProb;

  DARS** DirParamArray;
  // AlphaParameters is an array with 5 elements
  // element 0 is number of observations
  // element 1 is dispersion parameter
  // element 2 is 1 - last proportion parameter
  // element 3 is 
  // element 4 is 
  double AlphaParameters[5];

  MuSampler muSampler;

  void SampleEta(unsigned n, double *sumlogtheta, double *eta, double *mu);

  static double logf( const double*, const int*, const double* , double );
  
  static double dlogf( const double*, const int *, const double* , double );
  
  static double ddlogf( const double*, const int *, const double* , double );
};

#endif /* ! DirichletParamSampler_H */
