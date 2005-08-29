// *-*-C++-*-*
#ifndef DirichletParamSampler_H
#define DirichletParamSampler_H 1

#include <gsl/gsl_sf_gamma.h>
#include "DARS.h"
#include "StepSizeTuner.h"
#include "rand.h"

class DirichletParamSampler
{
public:
  DirichletParamSampler();
  DirichletParamSampler(unsigned int);
  ~DirichletParamSampler();
  
  void SetSize( unsigned int );
  void SetPriorEta( double, double );
  void SetPriorMu( double* );
  void Sample( unsigned int, double*, double*, double* );
  double getStepSize();
  double getExpectedAcceptanceRate();
  
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
  // element 0 is 
  // element 1 is 
  // 
  // element 4 is 
  double AlphaParameters[5];
  static double
  logf( const double*, const int*, const double* , double );
  
  static double
  dlogf( const double*, const int *, const double* , double );
  
  static double
  ddlogf( const double*, const int *, const double* , double );
};

#endif /* ! DirichletParamSampler_H */
