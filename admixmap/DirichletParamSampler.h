// *-*-C++-*-*
#ifndef DirichletParamSampler_H
#define DirichletParamSampler_H 1

#include <gsl/gsl_sf_gamma.h>
#include "TuneRW.h"
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
  
private:
  TuneRW TuneEta;
  TuneRW TuneMu;
  unsigned int d;
  double etanew;
  double *munew;
  double *gamma;
  double *alpha;
  double EtaAlpha;
  double EtaBeta;
};

#endif /* ! DirichletParamSampler_H */
