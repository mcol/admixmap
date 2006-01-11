// *-*-C++-*-*
#ifndef DirichletParamSampler_H
#define DirichletParamSampler_H 1

#define SAMPLERTYPE 1// 1 = ARS + MHRW
                     // 2 = Hamiltonian 


#include <gsl/gsl_sf_gamma.h>
#include "rand.h"

#if SAMPLERTYPE==1
#include "AdaptiveRejection.h"
#include "StepSizeTuner.h"
#elif SAMPLERTYPE==2
#include "HamiltonianMonteCarlo.h"

typedef struct{
  int dim;
  int n;
  double eps0;
  double eps1;
  const double* sumlogtheta;
}AlphaSamplerArgs;

#endif

class DirichletParamSampler
{
public:
  DirichletParamSampler();
  DirichletParamSampler(unsigned, unsigned, double, double);
  ~DirichletParamSampler();
  
  void SetSize( unsigned, unsigned, double, double );
  void SetPriorEta( double, double );
  void SetPriorMu( const double* const);
  void Sample( unsigned int, const double* const, std::vector<double> *alpha);
  double getEtaStepSize()const;
  double getEtaExpectedAcceptanceRate()const;
  double getMuStepSize()const;
  double getMuExpectedAcceptanceRate()const;
  
private:
#if SAMPLERTYPE==1
  StepSizeTuner TuneEta;
  unsigned int d;
  double eta;
  double *mu;
  double etanew;
  double *munew;
  double *gamma;
  double EtaAlpha;
  double EtaBeta;
  double step, step0;
  double LogAccProb;

  AdaptiveRejection** DirParamArray;
  // AlphaParameters is an array with 5 elements
  // element 0 is number of observations
  // element 1 is dispersion parameter
  // element 2 is 1 - last proportion parameter
  // element 3 is 
  // element 4 is 
  double AlphaParameters[5];

  //MuSampler muSampler;

  void SampleEta(unsigned n, const double* const sumlogtheta, double *eta, const double* const mu);

  static double logf( double, const void* const );
  
  static double dlogf( double, const void* const );
  
  static double ddlogf( double, const void* const );
#elif SAMPLERTYPE==2
  AlphaSamplerArgs AlphaArgs;
  double *logalpha;
  HamiltonianMonteCarlo AlphaSampler;
  double initialAlphaStepsize;
  float targetAlphaAcceptRate;

  static double findE(const double* const theta, const void* const args);
  static void gradE(const double* const theta, const void* const args, double *g);
#endif
};

#endif /* ! DirichletParamSampler_H */
