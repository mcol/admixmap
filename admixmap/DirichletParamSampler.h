// *-*-C++-*-*
#ifndef DirichletParamSampler_H
#define DirichletParamSampler_H 1

#define SAMPLERTYPE 1 // 1 = ARS + MHRW
                     // 2 = Hamiltonian 

#include <gsl/gsl_sf_gamma.h>
#include "rand.h"
//#include "StepSizeTuner.h"
#include "HamiltonianMonteCarlo.h"

#if SAMPLERTYPE==1
#include "AdaptiveRejection.h"
#endif

typedef struct{
  int dim;
  int n;
  double eps0;
  double eps1;
  const double* sumlogtheta;
}AlphaSamplerArgs;

typedef struct{
  unsigned numpops;
  unsigned numobs;
  const double* mu;
  const double* sumlogtheta;
  double priorshape;
  double priorrate;
} PopAdmixEtaSamplerArgs;

class DirichletParamSampler
{
public:
  DirichletParamSampler();
  DirichletParamSampler(unsigned, unsigned);
  ~DirichletParamSampler();
  
  void SetSize( unsigned, unsigned);
  void SetPriorEta( double, double );
  void SetPriorMu( const double* const);
  void Sample( const double* const, std::vector<double> *alpha);
  double getStepSize()const;
  double getExpectedAcceptanceRate()const;
    
private:
  unsigned int K;
#if SAMPLERTYPE==1
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

  static double etaEnergy( const double* const eta, const void* const vargs );
  static void etaGradient( const double* const eta, const void* const vargs, double* g );
#elif SAMPLERTYPE==2
  AlphaSamplerArgs AlphaArgs;
  double *logalpha;
  HamiltonianMonteCarlo AlphaSampler;
  double initialAlphaStepsize;
  float targetAlphaAcceptRate;

  static double findE(const double* const theta, const void* const args);
  static void gradE(const double* const theta, const void* const args, double *g);

#endif

  void Initialise();
};

#endif /* ! DirichletParamSampler_H */
