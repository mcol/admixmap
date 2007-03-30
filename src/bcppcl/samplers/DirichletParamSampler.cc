/** 
 *   DirichletParamSampler.cc 
 *   Class to sample parameters of a Dirichlet distribution
 *   Copyright (c) 2005, 2006 Clive Hoggart, David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "DirichletParamSampler.h"
#include <algorithm>
#include <numeric>
#include "utils/misc.h"
#include <gsl/gsl_math.h>

using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

DirichletParamSampler::DirichletParamSampler() {
  Initialise();
}

// DirichletParamSampler::DirichletParamSampler( unsigned numind, unsigned numpops) {
//   Initialise();
//   SetSize(numind, numpops);
// }

void DirichletParamSampler::Initialise() {
#if DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_ARS_SAMPLER
  mu = 0;
  munew = 0;
  muDirichletParams = 0;
#elif DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_HAMILTONIAN_SAMPLER
  logalpha = 0;
 
#endif
}

/// sets number of elements in Dirichlet parameter vector.
/// Instantiates an adaptive rejection sampler object for each element and 
/// sets up sampler objects. 
///  numobs = number of observations 
void DirichletParamSampler::SetSize( unsigned numobs, unsigned numStates, float InitialStepSize, unsigned NumLeapFrogSteps)
{
   K = numStates;
#if DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_ARS_SAMPLER
   AlphaParameters[0] = numobs;
   muDirichletParams = new double[K];
   mu = new double[K];
   munew = new double[K];
   for( unsigned int i = 0; i < K; i++ )
     muDirichletParams[i] = 1.0;
   DirParamArray.Initialise(true, true, 1.0, 0.0, logf, dlogf); // avoid singularity at mu[j]=0

   EtaArgs.priorshape = K; // for compatibility with gamma(1, 1) prior on alpha
   EtaArgs.priorrate = 1;
   EtaArgs.numpops = K;
   EtaArgs.numobs = numobs; 
   // use many small steps to ensure that leapfrog algorithm does not jump to minus infinity
   EtaSampler.SetDimensions(1, InitialStepSize, 0.0001/*min stepsize*/, 1.0/*max stepsize*/, 
			    NumLeapFrogSteps, 0.9/*target accept rate*/, etaEnergy, etaGradient);
#elif DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_HAMILTONIAN_SAMPLER
   logalpha = new double[K];
   AlphaArgs.n = numobs; //num individuals/gametes will be passed as arg to sampler
   AlphaArgs.dim = K;
   AlphaArgs.eps0 = 1.0; //Gamma(1, 1) prior on alpha
   AlphaArgs.eps1 = 1.0;
   initialAlphaStepsize = 0.05;//need a way of setting this without recompiling, or a sensible fixed value
   targetAlphaAcceptRate = 0.8;// use high value with many leapfrog steps
   AlphaSampler.SetDimensions(K, initialAlphaStepsize, 0.01, 100.0, 40, targetAlphaAcceptRate, alphaEnergy, alphaGradient);
#endif
}

DirichletParamSampler::~DirichletParamSampler()
{
#if DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_ARS_SAMPLER
  delete[] mu;
  delete[] munew;
  delete[] muDirichletParams;
#elif DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_HAMILTONIAN_SAMPLER
  delete[] logalpha;
#endif
}

#if DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_ARS_SAMPLER
void DirichletParamSampler::SetPriorEta( double inEtaAlpha, double inEtaBeta ) {
   EtaArgs.priorshape = inEtaAlpha; 
   EtaArgs.priorrate = inEtaBeta;
}
void DirichletParamSampler::SetPriorMu( const double* const ingamma ) {
   for( unsigned int i = 0; i < K; i++ ){
      muDirichletParams[i] = ingamma[i];
   }
}
#endif

/**
   Samples new values.
   sumlogtheta = sum log observed proportions. 
 
   Updates elements of mu with adaptive rejection sampler conditional on sumlogtheta
   Updates eta with Hamiltonian.  
*/
void DirichletParamSampler::Sample( const double* const sumlogtheta, std::vector<double>& alpha) {
#if DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_ARS_SAMPLER
  //separate alpha into mu and eta
  eta = accumulate(alpha.begin(), alpha.end(), 0.0, std::plus<double>());//eta = sum of alpha[0]
	
  for( unsigned i = 0; i < K; i++ ) {
    mu[i] = alpha[i]/eta;
  }

  //sample proportions
  AlphaParameters[1] = eta; // dispersion parameter
  SampleMu(mu, sumlogtheta);

  //sample dispersion
  eta = SampleEta(eta, mu, sumlogtheta);
  for( unsigned j = 0; j < K; j++ ) alpha[j] = mu[j]*eta;
  
#elif DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_HAMILTONIAN_SAMPLER
  // *** Hamiltonian sampler for alpha
  AlphaArgs.sumlogtheta = sumlogtheta;
  transform(alpha.begin(), alpha.end(), logalpha, xlog);//logalpha = log(alpha)
  AlphaSampler.Sample(logalpha, &AlphaArgs);//sample new values for logalpha
  transform(logalpha, logalpha+K, alpha.begin(), xexp);//alpha = exp(logalpha)
#endif
  
}

#if DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_ARS_SAMPLER
void DirichletParamSampler::SampleMu(double* mu, const double* const sumlogtheta){
  double b = 0.0; 
  for(int updates=0; updates < 2; ++ updates) { // loop twice 
    // loop over elements j,k of mu to update mu[j] conditional on (mu[j] + mu[k]), mu[i] where i neq j,k 
    for( unsigned int j = 1; j < K; ++j ) {
      for( unsigned int k = 0; k < j; ++k ) {
	b = mu[j] + mu[k]; 
	AlphaParameters[2] = b; 
	AlphaParameters[3] = sumlogtheta[j]; 
	AlphaParameters[4] = sumlogtheta[k];
	DirParamArray.setUpperBound(b); // avoid singularity at b
	mu[j] = DirParamArray.Sample(AlphaParameters, ddlogf);
	mu[k] = b - mu[j];
      }
    }
  }
}
//public function
void DirichletParamSampler::SampleEta(const double* const sumlogtheta, std::vector<double>& alpha){
  //separate alpha into mu and eta
  eta = accumulate(alpha.begin(), alpha.end(), 0.0, std::plus<double>());//eta = sum of alpha[0]
	
  for( unsigned i = 0; i < K; i++ ) {
    mu[i] = alpha[i]/eta;
  }

  eta = SampleEta(eta, mu, sumlogtheta);
  for( unsigned j = 0; j < K; j++ ) alpha[j] = mu[j]*eta;

}
void DirichletParamSampler::SampleEta(const double* const sumlogtheta, double* alpha){
  //separate alpha into mu and eta
  eta = accumulate(alpha, alpha+K, 0.0, std::plus<double>());//eta = sum of alpha[0]
	
  for( unsigned i = 0; i < K; i++ ) {
    mu[i] = alpha[i]/eta;
  }

  eta = SampleEta(eta, mu, sumlogtheta);
  for( unsigned j = 0; j < K; j++ ) alpha[j] = mu[j]*eta;

}

//private function
double DirichletParamSampler::SampleEta(double eta, const double* mu,  const double* const sumlogtheta){
  //SampleEta((unsigned)AlphaParameters[0], sumlogtheta, &eta, mu);//first arg is num obs
  EtaArgs.sumlogtheta = sumlogtheta;
  EtaArgs.mu = mu;
  // cout << "eta " << eta << endl << flush;
  etanew = log(eta);//sample for log of dispersion parameter
  EtaSampler.Sample(&etanew, &EtaArgs);
  return exp(etanew);
}

double DirichletParamSampler::etaEnergy( const double* const x, const void* const vargs )
{
  const PopAdmixEtaSamplerArgs* args = (const PopAdmixEtaSamplerArgs* )vargs;
  double eta = exp(*x);
  double E = 0.0;
  const double* mu = args->mu;
  const double* sumlogtheta = args->sumlogtheta;
  try{
    // minus log prior (in log eta basis)
    E +=  args->priorrate *  eta - args->priorshape * log(eta);
    // minus log likelihood
    E -= args->numobs * lngamma(eta);
    for( unsigned i = 0; i < args->numpops; ++i ) {
      E += args->numobs * lngamma(mu[i]*eta) - mu[i]*eta*sumlogtheta[i]; 
    }
  }
  catch(string s){
    throw string("Error in etaEnergy: " + s) ;
  }
  return E;
}

void DirichletParamSampler::etaGradient( const double* const x, const void* const vargs, double* g )
{
  const PopAdmixEtaSamplerArgs* args = (const PopAdmixEtaSamplerArgs* )vargs;
  const double* mu = args->mu;
  const double* sumlogtheta = args->sumlogtheta;
  double eta = exp(*x);
  g[0] = args->priorrate * eta - args->priorshape;  // derivative of minus log prior wrt log eta
  try{ // dE = derivative of minus log likelihood wrt eta
    double dE = -digamma(eta);
    for( unsigned i = 0; i < args->numpops; ++i ) {
      // as mu tends to 0, mu * ( digamma(eta*mu)) tends to -1/eta 
      dE += mu[i] * ( digamma(mu[i]*eta) - sumlogtheta[i] / args->numobs);
    }
    dE *= (double)args->numobs;
    //use chain rule 
    g[0] += eta * dE;
  } catch(string s) {
    throw string("Error in etaGradient: " + s) ;
  }
}

// these 3 functions calculate log density and derivatives for adaptive rejection sampling of 
// a pair of elements of the proportion parameter of the Dirichlet distribution
double DirichletParamSampler::logf( double muj, const void* const pars ) {
  const double* parameters = (const double*) pars;
   int n = (int)parameters[0];
   double eta = parameters[1], b = parameters[2], sumlogpj = parameters[3], sumlogpk = parameters[4];
   if(muj < 0 || (b - muj < 0)) {
     throw string("\nDirichletParamSampler: negative argument to lngamma function\n");
   }

   double f, y1, y2;
   try{
     y1 = lngamma(muj*eta);
     y2 = lngamma((b-muj)*eta);
     f = eta * muj * ( sumlogpj - sumlogpk ) - n * ( y1 + y2 );
   }
   catch(string s){
     throw string("\nERROR in DirichletParamSampler::logf " +s);
   }

   //cout << "\nlog density function passed muj " << muj << "\treturns logdensity " << f << endl;
   return f;
}

double DirichletParamSampler::dlogf( double muj, const void* const pars ) {
  const double* parameters = (const double*) pars;
  double f, x1, x2, y1, y2;
  int n = (int)parameters[0];
  double eta = parameters[1], b = parameters[2], sumlogpj = parameters[3], sumlogpk = parameters[4];
  if(muj < 0 || (b - muj < 0)) {
    throw string("\nDirichletParamSampler: negative argument to digamma function\n");
  } 
  x1 = eta*muj;
  x2 = eta*(b - muj);

  try{
    y1 = digamma(x1);
    y2 = digamma(x2);
    f =  eta * ( sumlogpj - sumlogpk - n*( y1 - y2) );
  }
   catch(string s){
     throw string("\nERROR in DirichletParamSampler::dlogf " +s);
   }
  //cout << "\ngradient function passed muj "<< muj << "\treturns gradient " << f << flush;
  return f;
}

double DirichletParamSampler::ddlogf( double muj, const void* const pars) {
  const double* parameters = (const double*) pars;
  double f, x1, x2;
  int n = (int)parameters[0];
  double eta = parameters[1], b = parameters[2];
  x1 = eta*muj;
  x2 = eta*(b - muj);
  try{
    f = -n * eta * eta *( trigamma(x1) + trigamma(x2) );
    if(f >= 0){
      throw string("DirichletParamSsampler: 2nd derivative non-negative\n");
    }
  }
  catch(string s){
     throw string("\nERROR in DirichletParamSampler::ddlogf " +s);
  }
  return(f);
}

#elif DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_HAMILTONIAN_SAMPLER

///calculate objective function (-log posterior) for log alpha, used in Hamiltonian Metropolis algorithm
double DirichletParamSampler::alphaEnergy(const double* const theta, const void* const vargs){
  /*
    theta = log dirichlet parameters (alpha)
    n = #individuals/gametes
    sumlogtheta (array, length dim) = sums of logs of individual admixture proportions
    eps0, eps1 = parameters of Gamma prior for alpha
  */
  const AlphaSamplerArgs* args = (const AlphaSamplerArgs*)vargs;

  double E = 0.0;
  double sumalpha = 0.0, sumgamma = 0.0, sumtheta = 0.0, sume = 0.0;

  try{
    for(int j = 0; j < args->dim; ++j){
      double alpha = exp(theta[j]);
      sumalpha += alpha;
      sumgamma += lngamma(alpha);
      sume += alpha * (args->eps1 - args->sumlogtheta[j]);
      sumtheta += theta[j];
    }
    E = args->n * (lngamma(sumalpha) - sumgamma) - sume + args->eps0 * sumtheta;
  }

  catch(string s){
    throw string( "Error in DirichletParamSampler::alphaEnergy, " + s) ;
  }
  return -E;

}

///calculate gradient for log alpha
void DirichletParamSampler::alphaGradient(const double* const theta, const void* const vargs, double *g){
  const AlphaSamplerArgs* args = (const AlphaSamplerArgs*)vargs;
  double sumalpha = 0.0, x, y1, y2;
  try{
    for(int j = 0; j < args->dim; ++j) {
      g[j] = 0.0;
      sumalpha += exp(theta[j]);
    }
    y1 = digamma(sumalpha);
    for(int j = 0; j < args->dim; ++j) {
      x = exp(theta[j]);
      y2 = digamma(x);
      if(x > 0.0 && gsl_finite(y1) && gsl_finite(y2)){//to avoid over/underflow problems
	g[j] = x *( args->n *(y2 - y1) + (args->eps1 - args->sumlogtheta[j])) - args->eps0;
      }
    }
  }
  catch(string s){
    throw string( "Error in DirichletParamSampler::alphaGradient, " + s) ;
  }
}
#endif

double DirichletParamSampler::getStepSize()const {
#if DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_ARS_SAMPLER
  //return TuneEta.getStepSize();
  return EtaSampler.getStepsize();
#elif DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_HAMILTONIAN_SAMPLER
    return AlphaSampler.getStepsize();
#endif
}

double DirichletParamSampler::getExpectedAcceptanceRate()const {
#if DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_ARS_SAMPLER
  //return TuneEta.getExpectedAcceptanceRate();
  return EtaSampler.getAcceptanceRate();
#elif DIRICHLETPARAM_SAMPLERTYPE==DIRICHLETPARAM_HAMILTONIAN_SAMPLER
    return AlphaSampler.getAcceptanceRate();
#endif
}

