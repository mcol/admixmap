/** 
 *   MuSampler.cc 
 *   Class to sample the proportion parameters of a multinomial-Dirichlet distribution
 *   parameterised as \mu_1, ..., \mu_k, \eta, where \eta = \sum \mu
 *   Copyright (c) 2005, 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "MuSampler.h"
#include "AdaptiveRejection.h"
#include "functions.h"
#include <numeric>
#include <gsl/gsl_linalg.h>

using namespace::std;

MuSampler::MuSampler(){
  params = 0;
}

MuSampler::~MuSampler(){
  delete[] params;
}

void MuSampler::setDimensions(unsigned inK, unsigned inH, double mustep0, double mumin, double mumax, double mutarget){
  //step0, min, max = initial, min and max stepsizes for sampler
  //target = target acceptance rate - maybe remove this and set to fixed value
  K = inK;
  H = inH;

  params = new double[H];
  muArgs.H = H;
  muArgs.K = K;

  muSampler.SetDimensions(H, mustep0, mumin, mumax, 30, mutarget, 
			  muEnergyFunction, muGradient);

}

void MuSampler::Sample(double* alpha, double eta, const int* const Counts){
  //alpha is the array (length H) of Dirichlet parameters
  //these are first transformed to logits of proportions
//Counts is an H*K array of counts
  if(H == 2)Sample1D(alpha, eta, Counts);//can sample directly from beta-binomial in one-dimensional case
  else
  {
    
    //transform alphas
    inv_softmax(H, alpha, params);//NOTE: inv_softmax function works with alpha as well as mu

    muArgs.eta = eta;
    muArgs.counts = Counts;
 
    muSampler.Sample(params, &muArgs);

    //set alpha by reversing transformation
    //NOTE: params may not sum to zero but the proportions will still be correct
    softmax(H, alpha, params);
    //alpha now holds proportions so multiply by eta
    for(unsigned t = 0; t < H; ++t){
      alpha[t] *= eta;
    } 
  }
}

void MuSampler::Sample1D(double* alpha, double eta, const int* const Counts )
{
  //Counts has dimension 2 * K
  //alpha has dimension 2
  double leftbound = 0.0001;
  MuSamplerArgs MuParameters;

  AdaptiveRejection SampleMu;
  SampleMu.Initialise( true, true, 1.0, 0.0, fMu, dfMu);

  SampleMu.setLowerBound( leftbound );

  MuParameters.eta = eta;
  MuParameters.K = K;
  MuParameters.counts = Counts;
  
  SampleMu.setUpperBound( 0.9999 );//set upper limit for sampler

  alpha[ 0 ] = eta * SampleMu.Sample(&MuParameters, ddfMu); //first state/allele
  // Last (second) prior frequency parameter is determined by sum of mu's = eta.
  alpha[ 1 ] = eta - alpha[ 0 ];
}

float MuSampler::getAcceptanceRate()const{
  return muSampler.getAcceptanceRate();
}
float MuSampler::getStepsize()const{
  return muSampler.getStepsize();
}

//**** End public interface *****************
//*******************************************

//******* Energy and gradient functions for Hamiltonian MC sampler

double MuSampler::muEnergyFunction(const double * const params, const void* const vargs){
  const MuSamplerArgs* args = (const MuSamplerArgs*)vargs;

  int H = args->H;// Number of haplotypes/alleles
  int K = args->K;// Number of populations
  double eta = args->eta;//dispersion parameter
  double E = 0.0;
  //params in softmax format , of length H, but may not sum to zero

  double* mu = new double[H];
  double* a = new double[H];
  softmax(H, mu, params);
  inv_softmax(H, mu, a); // normalized parameters; sum to zero

  E -= (double)(H); //logprior

  //compute z, needed for Jacobian below
  double z = 0.0;
  for(int h = 0; h < H; ++h){
    z+= exp(a[h]);
  }

  gsl_error_handler_t* old_handler =  gsl_set_error_handler_off();//disable default gsl error handler
  int status  = 0;
  gsl_sf_result psi1, psi2;

  for(int h = 0; h < H; ++h){
    if(status)break;
    double alpha = eta * mu[h];
  for(int k = 0; k < K; ++k){
      int offset = h*K +k;
      status = gsl_sf_lngamma_e(alpha, &psi1);if(status)break;
      status = gsl_sf_lngamma_e(args->counts[offset]+alpha, &psi2);if(status)break;
      E += psi1.val - psi2.val;//log likelihood
    }
  }
  gsl_set_error_handler (old_handler);//restore gsl error handler 
  E -= logJacobian(a, z, H);
  delete[] mu;
  delete[] a;
  if(status){
    throw string("lngamma error in muEnergy");
  }
  return E;
}

void MuSampler::muGradient(const double * const params, const void* const vargs, double *g){
  const MuSamplerArgs* args = (const MuSamplerArgs*)vargs;

  int H = args->H;// Number of haplotypes/alleles
  int K = args->K;// Number of populations
  double eta = args->eta;//dispersion parameter
  //params in softmax format , of length H, but may not sum to zero
  int status  = 0;
  gsl_sf_result psi1, psi2;

  double* mu = new double[H];
  double* a = new double[H];
  softmax(H, mu, params);
  inv_softmax(H, mu, a); // normalized parameters; sum to zero

  //prior term is zero
  double z = 0.0;
  for(int h = 0; h < H; ++h){
    z+= exp(a[h]);
  }
  delete[] a;

  double* dEdMu = new double[H];fill(dEdMu, dEdMu+H, 0.0);

  gsl_error_handler_t* old_handler =  gsl_set_error_handler_off();//disable default gsl error handler
  //first compute gradient wrt mu
  for(int h = 0; h < H; ++h){
    if(status)break;
    double alpha = eta * mu[h];

    for(int k = 0; k < K; ++k){
      int offset = h*K +k;
      status = gsl_sf_psi_e(alpha, &psi1);if(status)break;
      status = gsl_sf_psi_e(args->counts[offset]+alpha, &psi2);if(status)break;
      dEdMu[h] += psi1.val - psi2.val;//log likelihood term
    }
    dEdMu[h] *= eta;
  }
   //now use chain rule to obtain gradient wrt args
  for(int h = 0; h < H; ++h){
    g[h] = 0.0;
    for(int i = 0; i < H; ++i){
      if(i == h)g[h] += dEdMu[h] * mu[h] * (1.0 - mu[h]);
      else g[h] -= dEdMu[i] * exp(a[i]) * mu[h] * mu[i]; 
    }
  }

  delete[] mu;
  delete[] dEdMu;

  gsl_set_error_handler (old_handler);//restore gsl error handler 
  if(status){
    throw string("digamma error in muGradient");
  }
}
//****** Auxiliary function used in Energy function

double MuSampler::logJacobian(const double* a, const double z, unsigned H){
  //computes logJacobian for softmax transformation

  //construct matrix
  gsl_matrix *J = gsl_matrix_calloc(H-1, H-1);
  for(unsigned i = 0; i < H-1; ++i){
    gsl_matrix_set(J,i,i, a[i]*(z-a[i])/(z*z));//diagonal elements
    for(unsigned j = i+1; j < H-1; ++j){
      gsl_matrix_set(J, i, j, -exp(2*a[i]+a[j])/(z*z));//upper triangle
      gsl_matrix_set(J, j, i, -exp(2*a[j]+a[i])/(z*z));//lower triangle
    }
  }
  //LU decomposition
  gsl_permutation *p = gsl_permutation_alloc(H-1);
  gsl_permutation_init(p);
  int signum =1;
  
  //int status = 
  gsl_linalg_LU_decomp ( J , p, &signum);

  gsl_permutation_free(p);
  double logJ = gsl_linalg_LU_lndet(J); 
  gsl_matrix_free(J);
  return logJ; 
}

//******* Logdensity and derivatives for Adaptive rejection sampler
double MuSampler::fMu( double mu, const void* const args )
{
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int K = parameters->K;
  double eta = parameters->eta;
  const int *counts = parameters->counts;
  double alpha = mu * eta;
  double f = 0.0;

  try{
    double logprior = 0.1 * log( mu ) + 0.1 * log( 1 - mu  );//Beta(1.1, 1.1) prior
    f += logprior - K * (lngamma(alpha) + lngamma(eta-alpha));
    
    for(int k = 0; k < K; ++k){
      f += lngamma(alpha+counts[k]) + lngamma(eta - alpha+counts[K + k]);//state 1 + state 2
    }
  }
  catch(string s){
    throw string("Error in fMu " +s);
  }
  return f;
}

double MuSampler::dfMu( double mu, const void* const args )
{
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int K = parameters->K;
  double eta = parameters->eta;
  const int *counts = parameters->counts;
  double alpha = mu * eta;

  double logprior = 0.1 / mu - 0.1 / ( 1.0 - mu );//Beta(1.1, 1.1) prior
  double f = 0.0;

  try{
    f += K * ( digamma(alpha) - digamma(eta-alpha) );
    
    for(int k = 0; k < K; ++k){
      //first state/allele
      f += digamma(alpha+counts[k]);
      //second state/allele
      f -= digamma(eta - alpha+counts[K+k]);
    }
    f *= eta;
    f += logprior;
  }
  catch(string s){
    throw string("Error in dfMu " +s);
  }
  return f;
}

double MuSampler::ddfMu( double mu, const void* const args )
{
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int K = parameters->K;
  double eta = parameters->eta;
  const int *counts = parameters->counts;
  double alpha = mu * eta;

  if(alpha > eta) return 0.00001;// ??
  double logprior = -0.1 / ( mu * mu) - 0.1 / (( 1.0 - mu ) * ( 1.0 - mu ) );
  double f = 0.0;

  try{
    f -= K * ( trigamma(alpha) + trigamma(eta-alpha) );
    
    for(int k = 0; k < K; ++k){
      f += trigamma(alpha+counts[k]);
      f += trigamma(eta-alpha+counts[K+k]);
    }
    f*= eta*eta;
    f += logprior;
  }
  catch(string s){
    throw string("Error in ddfMu " + s);
  }
  return f;
}



