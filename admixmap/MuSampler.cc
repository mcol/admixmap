/** 
 *   ADMIXMAP
 *   MuSampler.cc 
 *   Class to sample the proportion parameters of a multinomial-Dirichlet distribution
 *   parameterised as \mu_1, ..., \mu_k, \eta, where \eta = \sum \mu
 *   Copyright (c) 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "MuSampler.h"
#include "AdaptiveRejection.h"
#include "functions.h"
#include <numeric>
#include <gsl/gsl_linalg.h>

using namespace::std;

MuSampler::MuSampler(){
  params = 0;
  muArgs = new double*[3];
  for(int i = 0; i < 4; ++i){muArgs[i] = 0;}
}

MuSampler::~MuSampler(){
  delete[] params;
  for(int i = 0; i < 4; ++i){delete[] muArgs[i];}
  delete[] muArgs;
}

void MuSampler::setDimensions(unsigned inK, unsigned inH, double mustep0, double mumin, double mumax, double mutarget){
  //step0, min, max = initial, min and max stepsizes for sampler
  //target = target acceptance rate - maybe remove this and set to fixed value
  K = inK;
  H = inH;

  params = new double[H];
  muArgs[0] = new double[2]; muArgs[0][0] = K;muArgs[0][1] = H;
  muArgs[1] = new double[H*K];//to hold counts
  muArgs[2] = new double[1]; //to hold eta
  muSampler.SetDimensions(H, mustep0, mumin, mumax, 20, mutarget, 
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

    //assign counts
    for(unsigned t = 0; t < K*H; ++t){
      muArgs[1][t] = (double)Counts[t];
    }
    muArgs[2][0] = eta;
 
    muSampler.Sample(params, muArgs);

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
  double leftbound = 0.1;
  MuSamplerArgs MuParameters;

  AdaptiveRejection SampleMu;
  SampleMu.Initialise( true, true, 1.0, 0.0, fMu, dfMu);

  SampleMu.setLowerBound( leftbound );

  MuParameters.eta = eta;
  MuParameters.K = K;
  MuParameters.counts = Counts;
  
  SampleMu.setUpperBound( eta - leftbound );//set upper limit for sampler

  alpha[ 0 ] = SampleMu.Sample(&MuParameters, ddfMu); //first state/allele
  // Last (second) prior frequency parameter is determined by sum of mu's = eta.
  alpha[ 1 ] = eta - alpha[ 0 ];
}

float MuSampler::getAcceptanceRate(){
  return muSampler.getAcceptanceRate();
}
float MuSampler::getStepsize(){
  return muSampler.getStepsize();
}

//**** End public interface *****************
//*******************************************

//******* Energy and gradient functions for Hamiltonian MC sampler

double MuSampler::muEnergyFunction(unsigned , const double * const params, const double* const *args){
  int H = (int)args[0][1]; // Number of haplotypes/alleles
  int K = (int)args[0][0];// Number of populations
  double eta = args[2][0];//dispersion parameter
  //Counts = args[2];//realized allele counts, stored as ints
  double E = 0.0;
  //params in softmax format , of length H, but may not sum to zero

  double mu[H];
  double a[H];
  softmax(H, mu, params);
  inv_softmax(H, mu, a); // normalized parameters; sum to zero

  E -= (double)(H); //logprior

  //compute z, needed for Jacobian below
  double z = 0.0;
  for(int h = 0; h < H; ++h){
    z+= exp(a[h]);
  }

  for(int h = 0; h < H; ++h){
    double alpha = eta * mu[h];
  for(int k = 0; k < K; ++k){
      int offset = h*K +k;
      E += gsl_sf_lngamma(alpha) - gsl_sf_lngamma(args[1][offset] + alpha);//log likelihood
    }
  }
  E -= logJacobian(a, z, H);
  return E;
}

void MuSampler::muGradient(unsigned , const double * const params, const double* const *args, double *g){
  int H = (int)args[0][1]; // Number of haplotypes/alleles
  int K = (int)args[0][0];// Number of populations
  double eta = args[2][0];//dispersion parameter
  //params in softmax format , of length H, but may not sum to zero

  double mu[H];
  double a[H];
  softmax(H, mu, params);
  inv_softmax(H, mu, a); // normalized parameters; sum to zero


  //prior term is zero
  double z = 0.0;
  for(int h = 0; h < H; ++h){
    z+= exp(a[h]);
  }
  double alphaH = eta * mu[H-1];
  double dEdMu[H-1];fill(dEdMu, dEdMu+H-1, 0.0);

  //first compute gradient wrt mu
  for(int h = 0; h < H-1; ++h){
    double alpha = eta * mu[h];

    for(int k = 0; k < K; ++k){
      
      int offset = h*K +k;
      
      dEdMu[h] += eta * (gsl_sf_psi(alpha) - gsl_sf_psi(args[1][offset]+alpha));//log likelihood term
      dEdMu[h] += eta * (gsl_sf_psi(alphaH) - gsl_sf_psi(args[1][offset]+alphaH));
    }
  }
   //now use chain rule to obtain gradient wrt args
  for(int h = 0; h < H-1; ++h){
    g[h] = 0.0;
    for(int i = 0; i < H-1; ++i){
      if(i == h)g[h] += dEdMu[h] * mu[h] * (1.0 - mu[h]);
      else g[h] -= dEdMu[i] * exp(a[i]) * mu[h] * mu[i]; 
    }
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
double MuSampler::fMu( double alpha, const void* const args )
{
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int K = parameters->K;
  double eta = parameters->eta;
  const int *counts = parameters->counts;
  //double mu = alpha / eta;
  double logprior = 0.0;//0.1 * log( mu ) + 0.1 * log( 1 - mu  );//Beta(1.1, 1.1) prior
  double f = logprior - 2 * gsl_sf_lngamma( alpha ) - 2 * gsl_sf_lngamma( eta - alpha );

  for(int k = 0; k < K; ++k){
    f += gsl_sf_lngamma( alpha+counts[k] ) + gsl_sf_lngamma( eta-alpha+counts[K + k] );//state 1 + state 2
  }

  return f;
}

double MuSampler::dfMu( double alpha, const void* const args )
{
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int K = parameters->K;
  double eta = parameters->eta;
  const int *counts = parameters->counts;

  double x, y1, y2;
  double logprior = 0.0;//0.1 / alpha - 0.1 / ( eta - alpha );//Beta(1.1, 1.1) prior
  double f = logprior;

  x = eta - alpha;
  if(alpha < 0)cout<<"\nError in dfMu in MuSampler.cc - arg alpha to ddigam is negative\n"; 
  ddigam( &alpha, &y1 );
  if(x < 0)cout<<"\nError in dfMu in MuSampler.cc - arg (eta-alpha) to ddigam is negative\n"; 
  ddigam( &x, &y2 );
  f += 2 * ( y2 - y1 );

  for(int k = 0; k < K; ++k){
    //first state/allele
    x = alpha + counts[k];
    ddigam( &x, &y2 );
    f += y2;
    //second state/allele
    x = eta - alpha + counts[K+k];
    ddigam( &x, &y2 );
    f -= y2;
  }

  return f;
}

double MuSampler::ddfMu( double alpha, const void* const args )
{
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int K = parameters->K;
  double eta = parameters->eta;
  const int *counts = parameters->counts;

  double x, y1, y2;
  if(alpha > eta) return 0.00001;
  double logprior = 0.0;//-0.1 / (alpha*alpha) - 0.1 / (( eta - alpha ) * ( eta - alpha ) );
  double f = logprior;
  x = eta - alpha;
  trigam( &alpha, &y1 );
  trigam( &x, &y2 );
  f -= 2 * ( y2 + y1 );

  for(int k = 0; k < K; ++k){
    x = alpha + counts[k];
    trigam( &x, &y2 );
    f += y2;
    x = eta - alpha + counts[K+k];
    trigam( &x, &y2 );
    f += y2;
  }
  return f;
}



