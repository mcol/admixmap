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
#include "functions.h"
using namespace::std;

MuSampler::MuSampler(){
  params = 0;
  muArgs = new double*[4];
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
  muArgs[1] = new double[1]; 
  muArgs[2] = new double[H*K];//to hold counts
  muArgs[3] = new double[1]; //to hold eta
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
    double sum = 0.0;
    for(unsigned t = 0; t < H-1; ++t){
      sum += log( alpha[t]/eta );
    }
    params[H-1] = -sum / H;
    for(unsigned t = 0; t < H-1; ++t){
      params[t] = params[H-1] + log(alpha[t]/eta) - log(alpha[H-1]/eta);
    }

    //assign counts
    for(unsigned t = 0; t < K*H; ++t){
      muArgs[2][t] = (double)Counts[t];
    }
    muArgs[3][0] = eta;
    muArgs[1][0] = 0.01 * muSampler.getStepsize();//"small" differential of 10% of current stepsize
    //used to compute derivative of Jacobian in gradient
 
    muSampler.Sample(params, muArgs);

    //set alpha by reversing transformation
    double z = 0.0;
    for(unsigned t = 0; t < H; ++t){
      z += exp(params[t]);
    } 
    for(unsigned t = 0; t < H; ++t){
      alpha[t] = eta * exp(params[t])/ z;
    } 

  }
}

void MuSampler::Sample1D(double* alpha, double eta, const int* const Counts )
{
  //Counts has dimension 2 * K
  //alpha has dimension 2
  double lefttruncation = 0.1;
  double MuParameters[2];

  DARS SampleMu( 0, 0, 0, MuParameters, fMu, dfMu, ddfMu, Counts, 0 );

  SampleMu.SetLeftTruncation( lefttruncation );

  MuParameters[0] = eta;
  MuParameters[1] = K;
  
  SampleMu.SetRightTruncation( eta - lefttruncation );//set upper limit for sampler
  SampleMu.UpdateParameters( MuParameters);
  alpha[ 0 ] = SampleMu.Sample(); //first state/allele
  // Last (second) prior frequency parameter is determined by sum of mu's = eta.
  alpha[ 1 ] = eta - alpha[ 0 ];
}

float MuSampler::getAcceptanceRate(){
  return muSampler.getAcceptanceRate();
}
float MuSampler::getStepsize(){
  return muSampler.getStepsize();
}

double MuSampler::muEnergyFunction(unsigned , const double * const a, const double* const *args){
  int H = (int)args[0][1]; // Number of haplotypes/alleles
  int K = (int)args[0][0];// Number of populations
  double eta = args[3][0];//dispersion parameter
  //Counts = args[2];//realized allele counts, stored as ints
  double E = 0.0;
  //alpha_jk = mu_jk *eta

  E -= (double)(H); //logprior

  double z = 0.0;
  for(int h = 0; h < H; ++h){
    z+= exp(a[h]);
  }

  for(int h = 0; h < H; ++h){
    double alpha = eta * exp(a[h])/z;
  for(int k = 0; k < K; ++k){
      int offset = h*K +k;
      E += gsl_sf_lngamma(alpha) - gsl_sf_lngamma(args[2][offset] + alpha);//log likelihood
    }
  }
  E -= logJacobian(a, z, H);
  return E;
}

void MuSampler::muGradient(unsigned , const double * const a, const double* const *args, double *g){
  int H = (int)args[0][1]; // Number of haplotypes/alleles
  int K = (int)args[0][0];// Number of populations
  double eta = args[3][0];//dispersion parameter

  double z = 0.0;
  for(int h = 0; h < H; ++h){
    z+= exp(a[h]);
  }
  double delta = args[1][0];
  double alphaH = eta * exp(a[H])/z;;

  for(int h = 0; h < H; ++h){
    double alpha = eta * exp(a[h])/z;

     for(int k = 0; k < K; ++k){

      int offset = h*K +k;

      g[h] = eta * (gsl_sf_psi(alpha) - gsl_sf_psi(args[2][offset]+alpha));//log likelihood term
      g[h] += eta * (gsl_sf_psi(alphaH) - gsl_sf_psi(args[2][offset]+alphaH));
    }
     g[h] -= DlogJacobian(a, z, H, h, delta);  //Jacobian term
  }
}

double MuSampler::fMu( const double* parameters, const int *counts,  const double *, double alpha )
{
  int K = (int)parameters[1];
  double eta = parameters[0];
  //double mu = alpha / eta;
  double logprior = 0.0;//0.1 * log( mu ) + 0.1 * log( 1 - mu  );//Beta(1.1, 1.1) prior
  double f = logprior - 2 * gsl_sf_lngamma( alpha ) - 2 * gsl_sf_lngamma( eta - alpha );

  for(int k = 0; k < K; ++k){
    f += gsl_sf_lngamma( alpha+counts[k] ) + gsl_sf_lngamma( eta-alpha+counts[K + k] );//state 1 + state 2
  }

  return f;
}

double MuSampler::dfMu( const double* parameters, const int *counts, const double *, double alpha )
{
  int K = (int)parameters[1];
  double eta = parameters[0], x, y1, y2;
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

double MuSampler::ddfMu( const double* parameters, const int *counts, const double *, double alpha )
{
  int K = (int)parameters[1];
  double eta = parameters[0], x, y1, y2;
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
  
  int status = gsl_linalg_LU_decomp ( J , p, &signum);

  gsl_permutation_free(p);
  double logJ = gsl_linalg_LU_lndet(J); 
  gsl_matrix_free(J);
  return logJ; 
}

double MuSampler::DlogJacobian(const double* const a, const double z, unsigned H, unsigned h, double delta){
  //computes numerical approximation to derivative of above Jacobian wrt a[h]
  //room for improvement here ?

  double b[H];
  copy(a, a+H, b);
  b[h] += delta;
  double d = logJacobian(b, z, H);
  d -= logJacobian(a, z, H);
  return d / delta;
}
