/** 
 *   DispersionSampler.cc 
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

#include "DispersionSampler.h"
#include "functions.h"
using namespace::std;

DispersionSampler::DispersionSampler(){
  Args.alpha = 0;
  Args.counts = 0;
}
DispersionSampler::~DispersionSampler(){
  delete[] Args.counts;
  delete[] Args.alpha;
}


void DispersionSampler::setDimensions(unsigned inL, unsigned inK, int* const inH, double step0, double min, double max, double target){
  //step0, min, max = initial, min and max stepsizes for sampler
  //target = target acceptance rate - maybe remove this and set to fixed value
  L = inL;
  K = inK;

  Args.L = L;
  Args.K = K;
  Args.H = inH;
  Args.alpha = new const double*[L];
  Args.counts = new const int*[L];

  //set up Hamiltonian sampler
  Sampler.SetDimensions(1, step0, min, max, 20, target, etaEnergyFunction, etaGradient);
  //set default gamma priors for eta
  Args.priorshape = 3.0;
  Args.priorrate = 0.01;//mean 300, variance 30 000
  logeta[0] = log(Args.priorshape/Args.priorrate); //initialise eta at prior mean
}

void DispersionSampler::setEtaPrior(double shape, double rate){
  Args.priorshape = shape;
  Args.priorrate = rate;
  logeta[0] = log(shape/rate);//log(shape/rate);
}

void DispersionSampler::addAlphas(unsigned i, const double* const alpha){
  Args.alpha[i] = alpha;
}
void DispersionSampler::addCounts(unsigned i, const int* const counts){
  Args.counts[i] = counts;
}

double DispersionSampler::Sample(){
  Sampler.Sample(logeta, &Args);

  return exp(logeta[0]);
}

float DispersionSampler::getAcceptanceRate()const{
  return Sampler.getAcceptanceRate();
}
float DispersionSampler::getStepsize()const{
  return Sampler.getStepsize();
}
double DispersionSampler::etaEnergyFunction(const double * const logeta, const void* const vargs){
  const EtaSamplerArgs* args = (const EtaSamplerArgs*)vargs;

  int K = args->K;// Number of populations
  int L = args->L;//number of loci

  double E = 0.0;
  double eta = exp(logeta[0]);

  E += (args->priorshape - 1.0)*logeta[0] - args->priorrate * eta;//log (gamma)prior
  E += logeta[0];//Jacobian

  for(int i = 0; i < L; ++i){
 
    for(int k = 0; k < K; ++k){

      double nik = 0.0;//sum of counts for locus i, pop k

       for(int h = 0; h < args->H[i]; ++h){

	double alpha = args->alpha[i][h];
	double count = args->counts[i][h*K +k];

	nik += count;
	E += gsl_sf_lngamma(count + alpha) - gsl_sf_lngamma(alpha);
      }
 
       E += gsl_sf_lngamma(eta) - gsl_sf_lngamma(eta+nik);
    }
  }
  return -E;
}
void DispersionSampler::etaGradient(const double * const logeta, const void* const vargs, double* g){
  const EtaSamplerArgs* args = (const EtaSamplerArgs*)vargs;

  int K = args->K;// Number of populations
  int L = args->L;//number of loci
  g[0] = 0.0;
  double eta = exp(logeta[0]);

  g[0] -= args->priorshape - args->priorrate*eta;//log prior term

  for(int i = 0; i < L; ++i){
    for(int k = 0; k < K; ++k){
      double nik = 0.0;//sum of counts for locus i, pop k
      for(int h = 0; h < args->H[i]; ++h){

	double alpha = args->alpha[i][h];
	double count = args->counts[i][h*K +k];
	nik += count;
	g[0] -= (alpha/eta)*( gsl_sf_psi(count + alpha) - gsl_sf_psi(alpha) );
      }
      g[0] -= gsl_sf_psi(eta) - gsl_sf_psi(eta+nik);
    }
  }
}
double DispersionSampler::getEnergy(double x){
  logeta[0] = x;
  return etaEnergyFunction(logeta, &Args);
}
double DispersionSampler::getGradient(double x){
  double g[1];
  logeta[0] = x;
  etaGradient(logeta, &Args, g);
  return g[0];
}
