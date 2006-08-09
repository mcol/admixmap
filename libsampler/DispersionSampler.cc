/** 
 *   DispersionSampler.cc 
 *   Class to sample the precision parameter eta of a multinomial-Dirichlet distribution
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
#include "misc.h"
using namespace::std;

DispersionSampler::DispersionSampler(){
  Args.alpha = 0;
  Args.counts = 0;
}
DispersionSampler::~DispersionSampler(){
  delete[] Args.counts;
  delete[] Args.alpha;
}


void DispersionSampler::setDimensions(unsigned inL, unsigned inK, int* const inH, double step0, 
				      double min, double max, double target){
  //step0, min, max = initial, min and max stepsizes for sampler
  //target = target acceptance rate 
  // eta is constant over L experiments (loci)
  // mu is constant over K experiments (populations) at each locus 
  L = inL;
  K = inK;
  Args.L = L; 
  Args.K = K; 
  Args.H = inH; // number of states at each locus
  Args.alpha = new const double*[L];
  Args.counts = new const int*[L];
  //set up Hamiltonian sampler
  Sampler.SetDimensions(1, step0, min, max, 20, target, etaEnergyFunction, etaGradient);
  //set default gamma priors for eta
  Args.priorshape = 1;
  Args.priorrate = 1; 
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
double DispersionSampler::etaEnergyFunction(const double* const logeta, const void* const vargs){
  const EtaSamplerArgs* args = (const EtaSamplerArgs*)vargs;
  int K = args->K; 
  int L = args->L; 
  double E = 0.0;
  double eta = exp(logeta[0]);
  E += args->priorshape * logeta[0] - args->priorrate * eta; // log prior in log eta basis
  for(int i = 0; i < L; ++i) {
    for(int k = 0; k < K; ++k) {
      double nik = 0.0;//sum of counts over locus i, pop k
      for(int h = 0; h < args->H[i]; ++h) {
	double alpha = args->alpha[i][h];
	double count = args->counts[i][h*K + k];
	if(count > 0) {
	  E += lngamma(count + alpha) - lngamma(alpha);
	  nik += count;
	}
      }
      if(nik > 0) {
	E += lngamma(eta) - lngamma(eta+nik);
      }
    }
  }
  return -E;
}

void DispersionSampler::etaGradient(const double* const logeta, const void* const vargs, double* g){
  const EtaSamplerArgs* args = (const EtaSamplerArgs*)vargs;
  // alpha[h] = mu[h] * eta
  // d_alpha[h]/d_eta = mu[h] = alpha[h] / eta,  d_eta / d_logeta = eta
  // so we have d_E[h]/ d_logeta = d_E / d_alpha[h] * alpha[h] 
  int K = args->K;// Number of populations
  int L = args->L;//number of loci
  double eta = exp(logeta[0]);
  g[0] =  args->priorrate * eta - -args->priorshape; // gradient wrt logeta of minus log prior 
  for(int i = 0; i < L; ++i){
    for(int k = 0; k < K; ++k){
      double nik = 0.0;//sum of counts for locus i, pop k
      for(int h = 0; h < args->H[i]; ++h){
	double alpha = args->alpha[i][h];
	double count = args->counts[i][h*K +k];
	nik += count;
	if(count > 0) {
	  g[0] += alpha * ( digamma(alpha) - digamma(count + alpha) );
	} 
      }
      if(nik > 0) {
	g[0] += digamma(eta + nik) - digamma(eta);
      }
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
