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
#include "utils/misc.h"
using namespace::std;

unsigned DispersionSampler::K;
unsigned DispersionSampler::L;
unsigned *DispersionSampler::NumStates;


DispersionSampler::DispersionSampler(){
  Args.alpha = 0;
  Args.counts = 0;
}
DispersionSampler::~DispersionSampler(){
  delete[] Args.counts;
  delete[] Args.alpha;
  if(NumStates) delete[] NumStates;
}

/**
   step0, min, max = initial, min and max stepsizes for sampler
   target = target acceptance rate 
   eta is constant over L experiments (loci)
   mu is constant over K experiments (populations) at each locus
   alpha[i] has H[i] elements
   counts[i] has offset H[i]*k + h for count of h th state in k th experiment 
*/

void DispersionSampler::setDimensions(unsigned inL, unsigned inK, int* const inH){
  L = inL;
  K = inK;
  NumStates = new unsigned[L];
  for(unsigned i = 0; i < L; ++i)
    NumStates[i] = inH[i];
}

void DispersionSampler::Initialise(double step0, double min, double max, double target) {
  Args.alpha = new const double*[L];
  Args.counts = new const int*[L]; 
  //set up Hamiltonian sampler
  Sampler.SetDimensions(1, step0, min, max, 20, target, etaEnergyFunction, etaGradient);
  //set default gamma prior for eta
  Args.priorshape = 10.0;
  Args.priorrate = 0.1; 
  logeta[0] = log(Args.priorshape/Args.priorrate); //initialise eta at prior mean
}

void DispersionSampler::setEtaPrior(double shape, double rate){
  Args.priorshape = shape;
  Args.priorrate = rate;
  logeta[0] = log(shape/rate); //initialize at prior mean;
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
double DispersionSampler::etaEnergyFunction(const double* const logeta, const void* const vargs) {
  const EtaSamplerArgs* args = (const EtaSamplerArgs*)vargs;

  double E = 0.0;
  double eta = exp(logeta[0]);
  E += args->priorrate * eta - args->priorshape * logeta[0];  // minus log prior in log eta basis
  for(unsigned i = 0; i < L; ++i) { // loop over elements with same eta
    for(unsigned k = 0; k < K; ++k) { // loop over elements with same alpha parameters

      double scalefactor = 0.0;
      for(unsigned h = 0; h < NumStates[i]; ++h) { // loop over states
	scalefactor += args->alpha[i][h];
      }
      scalefactor = eta / scalefactor;

      double nik = 0.0;//sum of counts over locus i, pop k
      for(unsigned h = 0; h < NumStates[i]; ++h) { // loop over states
	double alpha = args->alpha[i][h] * scalefactor; // rescale alpha to sum to eta
	double count = args->counts[i][K*h + k];
	if(count > 0) {
	  E += lngamma(alpha) - lngamma(count + alpha);
	  nik += count;
	}
      }
      if(nik > 0) {
	E +=  lngamma(eta+nik) - lngamma(eta);
      }
    }
  }
  return E;
}

void DispersionSampler::etaGradient(const double* const logeta, const void* const vargs, double* g) {
  const EtaSamplerArgs* args = (const EtaSamplerArgs*)vargs;
  // alpha[h] = mu[h] * eta
  // d_alpha[h]/d_eta = mu[h] = alpha[h] / eta,  d_eta / d_logeta = eta
  // so we have d_E(alpha[h]])/ d_logeta = d_E / d_alpha[h] * alpha[h] 

  double eta = exp(logeta[0]);
  g[0] =  args->priorrate * eta - args->priorshape; // gradient wrt logeta of minus log prior 
  try {
    for(unsigned i = 0; i < L; ++i) {
      for(unsigned k = 0; k < K; ++k) {
	double nik = 0.0;//sum of counts for locus i, pop k
	for(unsigned h = 0; h < NumStates[i]; ++h) {
	  double alpha = args->alpha[i][h];
	  double count = args->counts[i][K*h + k];
	  nik += count;
	  if(count > 0) {
	    g[0] += alpha * ( digamma(alpha) - digamma(count + alpha) );
	  }
	}
	if(nik > 0) {
	  g[0] += eta * ( digamma(eta + nik) - digamma(eta) );
	}
      }
    }
  } catch(string s) {
    throw string("Error in etaGradient " + s) ;
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
