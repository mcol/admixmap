/** 
 *   ADMIXMAP
 *   DispersionSampler.cc 
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

#include "DispersionSampler.h"
#include "functions.h"
using namespace::std;

DispersionSampler::DispersionSampler(){
  Args = new double*[5];
  Args[0] = new double[2];
  Args[4] = new double[2];//to hold prior params

}
DispersionSampler::~DispersionSampler(){
  for(int i = 0; i < 5; ++i)delete[] Args[i];
  delete[] Args;
}


void DispersionSampler::setDimensions(unsigned inL, unsigned inK, int* const inH, double step0, double min, double max, double target){
  //step0, min, max = initial, min and max stepsizes for sampler
  //target = target acceptance rate - maybe remove this and set to fixed value
  L = inL;
  K = inK;

  unsigned H = 0;
  Args[1] = new double[L]; 
  for(unsigned i = 0; i < L; ++i){
    Args[1][i] = (double)inH[i];
    H += inH[i];//accumulate total number of states/haplotypes
  }
  Args[2] = new double[H]; //to hold alpha params
  Args[0][0] = (double)K; Args[0][1] = (double)L;

  Args[3] = new double[H*K];//to hold counts

  //set up Hamiltonian sampler
  Sampler.SetDimensions(1, step0, min, max, 20, target, etaEnergyFunction, etaGradient);
  //set default gamma priors for eta
  Args[4][0] = 3.0;
  Args[4][1] = 0.01;//mean 300, variance 30 000
  logeta[0] = log(Args[4][0]/Args[4][1]); //initialise eta at prior mean
}

void DispersionSampler::setEtaPrior(double shape, double rate){
  Args[4][0] = shape;
  Args[4][1] = rate;
  logeta[0] = log(10.0);//log(shape/rate);
}

void DispersionSampler::addAlphas(unsigned i, const double* const alpha){
  //count how many have been added so far
  unsigned H = 0;
  if(i > 0)
    for(unsigned j = 0; j < i; ++j){
      H += (int)Args[1][j];
    }

  for(unsigned h = 0; h < Args[1][i]; ++h)
    Args[2][H+h] = alpha[h];
}
void DispersionSampler::addCounts(unsigned i, const int* const counts){
  //count how many have been added so far
  unsigned H = 0;
  if(i > 0)
    for(unsigned j = 0; j < i; ++j){
      H += (int)Args[1][j];
    }
  for(unsigned h = 0; h < K*Args[1][i]; ++h)
      Args[3][H*K + h] = (double)counts[h];
}

double DispersionSampler::Sample(){
  Sampler.Sample(logeta, Args);

  return exp(logeta[0]);
}

float DispersionSampler::getAcceptanceRate(){
  return Sampler.getAcceptanceRate();
}
float DispersionSampler::getStepsize(){
  return Sampler.getStepsize();
}
double DispersionSampler::etaEnergyFunction(unsigned , const double * const logeta, const double* const *args){
  int K = (int)args[0][0];// Number of populations
  int L = (int)args[0][1];//number of loci
  //args[1] are numbers of alleles/haplotypes/states
  //args[3] are counts, dimension L*H*K
  //args[2] are proportions, mu, dimension L*H
  double E = 0.0;
  double eta = exp(logeta[0]);
  double priorshape = args[4][0];
  double priorrate = args[4][1];

  E += (priorshape-1.0)*logeta[0] - priorrate*eta;//log (gamma)prior
  E += logeta[0];//Jacobian
  int H = 0;// counts how many states visited so far
  for(int i = 0; i < L; ++i){
 
    for(int k = 0; k < K; ++k){

      double nik = 0.0;//sum of counts for locus i, pop k

       for(int h = 0; h < args[1][i]; ++h){

	double alpha = args[2][H+h];
	double count = args[3][H*K + h*K +k];

	nik += count;
	E += gsl_sf_lngamma(count + alpha) - gsl_sf_lngamma(alpha);
      }
 
       E += gsl_sf_lngamma(eta) - gsl_sf_lngamma(eta+nik);
    }
    H += (int)args[1][i];
  }
  return -E;
}
void DispersionSampler::etaGradient(unsigned , const double * const logeta, const double* const *args, double* g){
  int K = (int)args[0][0];// Number of populations
  int L = (int)args[0][1];//number of loci
  //args[1] are numbers of alleles/haplotypes/states
  //args[3] are counts, dimension L*H*K
  //args[2] are proportions, mu, dimension L*H
  g[0] = 0.0;
  double eta = exp(logeta[0]);
  double priorshape = args[4][0];
  double priorrate = args[4][1];

  g[0] -= priorshape - priorrate*eta;//log prior term
  int H = 0;
  for(int i = 0; i < L; ++i){
    for(int k = 0; k < K; ++k){
      double nik = 0.0;//sum of counts for locus i, pop k
      for(int h = 0; h < args[1][i]; ++h){

	double alpha = args[2][H+h];
	double count = args[3][H*K + h*K +k];
	nik += count;
	g[0] -= (alpha/eta)*( gsl_sf_psi(count + alpha) - gsl_sf_psi(alpha) );
      }
      g[0] -= gsl_sf_psi(eta) - gsl_sf_psi(eta+nik);
    }
    H += (int)args[1][i];
  }
}
double DispersionSampler::getEnergy(double x){
  logeta[0] = x;
  return etaEnergyFunction(0, logeta, Args);
}
double DispersionSampler::getGradient(double x){
  double g[1];
  logeta[0] = x;
  etaGradient(0, logeta, Args, g);
  return g[0];
}
