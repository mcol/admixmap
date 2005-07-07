/** 
 *   ADMIXMAP
 *   HMCMC.cc 
 *   Class to implement a Hamiltonian (or hybrid )MCMC sampler
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
#include "HMCMC.h"

using namespace::std;

HMCMC::HMCMC(){
  n = 100;
  dim = 1;
  epsilon = 0.0;
  Tau = 1;
  eps0 = 1.0; 
  eps1 = 1.0;
  g = 0;
};

HMCMC::~HMCMC(){
  delete[] g;
}

//set dimensions
void HMCMC::SetDimensions(const unsigned pdim, const double pepsilon, const unsigned pTau, 
		   const unsigned pn, const double peps0, const double peps1){
  dim = pdim;
  epsilon = pepsilon;
  Tau = pTau;
  n = pn;
  eps0 = peps0;
  eps1 = peps1;

  g = new double[dim];
}

//calculate objective function
double HMCMC::findE(double *theta, unsigned n, double *sumlogtheta, double eps0, double eps1){
  /*
    theta = log dirichlet parameters (alpha)
    n = #individuals/gametes
    sumlogtheta = sums of logs of individual admixture proportions
    eps0, eps1 = parameters of Gamma prior for alpha
  */

  double E = 0.0;
  double sumalpha = 0.0, sumgamma = 0.0, sumtheta = 0.0, sume = 0.0;
  for(unsigned j = 0; j < dim;++j){
    sumalpha += exp(theta[j]);
    sumgamma += gsl_sf_lngamma(exp(theta[j]));
    sume += exp(theta[j]) * (eps1 - sumlogtheta[j]);
    sumtheta += theta[j];
  }
  
  E = n * (gsl_sf_lngamma(sumalpha) - sumgamma) - sume + (eps0 - 1.0) * sumtheta;
  return -E;
}

//calculate gradient
void HMCMC::gradE(double *theta, unsigned n, double *sumlogtheta, double eps0, double eps1, double *g){
  delete[] g;
  g = new double[dim];
  double sumalpha = 0.0, x, y1, y2;
  for(unsigned j = 0; j < dim; ++j) {
    sumalpha += exp(theta[j]);
  }
    ddigam(&sumalpha, &y1);
    for(unsigned j = 0; j < dim; ++j) {
      x = exp(theta[j]);
      ddigam(&x, &y2);
      g[j] = n * exp(theta[j]) * (y2 - y1) + exp(theta[j]) * (eps1 - sumlogtheta[j]) - eps0 + 1.0;
    }
}

void HMCMC::Initialise(double *x, double *sumlogtheta){
  gradE ( x, n, sumlogtheta, eps0, eps1, g ) ; // set gradient using initial x
  E = findE ( x, n, sumlogtheta, eps0, eps1 ) ;// set objective function too
}

void HMCMC::Sample(double *x, double *sumlogtheta){
  /*
    x = position
    p = momentum
    H = Hamiltonian
    epsilon = step size
    Tau = number of. leapfrog steps
    dim = dimension (= K)
    gradE = gradient function
    findE = objective function = -log density
  */

  bool accept = false;
  double *p = new double[dim];
  double H, Hnew, dH, sumpsq = 0.0;
  double *gnew, *xnew, Enew;
  xnew = new double[dim];
  gnew = new double[dim];
  
  for(unsigned i = 0; i < n; ++i)p[i] = gennor( 0.0, 1.0 ) ; // initial momentum is Normal(0,1)
  for(unsigned i = 0; i < dim; ++i)sumpsq += p[i]*p[i];
  H = 0.5 * sumpsq + E ; // evaluate H(x,p)
  for(unsigned i = 0; i < dim; ++i) {
    xnew[i] = x[i]; 
    gnew[i] = g[i];
  }
  for(unsigned tau = 0; tau < Tau; ++tau){ // make Tau `leapfrog' steps
    for(unsigned i = 0; i < dim; ++i) p[i] = p[i] - epsilon * gnew[i] * 0.5 ; // make half-step in p
    for(unsigned i = 0; i < dim; ++i) xnew[i] = xnew[i] + epsilon * p[i] ; // make step in x
    gradE ( xnew, n, sumlogtheta, eps0, eps1, gnew ) ; // find new gradient
    for(unsigned i = 0; i < dim; ++i) p[i] = p[i] - epsilon * gnew[i] * 0.5 ; // make half-step in p
  }
  Enew = findE ( xnew, n, sumlogtheta, eps0, eps1 ) ; // find new value of H
  Hnew = sumpsq *0.5 + Enew ;
  dH = Hnew - H ; // Decide whether to accept
  if ( dH < 0.0 ) accept = true ;
  else if ( myrand() < exp(-dH) ) accept = true;
  else accept = false ;
  
  if ( accept ){
    for(unsigned i = 0; i < dim; ++i){
      x[i] = xnew[i]; 
      g[i] = gnew[i];
    }
    E = Enew ;
  }
  
  delete[] p;
  delete[] xnew;
  delete[] gnew;  
}



