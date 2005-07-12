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
#include "gsl/gsl_math.h"

using namespace::std;

HMCMC::HMCMC(){
  dim = 1;
  epsilon = 0.0;
  Tau = 1;
  g = 0;
  accept_count = 0;
  findE = 0;
  gradE = 0;
};

HMCMC::~HMCMC(){
  delete[] g;
}

//set dimensions
void HMCMC::SetDimensions(const unsigned pdim, const double pepsilon, const unsigned pTau, 
			  double (*pfindE)(unsigned d,double *theta, double **args),
			  void (*pgradE)(unsigned d,double *theta, double **args, double *g)){
  dim = pdim;
  epsilon = pepsilon;
  Tau = pTau;
  findE = pfindE;
  gradE = pgradE;

  g = new double[dim];
}

void HMCMC::Sample(double *x, double **args){
  /*
    x = position
    p = momentum
    H = Hamiltonian
    epsilon = step size
    Tau = number of. leapfrog steps
    dim = dimension (= K)
    gradE = gradient function
    findE = objective function = -log density
    args = 2darray of arguments to gradE and findE, 2nd dimension is dim so effectively a vector of args for each scalar x 
  */

  bool accept = false;
  double *p = new double[dim];
  double H, Hnew, dH, sumpsq = 0.0;
  double *gnew, *xnew, Enew;
  xnew = new double[dim];
  gnew = new double[dim];

  gradE (dim, x, args, g ) ; // set gradient using initial x
  E = findE (dim, x, args ) ;// set objective function too
  
  for(unsigned i = 0; i < dim; ++i)p[i] = gennor( 0.0, 1.0 ) ; // initial momentum is Normal(0,1)
  for(unsigned i = 0; i < dim; ++i)sumpsq += p[i]*p[i];
  H = 0.5 * sumpsq + E ; // evaluate H(x,p)
  for(unsigned i = 0; i < dim; ++i) {//reset xnew and gnew
    xnew[i] = x[i]; 
    gnew[i] = g[i];
  }
  for(unsigned tau = 0; tau < Tau; ++tau){ // make Tau `leapfrog' steps
    for(unsigned i = 0; i < dim; ++i) p[i] = p[i] - epsilon * gnew[i] * 0.5 ; // make half-step in p
    for(unsigned i = 0; i < dim; ++i) xnew[i] = xnew[i] + epsilon * p[i] ; // make step in x
    gradE ( dim, xnew, args, gnew ) ; // find new gradient
    for(unsigned i = 0; i < dim; ++i) p[i] = p[i] - epsilon * gnew[i] * 0.5 ; // make half-step in p
  }
  sumpsq = 0.0;
  for(unsigned i = 0; i < dim; ++i){
     sumpsq += p[i]*p[i];
  }

  Enew = findE ( dim, xnew, args ) ; // find new value of H
    if(Enew !=-1.0){
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
	++accept_count;
	E = Enew ;
      }
    }
  delete[] p;
  delete[] xnew;
  delete[] gnew;  
}

float HMCMC::getAcceptanceCount(){
  return (float)accept_count;
}

