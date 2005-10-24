/** 
 *   ADMIXMAP
 *   HMC.cc 
 *   Class to implement a Hamiltonian (or hybrid )Monte Carlo sampler
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
#include "HamiltonianMonteCarlo.h"
#include "gsl/gsl_math.h"

using namespace::std;

HamiltonianMonteCarlo::HamiltonianMonteCarlo(){
  dim = 1;
  epsilon = 0.0;
  Tau = 1;
  g = 0;
  //overall_accept_count = 0;
  findE = 0;
  gradE = 0;
  //totalsamples = 0;
};

HamiltonianMonteCarlo::~HamiltonianMonteCarlo(){
  delete[] g;
}

//set dimensions
void HamiltonianMonteCarlo::SetDimensions(unsigned pdim, double pepsilon, double min, double max, unsigned pTau, float target,
			  double (*pfindE)(const double* const theta, const void* const args),
			  void (*pgradE)(const double* const theta, const void* const args, double *g)){
  dim = pdim;
  epsilon = pepsilon;
  Tau = pTau;
  findE = pfindE;
  gradE = pgradE;

  g = new double[dim];

  Tuner.SetParameters( epsilon, min, max, target);
}

void HamiltonianMonteCarlo::Sample(double* const x, const void* const args){
  /*
    x = position
    p = momentum
    H = Hamiltonian
    epsilon = step size
    Tau = number of. leapfrog steps
    dim = dimension (= K)
    gradE = gradient function
    findE = objective function = -log density
    args = pointer to object containing arguments to findE and gradE 
  */

  bool accept = false;
  double AccProb;
  double *p = new double[dim];
  double H, Hnew, dH, sumpsq = 0.0;
  double *gnew, *xnew, Enew;
  xnew = new double[dim];
  gnew = new double[dim];

  gradE (x, args, g ) ; // set gradient using initial x
  E = findE (x, args ) ;// set objective function too
  epsilon = Tuner.getStepSize();
  
  for(unsigned i = 0; i < dim; ++i)p[i] = gennor( 0.0, 1.0 ) ; // initial momentum is Normal(0,1)
  for(unsigned i = 0; i < dim; ++i)sumpsq += p[i]*p[i];
  H = 0.5 * sumpsq + E ; // evaluate H(x,p)
  for(unsigned i = 0; i < dim; ++i) {//reset xnew and gnew
    xnew[i] = x[i]; 
    gnew[i] = g[i];
  }
  for(unsigned tau = 0; tau < Tau; ++tau){ // make Tau `leapfrog' steps
    for(unsigned i = 0; i < dim; ++i) p[i] = p[i] - epsilon * gnew[i] * 0.5 ; // make half-step in p
    for(unsigned i = 0; i < dim; ++i) {xnew[i] = xnew[i] + epsilon * p[i] ; // make step in x
      //cout<<x[i]<<" "<<xnew[i]<<" "<<p[i]<<" "<<g[i]<<" "<<gnew[i]<<" "<<endl;
    }
    //cout<<endl<<endl;
    gradE ( xnew, args, gnew ) ; // find new gradient
    for(unsigned i = 0; i < dim; ++i) p[i] = p[i] - epsilon * gnew[i] * 0.5 ; // make half-step in p
  }
  //cout<<endl;
  sumpsq = 0.0;
  for(unsigned i = 0; i < dim; ++i){
     sumpsq += p[i]*p[i];
  }

  Enew = findE ( xnew, args ) ; // find new value of H

  AccProb = 0.0;
  if(Enew !=-1.0){// -1 means an error in calculation of energy function
      Hnew = sumpsq *0.5 + Enew ;
      dH = Hnew - H ; // Decide whether to accept
      if ( dH < 0.0 ) {accept = true ;AccProb = 1.0;}
      else {
	AccProb = exp(-dH);
	if ( myrand() < exp(-dH) ) accept = true;
	else accept = false ;
      }
      
      if ( accept ){
	for(unsigned i = 0; i < dim; ++i){
	  x[i] = xnew[i]; 
	  g[i] = gnew[i];
	}
	//++overall_accept_count;
	E = Enew ;
      }
  }
  Tuner.UpdateStepSize(AccProb);
  //++totalsamples;
  delete[] p;
  delete[] xnew;
  delete[] gnew;  
}

float HamiltonianMonteCarlo::getAcceptanceRate()const{
  //return (float)overall_accept_count/(float)totalsamples;
  return Tuner.getExpectedAcceptanceRate();
}

float HamiltonianMonteCarlo::getStepsize()const{
  return epsilon;
  //return Tuner.getStepsize();
}


