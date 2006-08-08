/** 
 *   HamiltonianMonteCarlo.cc 
 *   Class to implement a Hamiltonian (or hybrid )Monte Carlo sampler
 *   Copyright (c) 2005, 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "HamiltonianMonteCarlo.h"
#include "gsl/gsl_math.h"
#include <sstream>//for error handling

using namespace::std;

HamiltonianMonteCarlo::HamiltonianMonteCarlo(){
  dim = 1;
  epsilon = 0.0;
  Tau = 1;
  //overall_accept_count = 0;
  findE = 0;
  gradE = 0;
  //totalsamples = 0;
  xnew = 0;
  g = 0;
  gnew = 0;
  p = 0;
}

HamiltonianMonteCarlo::~HamiltonianMonteCarlo(){
  delete[] p;
  delete[] xnew;
  delete[] gnew;
  delete[] g;
}

///set dimensions
void HamiltonianMonteCarlo::SetDimensions(unsigned pdim, double pepsilon, double min, double max, unsigned pTau, float target,
			  double (*pfindE)(const double* const theta, const void* const args),
			  void (*pgradE)(const double* const theta, const void* const args, double *g)){
  dim = pdim;
  epsilon = pepsilon;
  Tau = pTau;
  findE = pfindE;
  gradE = pgradE;
  g = new double[dim];        //gradient (multidim)
  p = new double[dim];
  xnew = new double[dim];
  gnew = new double[dim];

  Tuner.SetParameters( epsilon, min, max, target);
}

void HamiltonianMonteCarlo::Sample(double* const x, const void* const args){
  /*
    x = position
    p = momentum
    H = Hamiltonian
    epsilon = step size
    Tau = number of leapfrog steps
    dim = dimension (= K)
    gradE = gradient function
    findE = objective function = -log density
    args = pointer to object containing arguments to findE and gradE 
  */
  double E = 0.0;         //value of objective function
  bool accept = false;
  double AccProb;
  double H, Hnew, dH, sumpsq = 0.0;
  double Enew;

  try{
    gradE (x, args, g ) ; // set gradient using initial x
    E = findE (x, args ) ;// set objective function too
    epsilon = Tuner.getStepSize();
  }
  catch(string s){
    std::ostringstream error_string;
    error_string << "Error in HamiltonianSampler:\n" << s
		 << " \nProbably passing bad arguments: x = ";
    for(unsigned i = 0; i < dim; ++i)error_string << x[i] <<" ";
  }
  try{
    for(unsigned i = 0; i < dim; ++i)p[i] = Rand::gennor( 0.0, 1.0 ) ; // initial momentum is Normal(0,1)
    for(unsigned i = 0; i < dim; ++i)sumpsq += p[i]*p[i];
    H = 0.5 * sumpsq + E ; // evaluate H(x,p)
    for(unsigned i = 0; i < dim; ++i) {//reset xnew and gnew
      xnew[i] = x[i]; 
      gnew[i] = g[i];
    }
    for(unsigned tau = 0; tau < Tau; ++tau){ // make Tau `leapfrog' steps
      for(unsigned i = 0; i < dim; ++i) p[i] = p[i] - epsilon * gnew[i] * 0.5 ; // make half-step in p
      for(unsigned i = 0; i < dim; ++i) {
	xnew[i] = xnew[i] + epsilon * p[i] ; // make step in x
	if( !gsl_finite(xnew[i]) ) {
	  throw string("\nleapfrog to infinity in Hamiltonian sampler - try using more small steps");
	}
      }
      gradE ( xnew, args, gnew ) ; // find new gradient
      for(unsigned i = 0; i < dim; ++i) {p[i] = p[i] - epsilon * gnew[i] * 0.5 ; // make half-step in p
	//cout<<x[i]<<" "<<xnew[i]<<" "<<p[i]<<" "<<g[i]<<" "<<gnew[i]<<" "<<endl;
      }
    }

    sumpsq = 0.0;
    for(unsigned i = 0; i < dim; ++i){
      sumpsq += p[i]*p[i];
    }
    
    Enew = findE ( xnew, args ) ; // find new value of H
  }
  catch(string s){
    std::ostringstream error_string;
    error_string << "Error in HamiltonianSampler:\n" << s
		 << " \nPassing argument x = ";
    for(unsigned i = 0; i < dim; ++i)error_string << xnew[i] <<" ";
    error_string << "\nCurrent stepsize = " << epsilon
		 << "\nCurrent momentum = ";
    for(unsigned i = 0; i < dim; ++i)error_string << p[i] <<" ";
    error_string << "\n and gradient = ";
    for(unsigned i = 0; i < dim; ++i)error_string << gnew[i] <<" ";

    throw(error_string.str());
  }

  AccProb = 0.0;
  if(Enew !=-1.0){// -1 means an error in calculation of energy function
      Hnew = sumpsq *0.5 + Enew ;
      dH = Hnew - H ; // Decide whether to accept
      //cout << dH << " " << E << " " << Enew << " ";
      if ( dH < 0.0 ) {accept = true ;AccProb = 1.0;}
      else {
	AccProb = exp(-dH);
	if ( Rand::myrand() < exp(-dH) ) accept = true;
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
  //cout << endl;

  Tuner.UpdateStepSize(AccProb);
  //++totalsamples;
  
}

float HamiltonianMonteCarlo::getAcceptanceRate()const{
  //return (float)overall_accept_count/(float)totalsamples;
  return Tuner.getExpectedAcceptanceRate();
}

float HamiltonianMonteCarlo::getStepsize()const{
  return epsilon;
  //return Tuner.getStepsize();
}

void HamiltonianMonteCarlo::resetStepSizeApproximator(int k) {
  Tuner.resetStepSizeApproximator(k);
}

