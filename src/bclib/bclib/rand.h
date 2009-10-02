// *-*-C++-*-*
/* 
 *   rand.h 
 *   Random number class
 *   Copyright (c) 2005, 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef RAND_H
#define RAND_H 
#include "bclib/bclib.h"
#include <iostream>
#include <vector>
extern "C" {
#include <gsl/gsl_rng.h>
}

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */



/// Random number generators
class Rand{
public:
  Rand();
  ~Rand();
  ///set random seed 
  static void setSeed( long seed );
  //standard uniform
  static double myrand();
  ///uniform
  static double myrandRange( double Min, double Max );
  ///gamma
  static double gengam(double aa,double bb);
  ///beta
  static double genbet(double aa,double bb);
  ///univariate normal
  static double gennor(double av,double sd);
  ///binomial
  static int genbinomial( int n, double p );
  ///Poisson
  static unsigned int genpoi( double );
  ///multinomial
  static std::vector<int> genmultinomial( int n, const std::vector<double> p );
  
  ///Poisson generator
  static long ignpoi(double mu);

  static int SampleFromDiscrete( const double probs[] , int numberofelements);
  ///Dirichlet
  template<typename ConstVecType, typename VecType> \
	static void gendirichlet( size_t K, const ConstVecType & alpha, VecType & theta );
  static void gendirichlet( size_t K, const double * alpha, double * theta )
      { gendirichlet<const double *, double *>( K, alpha, theta ); }

private:
  static gsl_rng *RandomNumberGenerator;

  Rand(const Rand&);
  Rand& operator=(const Rand);
};


/** @} */

END_BCLIB_NAMESPACE

#endif
