//=============================================================================
//
// Copyright (C) 2005, 2006  David O'Donnell and Paul McKeigue
// Portions Copyright (C) 2010  Marco Colombo
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file rand.h
/// Definition of the bclib::Rand class.
//=============================================================================

#ifndef __bclib_rand_h
#define __bclib_rand_h

#include "bclib/bclib.h"
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
  static std::vector<int> genmultinomial(int n, const std::vector<double>& p);
  
  ///Poisson generator
  static long ignpoi(double mu);

  static int SampleFromDiscrete( const double probs[] , int numberofelements);
  static int SampleFromDiscreteFast(const double *probs, int numberofelements);

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

#endif // ! __bclib_rand_h
