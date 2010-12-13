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
/// \file rand.cc
/// Implementation of the bclib::Rand class.
//=============================================================================

#include <cmath>
#include <iostream>
#include <vector>
#include "bclib/rand.h"
#include <limits>
#include <sstream>
extern "C" {
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#ifdef PARALLEL
#include "gsl-sprng.h"
#endif
}


#define INEXPLICABLY_THROW_AWAY_ERROR_INFO  0

#if 0
  #define OUTPUT_RNG(X) X
#else
  #define OUTPUT_RNG(X)
#endif


// For template instantiation:
#include "bclib/pvector.h"


using namespace std;

BEGIN_BCLIB_NAMESPACE

#ifdef PARALLEL
 // allocate a sprng generator
const gsl_rng_type *gsl_rng_sprng20 = &sprng_type;
gsl_rng *Rand::RandomNumberGenerator = gsl_rng_alloc( gsl_rng_sprng20 );

#else
//allocate a Tausworthe generator
gsl_rng *Rand::RandomNumberGenerator = gsl_rng_alloc( gsl_rng_taus );
#endif

Rand::Rand(){

}
Rand::~Rand(){
  gsl_rng_free(RandomNumberGenerator);
}

///standard uniform number generator
double Rand::myrand()
{
  double rv = gsl_rng_uniform( RandomNumberGenerator );
  OUTPUT_RNG(
    fprintf(stderr, "*** myrand: %.9lf\n", rv);
  )
  return rv;
  //?? should use gsl_rng_uniform_pos to exclude 0
}

/// ** U(Min, Max) 
double Rand::myrandRange( double Min, double Max )
{
  double Range = Max - Min;
  return( Min + Range * gsl_rng_uniform( RandomNumberGenerator ) );
}

///set random number seed
void Rand::setSeed( long seed )
{
#ifdef PARALLEL
    init_sprng(DEFAULT_RNG_TYPE, seed, SPRNG_DEFAULT);
#else
  gsl_rng_set(RandomNumberGenerator,
	      static_cast< unsigned long int >( seed ) );
#endif
  //in sprng20 case, calls init_sprng
}

/// ** Gamma distribution **
double Rand::gengam( double shape, double rate )
{
  double x = 0.0;
  do
    x = gsl_ran_gamma( RandomNumberGenerator, shape, 1.0 / rate ) ;
  while (x < numeric_limits<double>::min( ));
  OUTPUT_RNG(
    fprintf(stderr, "*** gengam(%.3lf, %.3lf): %.9lf\n", shape, rate, x);
  )
  return x;
}
/// ** Beta distribution **
double Rand::genbet( double aa, double bb )
{
  return( gsl_ran_beta( RandomNumberGenerator, aa, bb ) );
}
/// ** (univariate) Normal distribution **
double Rand::gennor( double av, double sd )
{
  //check for valid args
  if(isinf(av) || isnan(av) || !(sd > 0.0)){
    stringstream err;
    err <<"Bad arguments to gennor: "<< av << ", " << sd;
    throw err.str();
  }
  double rv = av + gsl_ran_gaussian( RandomNumberGenerator, sd );
  OUTPUT_RNG(
    fprintf(stderr, "*** gennor(%.3lf, %.3lf): %.9lf\n", av, sd, rv);
  )
  return rv;
}
/// ** Binomial distribution **
int Rand::genbinomial( int n, double p )
{
  return( gsl_ran_binomial( RandomNumberGenerator, p, n ) );
}
/// ** Poisson distribution **
unsigned int Rand::genpoi( double mu )
{
  return( gsl_ran_poisson( RandomNumberGenerator, mu ) );
}

long Rand::ignpoi( double mu )
{
   return( gsl_ran_poisson( RandomNumberGenerator, mu ) );
}

/// ** Multinomial distribution **
void Rand::genmultinomial(int N, const vector<double>& p, vector<int>& n) {

  int K = (int)p.size();
  n.resize(K);

  //gsl_ran_multinomial(RandomNumberGenerator, K, N, p, n);
  double norm = 0.0;
  double sum_p = 0.0;
  unsigned int sum_n = 0;

  /* p[k] may contain non-negative weights that do not sum to 1.0.
   * Even a probability distribution will not exactly sum to 1.0
   * due to rounding errors. 
   */
  
  for (int k = 0; k < K; k++)
    {
      norm += p[k];
    }
  
  for (int k = 0; k < K; k++)
    {
      if (p[k] > 0.0)
        {
          n[k] = gsl_ran_binomial (RandomNumberGenerator, p[k] / (norm - sum_p), N - sum_n);
        }
      else
        {
          n[k] = 0;
        }
      
      sum_p += p[k];
      sum_n += n[k];
    }
}

/// ** sample from discrete probability distribution given by probs **
int Rand::SampleFromDiscrete( const double probs[] , int numberofelements )
{
  double* cdf = new double[ numberofelements ];
  cdf[0] = probs[0];
  for( int i = 1; i < numberofelements; i++ )
    cdf[i] = (cdf[i-1] + probs[i]); 
  //   for( int i = 0; i < numberofelements; i++ )
  //     cdf[i] /= cdf[ numberofelements - 1 ];
  double u = myrand()* cdf[ numberofelements - 1 ];
  int k = 0;
  while( u > cdf[k] )
    ++k;
  delete[] cdf;
  return(k);
}

/// Sample from a discrete probability distribution given by probabilities,
/// specialised for the case in which the number of elements is small
int Rand::SampleFromDiscreteFast(const double *probs, int numberofelements) {

  // If the number of elements is small enough, the time spent recomputing
  // the cdf is smaller than the memory allocation overhead.

  int k = 0;
  double cdf = probs[0];
  while (++k < numberofelements)
    cdf += probs[k];

  double u = myrand() * cdf;
  k = 0;
  cdf = probs[0];
  while (u > cdf)
    cdf += probs[++k];

  return k;
}


/// ** Dirichlet distribution **
#if 0
    void Rand::gendirichlet(const size_t K, const double alpha[], double theta[] )
#else
    template<typename ConstVecType, typename VecType> \
	void Rand::gendirichlet( size_t K, const ConstVecType & alpha, VecType & theta )
#endif
  {
    #if 0 // We can't do this because the template is used with "double *" arrays.
	gp_assert_eq( theta.size(), K );
	gp_assert_eq( alpha.size(), K );
    #endif

    //bool invalid = false;
  //do{
    //invalid = false;
    double sum = 0.0;

    #if INEXPLICABLY_THROW_AWAY_ERROR_INFO
	try{
    #endif

      for( unsigned int i = 0; i < K; i++ ) {
	if( alpha[i] > 0 ){
	  theta[i] = gengam( alpha[i], 1.0 );
	}
	else theta[i] = 0.0;
	sum += theta[i]; 
      }

    #if INEXPLICABLY_THROW_AWAY_ERROR_INFO
	}
	catch(...){
	  throw string ("Error in gendirichlet");
	}
    #endif

    if( sum > 0.0 )
      for( unsigned int i = 0; i < K; i++ ){
	theta[i] /= sum;
	//if(theta[i] ==1 || theta[i] == 0){
	//invalid = true;break;
	  //problem here is that while theta[i] may be large enough to be considered nonzero, it may be small enough that sum 
	  //is not increased by enough to prevent one of the elements of theta being set to 1
	//}
      }

    else throw string("all gamma draws zero in function gendirichlet"); 
    //}
    //while(invalid);
  
}


// Instantiate templates:
template void Rand::gendirichlet<double *,double *>( size_t K, double * const & alpha, double * & theta );
template void Rand::gendirichlet<const double *,double *>( size_t K, const double * const & alpha, double * & theta );
template void Rand::gendirichlet<double *, bclib::pvector<double> >( size_t K, double * const & alpha, pvector<double> & theta );
template void Rand::gendirichlet< std::vector<double>, std::vector<double> >( size_t K, const std::vector<double> & alpha, std::vector<double> & theta );
template void Rand::gendirichlet< std::vector<double>, bclib::pvector<double> >( size_t K, const std::vector<double> & alpha, pvector<double> & theta );



END_BCLIB_NAMESPACE
