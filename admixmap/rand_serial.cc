/** 
 *   ADMIXMAP
 *   rand_serial.cc 
 *   Serial random number generators for admixmap
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include <cmath>
#include <iostream>
#include "rand.h"

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
}

using namespace std;

static gsl_rng *RandomNumberGenerator = gsl_rng_alloc( gsl_rng_taus );

double myrand()
{
  return( gsl_rng_uniform( RandomNumberGenerator ) );
}

double myrandRange( double Min, double Max )
{
  double Range = Max - Min;
  return( Min + Range * gsl_rng_uniform( RandomNumberGenerator ) );
}

void smyrand( long seed )
{
  gsl_rng_set(RandomNumberGenerator,
	      static_cast< unsigned long int >( seed ) );
}

// ** Gamma distribution **
double gengam( double shape, double rate )
{
  double x = 0.0;
  do
    x =  gsl_ran_gamma( RandomNumberGenerator, shape, 1.0 / rate ) ;
  while (x < 0.000001);
  return x;
}
// ** Beta distribution **
double genbet( double aa, double bb )
{
  return( gsl_ran_beta( RandomNumberGenerator, aa, bb ) );
}
// ** (univariate) Normal distribution **
double gennor( double av, double sd )
{
  return( av + gsl_ran_gaussian( RandomNumberGenerator, sd ) );
}
// ** Binomial distribution **
int genbinomial( int n, double p )
{
  return( gsl_ran_binomial( RandomNumberGenerator, p, n ) );
}
// ** Poisson distribution **
unsigned int genpoi( double mu )
{
  return( gsl_ran_poisson( RandomNumberGenerator, mu ) );
}

long ignpoi( double mu )
{
   return( gsl_ran_poisson( RandomNumberGenerator, mu ) );
}

// ** Multinomial distribution **
std::vector<int> genmultinomial2(int N, const std::vector<double> theta)
{
  int K = (int)theta.size();
  unsigned* n = new unsigned[ K ];
  double* p = new double[ K ];
  std::vector<int> sample( K );
  for(int i = 0; i < K; i++){
    p[i] = theta[i];
  }
  gsl_ran_multinomial(RandomNumberGenerator, K, N, p, n);
  for(int i = 0; i < K; i++){
    sample[i] = n[i];
  }
  delete[] n;
  delete[] p;
  return( sample );
}

// ** sample from discrete probability distribution gievn by probs **
int SampleFromDiscrete( const double probs[] , int numberofelements)
{
  double* cdf = new double[ numberofelements ];
  cdf[0] = probs[0];
  for( int i = 1; i < numberofelements; i++ )
    cdf[i] = (cdf[i-1] + probs[i]); 
  for( int i = 0; i < numberofelements; i++ )
    cdf[i] /= cdf[ numberofelements - 1 ];
  double u = myrand();
  int k = 0;
  while( u > cdf[k] )
    k++;
  delete[] cdf;
  return(k);
}
// ** Dirichlet distribution **
void gendirichlet(const size_t K, const double alpha[], double theta[] ) {
  double sum = 0.0;
  for( unsigned int i = 0; i < K; i++ ) {
    if( alpha[i] > 0 )
      theta[i] = gengam( alpha[i], 1.0 );
    else theta[i] = 0.0;
    sum += theta[i]; 
  }
  if( sum > 0.0 ) {
    for( unsigned int i = 0; i < K; i++ )
      theta[i] /= sum;
  } else throw string("all gamma draws zero in function gendirichlet"); 
}


