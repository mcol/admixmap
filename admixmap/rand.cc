/** 
 *   ADMIXMAP
 *   rand.cc 
 *   Random number generators for admixmap
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
#include <limits>

extern "C" {
  //#include <gsl/gsl_rng.h>
#ifdef PARALLEL
#include "gsl-sprng.h"//includes sprng.h and defines SIMPLE_SPRNG and USE_MPI
#endif
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
}

using namespace std;

#ifdef PARALLEL
  //allocate a sprng generator
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

//standard uniform number generator
double Rand::myrand()
{
  return( gsl_rng_uniform( RandomNumberGenerator ) );
  //?? should use gsl_rng_uniform_pos to exclude 0
}

// ** U(Min, Max) 
double Rand::myrandRange( double Min, double Max )
{
  double Range = Max - Min;
  return( Min + Range * gsl_rng_uniform( RandomNumberGenerator ) );
}

//set random number seed
void Rand::setSeed( long seed )
{
  gsl_rng_set(RandomNumberGenerator,
	      static_cast< unsigned long int >( seed ) );
  //in sprng20 case, calls init_sprng
}

// ** Gamma distribution **
double Rand::gengam( double shape, double rate )
{
  double x = 0.0;
  do
    x = gsl_ran_gamma( RandomNumberGenerator, shape, 1.0 / rate ) ;
  while (x < numeric_limits<double>::min( ));
  return x;
}
// ** Beta distribution **
double Rand::genbet( double aa, double bb )
{
  return( gsl_ran_beta( RandomNumberGenerator, aa, bb ) );
}
// ** (univariate) Normal distribution **
double Rand::gennor( double av, double sd )
{
  return( av + gsl_ran_gaussian( RandomNumberGenerator, sd ) );
}
// ** Binomial distribution **
int Rand::genbinomial( int n, double p )
{
  return( gsl_ran_binomial( RandomNumberGenerator, p, n ) );
}
// ** Poisson distribution **
unsigned int Rand::genpoi( double mu )
{
  return( gsl_ran_poisson( RandomNumberGenerator, mu ) );
}

long Rand::ignpoi( double mu )
{
   return( gsl_ran_poisson( RandomNumberGenerator, mu ) );
}

// ** Multinomial distribution **
std::vector<int> Rand::genmultinomial(int N, const std::vector<double> theta)
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
int Rand::SampleFromDiscrete( const double probs[] , int numberofelements)
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
void Rand::gendirichlet(const size_t K, const double alpha[], double theta[] ) {
  bool invalid = false;
  //do{
    invalid = false;
    double sum = 0.0;
    for( unsigned int i = 0; i < K; i++ ) {
      if( alpha[i] > 0 ){
	theta[i] = gengam( alpha[i], 1.0 );
      }
      else theta[i] = 0.0;
      sum += theta[i]; 
    }
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


