#include <cassert>
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

double gengam( double bb, double aa )
{
  double x = 0.0;
  do
   x =  gsl_ran_gamma( RandomNumberGenerator, aa, 1.0 / bb ) ;
   while (x < 0.000001);
  return x;
}

double genbet( double aa, double bb )
{
   return( gsl_ran_beta( RandomNumberGenerator, aa, bb ) );
}

double gennor( double av, double sd )
{
   return( av + gsl_ran_gaussian( RandomNumberGenerator, sd ) );
}

int genbinomial( int n, double p )
{
   return( gsl_ran_binomial( RandomNumberGenerator, p, n ) );
}

unsigned int genpoi( double mu )
{
   return( gsl_ran_poisson( RandomNumberGenerator, mu ) );
}

std::vector<int> genmultinomial2(int N, std::vector<double> theta)
{
  int K = (int)theta.size();
  unsigned int n[ K ];
  double p[ K ];
  std::vector<int> sample( K );
  for(int i = 0; i < K; i++){
      p[i] = theta[i];
  }
  gsl_ran_multinomial(RandomNumberGenerator, K, N, p, n);
  for(int i = 0; i < K; i++){
    sample[i] = n[i];
  }
  return( sample );
}

double MultinomialLikelihood( std::vector<int> r, std::vector<double> theta )
{
   if( r.size() != theta.size() ){
      cout << "Length of vector arguments to MultinomialLikelihood not equal." << endl;
      exit(0);
   }
   unsigned K = (int)r.size();
   unsigned int n[ K ];
   double p[ K ];
   for( unsigned i = 0; i < K; i++ ){
      p[i] = theta[i];
      n[i] = r[i];
   }
   return( gsl_ran_multinomial_pdf( K, p , n ) );
}

long ignpoi( double mu )
{
   return( gsl_ran_poisson( RandomNumberGenerator, mu ) );
}

int SampleFromDiscrete( double probs[] , int numberofelements)
{
   double cdf[ numberofelements ];
   cdf[0] = probs[0];
   for( int i = 1; i < numberofelements; i++ )
     cdf[i] = (cdf[i-1] + probs[i]); 
  for( int i = 0; i < numberofelements; i++ )
    cdf[i] /= cdf[ numberofelements - 1 ];
   double u = myrand();
   int k = 0;
   while( u > cdf[k] )
      k++;
    return(k);
}

void gendirichlet(const size_t K, const double alpha[], double theta[] )
{
   double sum = 0.0;

   for( unsigned int i = 0; i < K; i++ )
     {
       if( alpha[i] > 0 )
	 theta[i] = gengam( 1.0, alpha[i] );
	 
       else
         theta[i] = 0.0;
       sum += theta[i]; 
     }
   //need to do proper error catching here
   assert( sum > 0.0 );
   for( unsigned int i = 0; i < K; i++ )
     theta[i] /= sum;
   
}


void ddigam(  double *X, double *ddgam  )
{
//  FOLLOWS CLOSELY ALG AS 103 APPL.STATS.(1976) VOL 25
//  CALCS DIGAMMA(X)=D(LOG(GAMMA(X)))/DX
//
//  SET CONSTANTS.SN=NTH STIRLING COEFFICIENT,D1=DIGAMMA(1.)
//
   double S, C, S3, S4, S5, D1, Y, R;

   S = 1.0e-5;
   C = 8.5e0;
   S3 = 8.333333333e-2; 
   S4 = 8.333333333e-3;
   S5 = 3.968253968e-3;
   D1 = -0.5772156649;

//      DATA S,C,S3,S4,S5,D1/1.0D-5,8.5D0,8.333333333D-2,
//    1  8.333333333D-3,3.968253968D-3,-0.5772156649D0/


//  CHECK ARGUMENT IS POSITIVE

   *ddgam=0.0;
   Y=*X;
   if(Y < 0.0){
      std::cout << " ARG OF DIGAMMA FN = " << Y << ". ARG MUST BE POSITIVE " << std::endl;
      exit(1);}

//  USE APPROXIMATION IF ARGUMENT .LE.S

   if(Y > S){

//  REDUCE TO DIGAMMA(X+N),(X+N).GE.C

      while( Y < C ){
         *ddgam=*ddgam-(1.0/Y);
         Y=Y+1.0;}
   
//  USE STIRLING IF ARGUMENT .GE.C

      R=1.0/Y;
      *ddgam=*ddgam+log(Y)-0.5*R;
      R=R*R;
      *ddgam=*ddgam-R*(S3-R*(S4-R*S5));}
   else
      *ddgam=D1-1.0/Y;
}

void trigam( double *x, double *trgam )
{
/*
 * closely follows alg. as 121 appl.stats. (1978) 
 * vol 27, 97-99. (b.e. schneider)
 *
 * calculates trigamma(x)=d**2(log(gamma(x)))/dx**2
 */
   double a=1.0e-4,b=5.0,one=1.0,half=0.5,y,z,trigam1=0.0;
   double b2=0.1666666667,b4=-0.03333333333;
   double b6=0.02380952381,b8=-0.0333333333;
/*
 *  b2,b4,b6,b8 are bernoulli numbers
 *
 *  check that argument is positive
 *
 */ 
   z=*x;
/*
 *  use small value approximation if x.le.a
 */
   if(z<=a){ 
      trigam1=1.0/(z*z);
      *trgam=trigam1;}
   else{
/*
 *  increase argument to (x+i).ge.b
 */
      while(z<b){
         trigam1=trigam1+1.0/(z*z);
         z=z+1.0;
      }
/*
 *  apply asymptotic formula if argument.ge.b
 */
      y=1.0/(z*z);
      trigam1=trigam1+half*y+(one+y*(b2+y*(b4+y*(b6+y*b8))))/z;
      *trgam=trigam1;
   }
}
