#include "functions.h"
#include <cassert>
#include <gsl/gsl_sf_gamma.h>

double getGammaLogDensity(double alpha, double beta, double x)
{
  return alpha * log(beta) - gsl_sf_lngamma(alpha) + (alpha-1.0)*log(x) - beta*x;
}

double getDirichletLogDensity(const Vector_d& a, const Vector_d& x)
{
  double f;
  Vector_d xdash;    //This adds the implied last (kth) element of x
  xdash = x;                              //when x has length (k-1)
  int k = a.GetNumberOfElements();        //x must sum to 1
  if (k - 1 == x.GetNumberOfElements()) {
    xdash.AddElement( k - 1 );
    xdash( k - 1 ) = 1 - x.Sum();
  }
  assert( k == xdash.GetNumberOfElements() ); //Error in getDirichletLogDensity - lengths of vector arguments do not match.\n";

  f = gsl_sf_lngamma( a.Sum() );
  for( int i = 0; i < k; i++ )
    if( a(i) > 0.0 )
      f += ( a(i) - 1 ) * log( xdash(i) ) - gsl_sf_lngamma( a(i) );

  return f;
}

double AverageOfLogs(const std::vector<double>& vec, double max)
{
  double sum = 0;

  for ( unsigned int i = 0; i < vec.size(); i++ )
    sum += exp( vec[i] - max );

  sum /= vec.size();

  return log(sum) + max;
}
