#include "DirichletParamSampler.h"

using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

DirichletParamSampler::DirichletParamSampler()
{
   TuneEta.SetParameters( 10, .1, 0.1, 100, 0.44 );
   EtaAlpha = 1;
   EtaBeta = 1;
}

DirichletParamSampler::DirichletParamSampler( unsigned int ind )
{
   Matrix_i empty_i(1,1);
   Matrix_d empty_d(1,1);
   d = ind;
   gamma = new double[d];
   munew = new double[d];
   TuneEta.SetParameters( 10, 0.1, 0.01, 100, 0.44 );
   EtaAlpha = 1;
   EtaBeta = 1;
   for( unsigned int i = 0; i < d; i++ )
      gamma[i] = 1.0;
   
  DirParamArray = new DARS*[ d ];
  for( unsigned int j = 0; j < d; j++ ){
     DirParamArray[j] = new DARS();
     DirParamArray[j]->SetParameters( 0, 0, 0.1, AlphaParameters,5,
                                      logf, dlogf, ddlogf, empty_i, empty_d );
  }
}

void DirichletParamSampler::SetSize( unsigned int ind )
{
   Matrix_i empty_i(1,1);
   Matrix_d empty_d(1,1);
   d = ind;
   gamma = new double[d];
   munew = new double[d];
   for( unsigned int i = 0; i < d; i++ )
      gamma[i] = 1.0;

   DirParamArray = new DARS*[ d ];
   for( unsigned int j = 0; j < d; j++ ){
      DirParamArray[j] = new DARS();
      DirParamArray[j]->SetParameters( 0, 0, 0.1, AlphaParameters,5,
                                       logf, dlogf, ddlogf, empty_i, empty_d );
   }
}

DirichletParamSampler::~DirichletParamSampler()
{
   delete [] munew;
   delete [] gamma;
   for(unsigned int i=0; i<d; i++){
      delete DirParamArray[i];
   }
}

void DirichletParamSampler::SetPriorEta( double inEtaAlpha, double inEtaBeta )
{
   EtaAlpha = inEtaAlpha;
   EtaBeta = inEtaBeta;
}

void DirichletParamSampler::SetPriorMu( double *ingamma )
{
   for( unsigned int i = 0; i < d; i++ ){
      gamma[i] = ingamma[i];
   }
}

void DirichletParamSampler::Sample( unsigned int n, double *sumlogtheta, double *eta, double *mu )
/*
n = number of gametes/individuals
*/
{
   unsigned int i;
   double L1=0, P1=0, Proposal1=0;
   double b = mu[d-1] + mu[0] - 0.005;
   double summu = 1.0 - mu[d-1];
   AlphaParameters[0] = n;
   for( unsigned int j = 0; j < d-1; j++ ){
      AlphaParameters[1] = *eta;
      AlphaParameters[2] = summu - mu[j];
      AlphaParameters[3] = sumlogtheta[d-1];
      AlphaParameters[4] = sumlogtheta[j];
      DirParamArray[j]->SetLeftTruncation( 0.005 );
      DirParamArray[j]->SetRightTruncation( b );
      // elements of Dirichlet parameter vector are updated one at a time
      DirParamArray[j]->UpdateParameters( AlphaParameters, 5 );
      mu[j] = DirParamArray[j]->Sample();
      b = b - mu[j] + mu[j+1];
      summu = AlphaParameters[2] + mu[j];
   }
   mu[d-1] = 1.0 - summu;
   
   etanew = exp( gennor( log( *eta ), TuneEta.GetSigma() ) );
   Proposal1 = log(etanew) - log(*eta);
   P1 = ( EtaAlpha - 1.0 ) * ( log(etanew) - log(*eta) ) - EtaBeta * ( etanew - *eta );
   L1 = n * ( gsl_sf_lngamma( etanew ) - gsl_sf_lngamma( *eta ) );
   for( i = 0; i < d; i++ )
      L1 += mu[i] * (etanew - *eta) * sumlogtheta[i] - n*gsl_sf_lngamma( etanew * mu[i] ) + n*gsl_sf_lngamma( *eta * mu[i] );

   if( log(myrand()) < L1 + P1 + Proposal1 ){
      *eta = etanew;
      TuneEta.Event(true);
   }
   else
      TuneEta.Event(false);
}

// these 3 functions calculate log-likelihood and derivatives for adaptive rejection sampling of 
// Dirichlet population admixture parameters
double
DirichletParamSampler::logf( Vector_d &parameters , Matrix_i&, Matrix_d&, double x )
{
   int n = (int)parameters(0);
   double eta = parameters(1), summu = parameters(2), sumlj = parameters(4), sumln = parameters(3);
   double f = eta * ( x*sumlj + (1.0-summu-x)*sumln )
      - n * ( gsl_sf_lngamma(x*eta) + gsl_sf_lngamma((1.0-summu-x)*eta) );
   
  return f;
}

double
DirichletParamSampler::dlogf( Vector_d &parameters, Matrix_i&, Matrix_d&, double x )
{
  double f,x2,y1,y2;
  int n = (int)parameters(0);
  double eta = parameters(1), summu = parameters(2), sumlj = parameters(4), sumln = parameters(3);
  
  x2 = eta*x;
  if(x2 < 0)cout<<"\nError in  DirichletParamSampler::dlogf - arg x to ddigam is negative\n";   
  ddigam( &x2 , &y1 );

  x2 = eta*(1.0-x-summu);
  if(x2 < 0)cout<<"\nError in  DirichletParamSampler::dlogf - arg x2 to ddigam is negative\n";   
  ddigam( &x2 , &y2 );
  
  f =  eta * ( sumlj - sumln ) - n * eta * ( y1 - y2 );
  
  return f;
}

double
DirichletParamSampler::ddlogf( Vector_d &parameters, Matrix_i&, Matrix_d&, double x )
{
  double f,x2,y1,y2;
  int n = (int)parameters(0);
  double eta = parameters(1), summu = parameters(2);
  
  x2 = eta*x;
  trigam( &x2, &y1 );
  x2 = eta*(1.0-x-summu);
  trigam( &x2, &y2 );
  
  f = -n*eta*eta*( y2+y1 );
  
  return(f);
}
