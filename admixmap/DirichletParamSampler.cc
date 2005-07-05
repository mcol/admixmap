#include "DirichletParamSampler.h"

using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

DirichletParamSampler::DirichletParamSampler()
{
   TuneEta.SetParameters( 10, .1, 0.1, 100, 0.44 );
   TuneMu.SetParameters( 10, 0.002, 0.0001, 0.1, 0.23 );
   EtaAlpha = .1;
   EtaBeta = .1;
}

DirichletParamSampler::DirichletParamSampler( unsigned int ind )
{
   d = ind;
   gamma = new double[d];
   munew = new double[d];
   TuneEta.SetParameters( 10, 0.1, 0.01, 100, 0.44 );
   TuneMu.SetParameters( 10, 0.001, 0.0001, 0.1, 0.44 );
   EtaAlpha = .1;
   EtaBeta = .1;
   for( unsigned int i = 0; i < d; i++ )
      gamma[i] = 1.0;
}

void DirichletParamSampler::SetSize( unsigned int ind )
{
   d = ind;
   gamma = new double[d];
   munew = new double[d];
   for( unsigned int i = 0; i < d; i++ )
      gamma[i] = 1.0;
}

DirichletParamSampler::~DirichletParamSampler()
{
   delete [] munew;
   delete [] gamma;
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
   double *alpha, L1=0, L2=0, P1=0, P2=0, Proposal1=0, Proposal2=0;
   alpha = new double[d];
   double sigma = TuneMu.GetSigma();
   for( i = 0; i < d; i++ )
      alpha[i] = mu[i] / sigma;
   gendirichlet( d, alpha, munew );

   for( i = 0; i < d; i++ ){
      Proposal1 += (munew[i]/sigma - 1) * log( mu[i] ) - gsl_sf_lngamma( munew[i]/sigma );
      Proposal2 += (alpha[i] - 1) * log( munew[i] ) - gsl_sf_lngamma( alpha[i] );
      L1 += (*eta * mu[i] - 1) * sumlogtheta[i] - n*gsl_sf_lngamma( *eta * mu[i] );
      L2 += (*eta * munew[i] - 1) * sumlogtheta[i] - n*gsl_sf_lngamma( *eta * munew[i] );
      P1 += (gamma[i] - 1) * log( mu[i] );
      P2 += (gamma[i] - 1) * log( munew[i] );
//      PR(i);PR(mu[i]);PR(munew[i]);PR(Proposal1);PR(Proposal2);PR(L1);PR(L2);PR(P1);PR(P2);
   }
   if( log(myrand()) <  P2 + L2 - P1 - L1 - Proposal2 + Proposal1 ){
      TuneMu.Event(true);
      for( i = 0; i < d; i++ )
         mu[i] = munew[i];
   }
   else
      TuneMu.Event(false);
   
   etanew = exp( gennor( log( *eta ), TuneEta.GetSigma() ) );
   Proposal1 = log(etanew) - log(*eta);
   P1 = ( EtaAlpha - 1.0 ) * ( log(etanew) - log(*eta) ) - EtaBeta * ( etanew - *eta );
   L1 = n * ( gsl_sf_lngamma( etanew ) - gsl_sf_lngamma( *eta ) );
   for( i = 0; i < d; i++ )
      L1 += mu[i] * (etanew - *eta) * sumlogtheta[i] - n*gsl_sf_lngamma( etanew * mu[i] ) + n*gsl_sf_lngamma( *eta * mu[i] );

//   PR(etanew);PR(*eta);PR(L1);PR(P1);PR(Proposal1);

   if( log(myrand()) < L1 + P1 + Proposal1 ){
      *eta = etanew;
      TuneEta.Event(true);
   }
   else
      TuneEta.Event(false);
//   PR(*eta);
//   cout << endl;
}
