#include "DirichletParamSampler.h"

using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

DirichletParamSampler::DirichletParamSampler()
{
  step0 = 0.1; //sd of proposal distribution for log eta
  // need to choose sensible value for this initial RW sd
  step = step0;
  TuneEta.SetParameters( step0, 0.01, 10, 0.44); 
  EtaAlpha = 1;
  EtaBeta = 1;
}

DirichletParamSampler::DirichletParamSampler( unsigned numind, unsigned numpops )
{
  SetSize(numind, numpops);
}

void DirichletParamSampler::SetSize( unsigned numind, unsigned numpops )
  // sets number of elements in Dirichlet parameter vector
  // instantiates an adaptive rejection sampler object for each element
//sets up MuSampler object 
{
   d = numpops;
   gamma = new double[d];
   munew = new double[d];
   for( unsigned int i = 0; i < d; i++ )
      gamma[i] = 1.0;

   DirParamArray = new NewDARS*[ d ];
   for( unsigned int j = 0; j < d; j++ ){
      DirParamArray[j] = new NewDARS();
      //DirParamArray[j]->SetParameters( 0, 0, 0.1, AlphaParameters, logf, dlogf, ddlogf, 0, 0 );
      DirParamArray[j]->Initialise(true, true, 0.0, 1.0, logf, dlogf);
   }
   //muSampler.setDimensions(numind, numpops, 0.00001 , 0.0, 1.0, 0.44);//may need to modify initial stepsize (arg 3)
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

//sample mus with DARS, conditional on frequencies sumlogtheta, and eta with RW
void DirichletParamSampler::Sample( unsigned int n, double *sumlogtheta, double *eta, double *mu )
/*
  n = number of observations
*/
{

  double b = mu[d-1] + mu[0];//upper bound for sampler
   double summu = 1.0 - mu[d-1];
   AlphaParameters[0] = n;
   for( unsigned int j = 0; j < d-1; j++ ){
     AlphaParameters[1] = *eta; // dispersion parameter
     AlphaParameters[2] = summu - mu[j]; // 1 - last proportion parameter
     AlphaParameters[3] = sumlogtheta[d-1]; 
     AlphaParameters[4] = sumlogtheta[j];
     //DirParamArray[j]->SetLeftTruncation( 0.005 );
     DirParamArray[j]->setLowerBound(0.00);
     //DirParamArray[j]->SetRightTruncation( b );
     DirParamArray[j]->setUpperBound(b);
     // Dirichlet proportion parameters mu[j] are updated one at a time
     //DirParamArray[j]->UpdateParameters( AlphaParameters );
     mu[j] = DirParamArray[j]->Sample(AlphaParameters, ddlogf);
     b = b - mu[j] + mu[j+1];
     summu = AlphaParameters[2] + mu[j];
   }
   mu[d-1] = 1.0 - summu;
   
   SampleEta(n, sumlogtheta, eta, mu);

}

// //sample mu with MuSampler (DARS for dim==2, HMC otherwise), conditional on counts, and eta with RW
// void DirichletParamSampler::Sample2(unsigned n, double *sumlogtheta, double *eta, double *mu, int *counts){
//   for(unsigned k=0; k < d; ++k)mu[k] *= *eta;
//   muSampler.Sample(mu, *eta, counts);
//   for(unsigned k=0; k < d; ++k)mu[k] /= *eta;
//   SampleEta(n, sumlogtheta, eta, mu);
// }

void DirichletParamSampler::SampleEta(unsigned n, double *sumlogtheta, double *eta, double *mu){
  // Dirichlet dispersion parameter eta is updated with a Metropolis random walk
  unsigned int i;
  double L1=0, P1=0, Proposal1=0;

  etanew = exp( gennor( log( *eta ), step ) );
  Proposal1 = log(etanew) - log(*eta);
  // log prior ratio P1 
  P1 = ( EtaAlpha - 1.0 ) * ( log(etanew) - log(*eta) ) - EtaBeta * ( etanew - *eta );
  // log likelihood ratio L1
  L1 = n * ( gsl_sf_lngamma( etanew ) - gsl_sf_lngamma( *eta ) );
  for( i = 0; i < d; i++ )
    L1 += mu[i] * (etanew - *eta) * sumlogtheta[i] - n*gsl_sf_lngamma( etanew * mu[i] ) + n*gsl_sf_lngamma( *eta * mu[i] );
  // calculate log acceptance probability ratio
  LogAccProb = 0.0;
  if(P1 + L1 + Proposal1 < 0.0)  
    LogAccProb = P1 + L1 + Proposal1; 
  //accept/reject proposal
  if( log(myrand()) < LogAccProb ){
    *eta = etanew;
  }
  //update step size
  step = TuneEta.UpdateStepSize( exp(LogAccProb) );
}

double DirichletParamSampler::getEtaStepSize()
{
    return TuneEta.getStepSize();
}

double DirichletParamSampler::getEtaExpectedAcceptanceRate()
{
    return TuneEta.getExpectedAcceptanceRate();
}
// double DirichletParamSampler::getMuStepSize()
// {
//     return muSampler.getStepsize();
// }

// double DirichletParamSampler::getMuExpectedAcceptanceRate()
// {
//     return muSampler.getAcceptanceRate();
// }


// these 3 functions calculate log-likelihood and derivatives for adaptive rejection sampling of 
// Dirichlet proportion parameters
double
DirichletParamSampler::logf( double x, const double* const parameters )
{
   int n = (int)parameters[0];
   double eta = parameters[1], summu = parameters[2], sumlj = parameters[4], sumln = parameters[3];
   double f = eta * ( x*sumlj + (1.0-summu-x)*sumln )
      - n * ( gsl_sf_lngamma(x*eta) + gsl_sf_lngamma((1.0-summu-x)*eta) );
   
  return f;
}

double
DirichletParamSampler::dlogf( double x, const double* const parameters )
{
  double f,x2,y1,y2;
  int n = (int)parameters[0];
  double eta = parameters[1], summu = parameters[2], sumlj = parameters[4], sumln = parameters[3];
  
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
DirichletParamSampler::ddlogf( double x, const double* const parameters)
{
  double f,x2,y1,y2;
  int n = (int)parameters[0];
  double eta = parameters[1], summu = parameters[2];
  
  x2 = eta*x;
  trigam( &x2, &y1 );
  x2 = eta*(1.0-x-summu);
  trigam( &x2, &y2 );
  
  f = -n*eta*eta*( y2+y1 );
  
  return(f);
}
