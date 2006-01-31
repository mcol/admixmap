#include "DirichletParamSampler.h"
#include <algorithm>
#include <numeric>

using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

DirichletParamSampler::DirichletParamSampler()
{
#if SAMPLERTYPE==1
  step0 = 0.1; //sd of proposal distribution for log eta
  // need to choose sensible value for this initial RW sd
  step = step0;
  TuneEta.SetParameters( step0, 0.01, 10, 0.44); 
  mu = 0;
  munew = 0;
  muDirichletParams = 0;
  DirParamArray = 0;
#elif SAMPLERTYPE==2
  logalpha = 0;
  initialAlphaStepsize = 0.05;//need a way of setting this without recompiling, or a sensible fixed value
  targetAlphaAcceptRate = 0.44;//need to choose suitable value for this
#endif
}

DirichletParamSampler::DirichletParamSampler( unsigned numind, unsigned numpops)
{
  SetSize(numind, numpops);
}

void DirichletParamSampler::SetSize( unsigned numobs, unsigned numpops)
  // sets number of elements in Dirichlet parameter vector
  // instantiates an adaptive rejection sampler object for each element
//sets up MuSampler object 
{
   K = numpops;
#if SAMPLERTYPE==1
   EtaAlpha = K; // for compatibility with gamma(1, 1) prior on alpha
   EtaBeta = 1;
   //Dirichlet(1, ..., 1) prior on proportions mu
   muDirichletParams = new double[K];
   mu = new double[K];
   munew = new double[K];
   for( unsigned int i = 0; i < K; i++ )
     muDirichletParams[i] = 1.0;
   DirParamArray = new AdaptiveRejection*[ K ];
   for( unsigned int j = 0; j < K; j++ ){
     DirParamArray[j] = new AdaptiveRejection();
     DirParamArray[j]->Initialise(true, true, 0.0, 1.0, logf, dlogf);
     DirParamArray[j]->setLowerBound(0.00);
   }
   
#elif SAMPLERTYPE==2
   logalpha = new double[K];
   //transform(alpha[0].begin(), alpha[0].end(), logalpha, xlog);//logalpha = log(alpha)
   
   //elem 0 is sum of log admixture props
   AlphaArgs.n = numobs; //num individuals/gametes will be passed as arg to sample method
   AlphaArgs.dim = K;
   //if( options->isRandomMatingModel() )AlphaArgs.n *= 2;
   AlphaArgs.eps0 = 1.0; //Gamma(1, 1) prior on alpha
   AlphaArgs.eps1 = 1.0;
   AlphaSampler.SetDimensions(K, initialAlphaStepsize, 0.01, 10.0, 20, targetAlphaAcceptRate, findE, gradE);
#endif
}

DirichletParamSampler::~DirichletParamSampler()
{
#if SAMPLERTYPE==1
  delete[] mu;
  delete [] munew;
  delete [] muDirichletParams;
  if(DirParamArray){ //in case not allocated
    for(unsigned int i = 0; i < K; i++){
      if(DirParamArray[i])
	delete DirParamArray[i];
    }
    delete[] DirParamArray;
  }

#elif SAMPLERTYPE==2
  delete[] logalpha;
#endif
}

#if SAMPLERTYPE==1
void DirichletParamSampler::SetPriorEta( double inEtaAlpha, double inEtaBeta ) {
  EtaAlpha = inEtaAlpha;
  EtaBeta = inEtaBeta;
}
void DirichletParamSampler::SetPriorMu( const double* const ingamma ) {
   for( unsigned int i = 0; i < K; i++ ){
      muDirichletParams[i] = ingamma[i];
   }
}
#endif

void DirichletParamSampler::Sample( unsigned int n, const double* const sumlogtheta, std::vector<double> *alpha )
/*
  n = number of observations
  sumlogtheta = sum log theta from Individuals
*/
{
  // *** sample elements of mu with adaptive rejection sampler, conditional on frequencies sumlogtheta, and eta with RW
#if SAMPLERTYPE==1
  eta = accumulate(alpha->begin(), alpha->end(), 0.0, std::plus<double>());//eta = sum of alpha[0]
  for( unsigned i = 0; i < K; i++ ){
    mu[i] = (*alpha)[i]/eta;
  }
  
  double b = mu[K-1] + mu[0];//upper bound for sampler
  double summu = 1.0 - mu[K-1];
  AlphaParameters[0] = n;
  for( unsigned int j = 0; j < K-1; j++ ){
    AlphaParameters[1] = eta; // dispersion parameter
    AlphaParameters[2] = summu - mu[j]; // 1 - last proportion parameter
    AlphaParameters[3] = sumlogtheta[K-1]; 
    AlphaParameters[4] = sumlogtheta[j];
    
    DirParamArray[j]->setUpperBound(b);
    // Dirichlet proportion parameters mu[j] are updated one at a time
    mu[j] = DirParamArray[j]->Sample(AlphaParameters, ddlogf);
    b = b - mu[j] + mu[j+1];
    summu = AlphaParameters[2] + mu[j];
  }
  mu[K-1] = 1.0 - summu;
  
  SampleEta(n, sumlogtheta, &eta, mu);
  
  for( unsigned j = 0; j < K; j++ )
    (*alpha)[j] = mu[j]*eta;
  
  
#elif SAMPLERTYPE==2
  // *** Hamiltonian sampler for alpha
  AlphaArgs.n = n;
  for (unsigned int k = 0; k < K; ++j) {
    AlphaArgs.sumlogtheta[k] = sumlogtheta[k];
  }
  AlphaSampler.Sample(logalpha, &AlphaArgs);//sample new values for logalpha
  transform(logalpha, logalpha+K, alpha[0].begin(), xexp);//alpha = exp(logalpha)
#endif
  
}

#if SAMPLERTYPE==1
void DirichletParamSampler::SampleEta(unsigned n, const double* const sumlogtheta, double *eta, const double* const mu){
  // Dirichlet dispersion parameter eta is updated with a Metropolis random walk
  unsigned int i;
  double L1=0, P1=0, Proposal1=0;

  etanew = exp( gennor( log( *eta ), step ) );
  Proposal1 = log(etanew) - log(*eta);
  // log prior ratio P1 
  P1 = ( EtaAlpha - 1.0 ) * ( log(etanew) - log(*eta) ) - EtaBeta * ( etanew - *eta );
  // log likelihood ratio L1
  L1 = n * ( gsl_sf_lngamma( etanew ) - gsl_sf_lngamma( *eta ) );
  for( i = 0; i < K; i++ )
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

double DirichletParamSampler::getEtaStepSize()const
{
    return TuneEta.getStepSize();
}

double DirichletParamSampler::getEtaExpectedAcceptanceRate()const
{
    return TuneEta.getExpectedAcceptanceRate();
}

// these 3 functions calculate log-likelihood and derivatives for adaptive rejection sampling of 
// Dirichlet proportion parameters
double DirichletParamSampler::logf( double x, const void* const pars )
{
  const double* parameters = (const double*) pars;
   int n = (int)parameters[0];
   double eta = parameters[1], summu = parameters[2], sumlj = parameters[4], sumln = parameters[3];
   double f = eta * ( x*sumlj + (1.0-summu-x)*sumln )
      - n * ( gsl_sf_lngamma(x*eta) + gsl_sf_lngamma((1.0-summu-x)*eta) );
  return f;
}

double DirichletParamSampler::dlogf( double x, const void* const pars )
{
  const double* parameters = (const double*) pars;
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

double DirichletParamSampler::ddlogf( double x, const void* const pars)
{
  const double* parameters = (const double*) pars;
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

#elif SAMPLERTYPE==2
double DirichletParamSampler::getAlphaSamplerStepsize()const {
    return AlphaSampler.getStepSize();
}

double DirichletParamSampler::getAlphaSamplerAcceptanceRate()const {
    return AlphaSampler.getExpectedAcceptanceRate();
}

//calculate objective function (-log posterior) for log alpha, used in Hamiltonian Metropolis algorithm
double DirichletParamSampler::findE(const double* const theta, const void* const vargs){
  /*
    theta = log dirichlet parameters (alpha)
    n = #individuals/gametes
    sumlogtheta (array, length dim) = sums of logs of individual admixture proportions
    eps0, eps1 = parameters of Gamma prior for alpha
  */
  const AlphaSamplerArgs* args = (const AlphaSamplerArgs*)vargs;

  double E = 0.0;
  double sumalpha = 0.0, sumgamma = 0.0, sumtheta = 0.0, sume = 0.0;
  bool flag = true;
  for(int j = 0; j < args->dim; ++j){
    if(exp(theta[j]) == 0.0){flag = false;break;} //to avoid underflow problems
    sumalpha += exp(theta[j]);
    sumgamma += gsl_sf_lngamma(exp(theta[j]));
    sume += exp(theta[j]) * (args->eps1 - args->sumlogtheta[j]);
    sumtheta += theta[j];
  }
  if(flag){
    E = args->n * (gsl_sf_lngamma(sumalpha) - sumgamma) - sume + args->eps0 * sumtheta;
    return -E;
  }
  else return -1.0;//is there a better return value? possibly use flag pointer
}

//calculate gradient for log alpha
void DirichletParamSampler::gradE(const double* const theta, const void* const vargs, double *g){
  const AlphaSamplerArgs* args = (const AlphaSamplerArgs*)vargs;
  double sumalpha = 0.0, x, y1, y2;
  for(int j = 0; j < args->dim; ++j) {
    g[j] = 0.0;
    sumalpha += exp(theta[j]);
  }
    ddigam(&sumalpha, &y1);
    for(int j = 0; j < args->dim; ++j) {
      x = exp(theta[j]);
      ddigam(&x, &y2);
      if(x > 0.0 && gsl_finite(y1) && gsl_finite(y2)){//to avoid over/underflow problems
	g[j] = x *( args->n *(y2 - y1) + (args->eps1 - args->sumlogtheta[j])) - args->eps0;
      }
    }
}


#endif
