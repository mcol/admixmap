#include "DirichletParamSampler.h"
#include <algorithm>
#include <numeric>
#include "functions.h"
#include <gsl/gsl_math.h>

using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

DirichletParamSampler::DirichletParamSampler() {
  Initialise();
}

DirichletParamSampler::DirichletParamSampler( unsigned numind, unsigned numpops) {
  Initialise();
  SetSize(numind, numpops);
}

void DirichletParamSampler::Initialise() {
#if SAMPLERTYPE==1
//   step0 = 0.1; //sd of proposal distribution for log eta
//   // need to choose sensible value for this initial RW sd
//   step = step0;
//   TuneEta.SetParameters( step0, 0.01, 10, 0.44); 
  mu = 0;
  munew = 0;
  muDirichletParams = 0;
  //  DirParamArray = 0;
#elif SAMPLERTYPE==2
  logalpha = 0;
  initialAlphaStepsize = 0.02;//need a way of setting this without recompiling, or a sensible fixed value
  targetAlphaAcceptRate = 0.44;//need to choose suitable value for this
#endif
}

void DirichletParamSampler::SetSize( unsigned numobs, unsigned numpops)
// sets number of elements in Dirichlet parameter vector
// instantiates an adaptive rejection sampler object for each element
//sets up Sampler object
//  numobs = number of observations 
{
   K = numpops;
#if SAMPLERTYPE==1
   AlphaParameters[0] = numobs;
   muDirichletParams = new double[K];
   mu = new double[K];
   munew = new double[K];
   for( unsigned int i = 0; i < K; i++ )
     muDirichletParams[i] = 1.0;
   //DirParamArray = new AdaptiveRejection();
   DirParamArray.Initialise(true, true, 1.0, 0.00001, logf, dlogf);
//      DirParamArray[j]->setLowerBound(0.0);
//   for( unsigned int j = 0; j < K; j++ ){
//      DirParamArray[j] = new AdaptiveRejection();
//      DirParamArray[j]->Initialise(true, true, 0.0, 1.0, logf, dlogf);
//      DirParamArray[j]->setLowerBound(0.0);
//   }

   EtaArgs.priorshape = K; // for compatibility with gamma(1, 1) prior on alpha
   EtaArgs.priorrate = 1;
   EtaArgs.numpops = numpops;
   EtaArgs.numobs = numobs; 
   EtaSampler.SetDimensions(1, 0.01/*initial stepsize*/, 0.01/*min stepsize*/, 100.0/*max stepsize*/, 20/*num leapforgs*/, 
			    0.44/*target accept rate*/, etaEnergy, etaGradient);  
#elif SAMPLERTYPE==2
   logalpha = new double[K];
   
   //elem 0 is sum of log admixture props
   AlphaArgs.n = numobs; //num individuals/gametes will be passed as arg to sampler
   AlphaArgs.dim = K;
   //if( options->isRandomMatingModel() )AlphaArgs.n *= 2;
   AlphaArgs.eps0 = 1.0; //Gamma(1, 1) prior on alpha
   AlphaArgs.eps1 = 1.0;
   AlphaSampler.SetDimensions(K, initialAlphaStepsize, 0.01, 100.0, 20, targetAlphaAcceptRate, findE, gradE);
#endif
}

DirichletParamSampler::~DirichletParamSampler()
{
#if SAMPLERTYPE==1
  delete[] mu;
  delete [] munew;
  delete [] muDirichletParams;
  //if(DirParamArray) delete[] DirParamArray;

#elif SAMPLERTYPE==2
  delete[] logalpha;
#endif
}

#if SAMPLERTYPE==1
void DirichletParamSampler::SetPriorEta( double inEtaAlpha, double inEtaBeta ) {
//   EtaAlpha = inEtaAlpha;
//   EtaBeta = inEtaBeta;
   EtaArgs.priorshape = inEtaAlpha; 
   EtaArgs.priorrate = inEtaBeta;
}
void DirichletParamSampler::SetPriorMu( const double* const ingamma ) {
   for( unsigned int i = 0; i < K; i++ ){
      muDirichletParams[i] = ingamma[i];
   }
}
#endif

void DirichletParamSampler::Sample( const double* const sumlogtheta, std::vector<double> *alpha )
/*
  sumlogtheta = sum log theta from Individuals
*/
{
  // *** sample elements of mu with adaptive rejection sampler, conditional on frequencies sumlogtheta, and eta with RW
#if SAMPLERTYPE==1
  eta = accumulate(alpha->begin(), alpha->end(), 0.0, std::plus<double>());//eta = sum of alpha[0]
  for( unsigned i = 0; i < K; i++ ){
    mu[i] = (*alpha)[i]/eta;
  }
  
  double b = 0.0; // mu[K-1] + mu[0];//upper bound for sampler
  //double summu = 1.0 - mu[K-1];

  // loop over elements j,k of mu to update mu[j] conditional on (mu[j] + mu[k]), mu[i] where i neq j,k 
  for( unsigned int j = 1; j < K; ++j ) {
    for( unsigned int k = 0; k < j; ++k ) {
      b = mu[j] + mu[k]; 
      AlphaParameters[1] = eta; // dispersion parameter
      AlphaParameters[2] = b; 
      AlphaParameters[3] = sumlogtheta[j]; 
      AlphaParameters[4] = sumlogtheta[k];
      DirParamArray.setUpperBound(b-0.00001);
//       cout << "AlphaParameters " << AlphaParameters[0] << " " <<  AlphaParameters[1] << " " << 
// 	AlphaParameters[2] << " " <<  AlphaParameters[3] << " " <<  AlphaParameters[4] << " " << endl;
      try {
	mu[j] = DirParamArray.Sample(AlphaParameters, ddlogf);
      } catch(string msg) {
	cout << msg << endl;
	exit(1);
      }
      mu[k] = b - mu[j];
    }
    //cout << "mu" << j << " " << mu[j] << "\t";
  }
  //cout << endl << endl;

  //SampleEta((unsigned)AlphaParameters[0], sumlogtheta, &eta, mu);//first arg is num obs
  EtaArgs.sumlogtheta = sumlogtheta;
  EtaArgs.mu = mu;
  etanew = log(eta);//sample for log of dispersion parameter
  EtaSampler.Sample(&etanew, &EtaArgs);
  eta = exp(etanew);
  for( unsigned j = 0; j < K; j++ )
    (*alpha)[j] = mu[j]*eta;
  
  
#elif SAMPLERTYPE==2
  // *** Hamiltonian sampler for alpha
  AlphaArgs.sumlogtheta = sumlogtheta;
//   for (unsigned int k = 0; k < K; ++k) {
//     AlphaArgs.sumlogtheta[k] = sumlogtheta[k];
//   }
  transform(alpha[0].begin(), alpha[0].end(), logalpha, xlog);//logalpha = log(alpha)
  AlphaSampler.Sample(logalpha, &AlphaArgs);//sample new values for logalpha
  transform(logalpha, logalpha+K, alpha[0].begin(), xexp);//alpha = exp(logalpha)
#endif
  
}

#if SAMPLERTYPE==1
// void DirichletParamSampler::SampleEta(unsigned n, const double* const sumlogtheta, double *eta, const double* const mu){
//   // Dirichlet dispersion parameter eta is updated with a Metropolis random walk
//   unsigned int i;
//   double L1=0, P1=0, Proposal1=0;

//   etanew = exp( gennor( log( *eta ), step ) );
//   Proposal1 = log(etanew) - log(*eta);
//   // log prior ratio P1 
//   P1 = ( EtaAlpha - 1.0 ) * ( log(etanew) - log(*eta) ) - EtaBeta * ( etanew - *eta );
//   // log likelihood ratio L1
//   L1 = n * ( gsl_sf_lngamma( etanew ) - gsl_sf_lngamma( *eta ) );
//   for( i = 0; i < K; i++ )
//     L1 += mu[i] * (etanew - *eta) * sumlogtheta[i] - n*gsl_sf_lngamma( etanew * mu[i] ) + n*gsl_sf_lngamma( *eta * mu[i] );
//   // calculate log acceptance probability ratio
//   LogAccProb = 0.0;
//   if(P1 + L1 + Proposal1 < 0.0)  
//     LogAccProb = P1 + L1 + Proposal1; 
//   //accept/reject proposal
//   if( log(myrand()) < LogAccProb ){
//     *eta = etanew;
//   }
//   //update step size
//   step = TuneEta.UpdateStepSize( exp(LogAccProb) );
// }

double DirichletParamSampler::etaEnergy( const double* const x, const void* const vargs )
{
  const PopAdmixEtaSamplerArgs* args = (const PopAdmixEtaSamplerArgs* )vargs;
  double eta = exp(*x);
  double E = 0.0;
  const double*mu = args->mu;
  const double* sumlogtheta = args->sumlogtheta;

  // log prior  
  E += ( args->priorshape - 1.0 ) * ( log(eta) ) - args->priorrate *  eta ;
  // log likelihood 
  E += args->numobs * ( gsl_sf_lngamma( eta ) );
  for( unsigned i = 0; i < args->numpops; i++ )
    E += mu[i] * eta  * sumlogtheta[i] -  args->numobs * gsl_sf_lngamma( eta * mu[i] );
  //Jacobian
  E += eta;
  return -E;
}

void DirichletParamSampler::etaGradient( const double* const x, const void* const vargs, double* g )
{
  const PopAdmixEtaSamplerArgs* args = (const PopAdmixEtaSamplerArgs* )vargs;
  const double* mu = args->mu;
  const double* sumlogtheta = args->sumlogtheta;
  double eta = exp(*x);
  double psi;

  g[0] = 0;
  // log prior  
  g[0] -= ( args->priorshape - 1.0 ) / eta - args->priorrate;
  // log likelihood
  ddigam(&eta, &psi); 
  g[0] -= args->numobs * psi;//( gsl_sf_psi( eta ) );
  for( unsigned i = 0; i < args->numpops; i++ ){
    double alpha = eta*mu[i];
    ddigam(&alpha, &psi);
    g[0] -= mu[i] * sumlogtheta[i] -  args->numobs * mu[i] * psi;//gsl_sf_psi( eta * mu[i] );
  }
  //Jacobian
  g[0] *= eta;
}

// these 3 functions calculate log density and derivatives for adaptive rejection sampling of 
// a pair of elements of the proportion parameter of the Dirichlet distribution
double DirichletParamSampler::logf( double muj, const void* const pars ) {
  const double* parameters = (const double*) pars;
   int n = (int)parameters[0];
   double eta = parameters[1], b = parameters[2], sumlogpj = parameters[3], sumlogpk = parameters[4];
   if(muj < 0 || (b - muj < 0)) {
     throw string("\nDirichletParamSampler: negative argument to lngamma function\n");
   } 
   double f = eta * muj * ( sumlogpj - sumlogpk )
     - n * ( gsl_sf_lngamma(muj*eta) + gsl_sf_lngamma( (b - muj)*eta) );
   //cout << "\nlog density function passed muj " << muj<< "\treturns logdensity " << f << endl;
  return f;
}

double DirichletParamSampler::dlogf( double muj, const void* const pars ) {
  const double* parameters = (const double*) pars;
  double f, x1, x2, y1, y2;
  int n = (int)parameters[0];
  double eta = parameters[1], b = parameters[2], sumlogpj = parameters[3], sumlogpk = parameters[4];
  if(muj < 0 || (b - muj < 0)) {
    throw string("\nDirichletParamSampler: negative argument to digamma function\n");
  } 
  x1 = eta*muj;
  x2 = eta*(b - muj);
  ddigam( &x1, &y1 );
  ddigam( &x2, &y2 );
  f =  eta * ( sumlogpj - sumlogpk - n*( y1 - y2) );
  //cout << "\ngradient function passed muj "<< muj << "\treturns gradient " << f << flush;
  return f;
}

double DirichletParamSampler::ddlogf( double muj, const void* const pars) {
  const double* parameters = (const double*) pars;
  double f, x1, x2, y1, y2;
  int n = (int)parameters[0];
  double eta = parameters[1], b = parameters[2];
  x1 = eta*muj;
  x2 = eta*(b - muj);
  trigam( &x1, &y1 );
  trigam( &x2, &y2 );
  f = -n * eta * eta *( y1 + y2 );
  if(f >= 0) {
    throw string("DirichletParamSsampler: 2nd derivative non-negative\n");
  }
  return(f);
}

#elif SAMPLERTYPE==2

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

double DirichletParamSampler::getStepSize()const {
#if SAMPLERTYPE==1
  //return TuneEta.getStepSize();
  return EtaSampler.getStepsize();
#elif SAMPLERTYPE==2
    return AlphaSampler.getStepsize();
#endif
}

double DirichletParamSampler::getExpectedAcceptanceRate()const {
#if SAMPLERTYPE==1
  //return TuneEta.getExpectedAcceptanceRate();
  return EtaSampler.getAcceptanceRate();
#elif SAMPLERTYPE==2
    return AlphaSampler.getAcceptanceRate();
#endif
}

