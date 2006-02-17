#include "DirichletParamSampler.h"
#include <algorithm>
#include <numeric>
#include "functions.h"
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_psi.h>
#include<gsl/gsl_sf_result.h>

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
  mu = 0;
  munew = 0;
  muDirichletParams = 0;
#elif SAMPLERTYPE==2
  logalpha = 0;
  initialAlphaStepsize = 0.02;//need a way of setting this without recompiling, or a sensible fixed value
  targetAlphaAcceptRate = 0.7;// use high value with many leapfrog steps
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
   DirParamArray.Initialise(true, true, 1.0, 0.0, logf, dlogf); // avoid singularity at mu[j]=0

   EtaArgs.priorshape = K; // for compatibility with gamma(1, 1) prior on alpha
   EtaArgs.priorrate = 1;
   EtaArgs.numpops = numpops;
   EtaArgs.numobs = numobs; 
   // use many small steps to ensure that leapfrog algorithm does not jump to minus infinity
   EtaSampler.SetDimensions(1, 0.005/*initial stepsize*/, 0.001/*min stepsize*/, 1.0/*max stepsize*/, 100/*num leapfrogs*/,
			    0.8/*target accept rate*/, etaEnergy, etaGradient);
#elif SAMPLERTYPE==2
   logalpha = new double[K];
   AlphaArgs.n = numobs; //num individuals/gametes will be passed as arg to sampler
   AlphaArgs.dim = K;
   AlphaArgs.eps0 = 1.0; //Gamma(1, 1) prior on alpha
   AlphaArgs.eps1 = 1.0;
   AlphaSampler.SetDimensions(K, initialAlphaStepsize, 0.01, 100.0, 50, targetAlphaAcceptRate, findE, gradE);
#endif
}

DirichletParamSampler::~DirichletParamSampler()
{
#if SAMPLERTYPE==1
  delete[] mu;
  delete[] munew;
  delete[] muDirichletParams;
#elif SAMPLERTYPE==2
  delete[] logalpha;
#endif
}

#if SAMPLERTYPE==1
void DirichletParamSampler::SetPriorEta( double inEtaAlpha, double inEtaBeta ) {
   EtaArgs.priorshape = inEtaAlpha; 
   EtaArgs.priorrate = inEtaBeta;
}
void DirichletParamSampler::SetPriorMu( const double* const ingamma ) {
   for( unsigned int i = 0; i < K; i++ ){
      muDirichletParams[i] = ingamma[i];
   }
}
#endif

void DirichletParamSampler::Sample( const double* const sumlogtheta, std::vector<double> *alpha ) {
  // sumlogtheta = sum log observed proportions
  // update elements of mu with adaptive rejection sampler conditional on sumlogtheta
  // update eta with Hamiltonian  
#if SAMPLERTYPE==1
  gsl_set_error_handler_off();
  eta = accumulate(alpha->begin(), alpha->end(), 0.0, std::plus<double>());//eta = sum of alpha[0]
  for( unsigned i = 0; i < K; i++ ) {
    mu[i] = (*alpha)[i]/eta;
    //cout << (*alpha)[i] << " ";
  }
  double b = 0.0; 

  for(int updates=0; updates < 2; ++ updates) { // loop twice 
    // loop over elements j,k of mu to update mu[j] conditional on (mu[j] + mu[k]), mu[i] where i neq j,k 
    for( unsigned int j = 1; j < K; ++j ) {
      for( unsigned int k = 0; k < j; ++k ) {
	b = mu[j] + mu[k]; 
	AlphaParameters[1] = eta; // dispersion parameter
	AlphaParameters[2] = b; 
	AlphaParameters[3] = sumlogtheta[j]; 
	AlphaParameters[4] = sumlogtheta[k];
	DirParamArray.setUpperBound(b); // avoid singularity at b
	try {
	  mu[j] = DirParamArray.Sample(AlphaParameters, ddlogf);
	} catch(string msg) {
	  cout << msg << endl;
	  exit(1);
	}
	mu[k] = b - mu[j];
      }
    }
  }

  //SampleEta((unsigned)AlphaParameters[0], sumlogtheta, &eta, mu);//first arg is num obs
  EtaArgs.sumlogtheta = sumlogtheta;
  EtaArgs.mu = mu;
  // cout << "eta " << eta << endl << flush;
  etanew = log(eta);//sample for log of dispersion parameter
  try {
    EtaSampler.Sample(&etanew, &EtaArgs);
  } catch(string msg) {
    cout << msg << endl;
    exit(1);
  }
  eta = exp(etanew);
  for( unsigned j = 0; j < K; j++ ) (*alpha)[j] = mu[j]*eta;
  
  
#elif SAMPLERTYPE==2
  // *** Hamiltonian sampler for alpha
  AlphaArgs.sumlogtheta = sumlogtheta;
  transform(alpha[0].begin(), alpha[0].end(), logalpha, xlog);//logalpha = log(alpha)
  AlphaSampler.Sample(logalpha, &AlphaArgs);//sample new values for logalpha
  transform(logalpha, logalpha+K, alpha[0].begin(), xexp);//alpha = exp(logalpha)
#endif
  
}

#if SAMPLERTYPE==1

double DirichletParamSampler::etaEnergy( const double* const x, const void* const vargs )
{
  const PopAdmixEtaSamplerArgs* args = (const PopAdmixEtaSamplerArgs* )vargs;
  double eta = exp(*x);
  double E = 0.0;
  const double*mu = args->mu;
  const double* sumlogtheta = args->sumlogtheta;

  // log prior (on log eta scale)
  E += ( args->priorshape ) * ( log(eta) ) - args->priorrate *  eta ;
  // log likelihood 
  int status = 0;
  gsl_sf_result lngamma_result;
  status = gsl_sf_lngamma_e( eta, &lngamma_result );
  if(status) throw string("gsl lngamma error\n");
  E += args->numobs * lngamma_result.val;
  for( unsigned i = 0; i < args->numpops; ++i ) {
    status = gsl_sf_lngamma_e( eta*mu[i], &lngamma_result );
    if(status) throw string("gsl lngamma error\n");
    E += mu[i] * eta  * sumlogtheta[i] -  args->numobs * lngamma_result.val; // gsl_sf_lngamma( eta * mu[i] );
  }

  return -E;
}

void DirichletParamSampler::etaGradient( const double* const x, const void* const vargs, double* g )
{
  const PopAdmixEtaSamplerArgs* args = (const PopAdmixEtaSamplerArgs* )vargs;
  const double* mu = args->mu;
  const double* sumlogtheta = args->sumlogtheta;
  double eta = exp(*x);
  //double psi;

  g[0] = 0;
  // log prior  
  g[0] -= ( args->priorshape ) / eta - args->priorrate;
  // log likelihood
  int status = 0;
  gsl_sf_result psi_result;
  status = gsl_sf_psi_e(eta, &psi_result);
  if(status) {
    cout << "\n" << eta <<" as argument eta to digamma function";
    throw string("\nERROR in etaGradient: gsl digamma function\n");
  }
  g[0] -= args->numobs * psi_result.val;//( gsl_sf_psi( eta ) );
  for( unsigned i = 0; i < args->numpops; ++i ){
    //double alpha = eta*mu[i];
    status = gsl_sf_psi_e(eta*mu[i], &psi_result);
    if(status) {
      cout << "\n" << eta*mu[i] <<" as argument eta*mu to digamma function";
      throw string("\nERROR in etaGradient: gsl digamma error\n");
    }
    g[0] -= mu[i] * sumlogtheta[i] -  args->numobs * mu[i] * psi_result.val;//gsl_sf_psi( eta * mu[i] );
  }
  //use chain rule 
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
   int status = 0;
   double y1, y2;
   gsl_sf_result lngamma_result;
   status = gsl_sf_lngamma_e( muj*eta, &lngamma_result );
   if(status) throw string("\nERROR in DirichletParamSampler::logf - gsl lngamma error\n");
   y1 = lngamma_result.val;
   status = gsl_sf_lngamma_e( (b-muj)*eta, &lngamma_result );
   if(status) throw string("\nERROR in DirichletParamSampler::logf - gsl lngamma error\n");
   y2 = lngamma_result.val;
   double f = eta * muj * ( sumlogpj - sumlogpk ) - n * ( y1 + y2 );
   //cout << "\nlog density function passed muj " << muj << "\treturns logdensity " << f << endl;
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
  int status = 0;
  gsl_sf_result psi_result;
  status = gsl_sf_psi_e(x1, &psi_result);
  if(status) throw string("gsl digamma error\n");
  y1 = psi_result.val;
  status = gsl_sf_psi_e(x2, &psi_result);
  if(status) throw string("gsl digamma error\n");
  y2 = psi_result.val;
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
  int status = 0;
  gsl_sf_result psi1_result;
  status = gsl_sf_psi_n_e(1, x1, &psi1_result);
  if(status) throw string("gsl trigamma error\n");
  y1 = psi1_result.val;
  status = gsl_sf_psi_n_e(1, x2, &psi1_result);
  if(status) throw string("gsl trigamma error\n");
  y2 = psi1_result.val;
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

