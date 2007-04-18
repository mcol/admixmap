#include "bcppcl/CoxRegression.h"
#include <algorithm>
#include "bcppcl/misc.h"
#include "bcppcl/rand.h"
#include <math.h>

using namespace::std;
const double* CoxRegression::EY;

CoxRegression::CoxRegression(){
  BetaSampler = 0;
  acceptbeta = 0;
  RegType = Cox;
}

CoxRegression::~CoxRegression(){
  delete BetaSampler;
}

void CoxRegression::Initialise(unsigned Number, double priorPrecision, const DataMatrix& Covars, const DataMatrix& Outcome, 
			       LogWriter &Log){
  ReadData(Outcome);
  Regression::Initialise(Number, Covars.nCols()-1, Covars.nRows(), Covars.getData());

  //passing numcovariates-1 to exclude intercept
  // NumCovariates in this class is now one less

  //TODO: ?? rescale precisions somehow 
  Log.setDisplayMode(Quiet);
  Log << "\nGaussian priors on Cox regression parameters with zero means and precisions\n (";
  
  for(int j = 0; j < NumCovariates; ++j){
    betaprecision[j] = priorPrecision ;
    Log << ", " << betaprecision[j];
  }
  Log << ")\n";
  // set expected outcome
  //SetExpectedY(Covariates, beta, ExpectedY);  
  for(int i = 0; i < NumIndividuals; ++i){
    ExpectedY[i] = 0.0;
    for(int j = 0; j < NumCovariates; ++j)
      ExpectedY[i] += Covariates[i*(NumCovariates+1)+j+1]*beta[j]; 
    ExpectedY[i] = myexp(ExpectedY[i]);
  }
  EY = ExpectedY;

  //initialisation specific to Cox Regression
  c = 0.01;//TODO: decide sensible defaults for these, add options to allow user to set values
  mu = 1;

  //initialise hazard rates to their prior means
  for(int j = 0; j < BetaParameters.NumIntervals; ++j){
    BetaParameters.HazardRates.push_back( 1.0/*Rand::gengam(mu*c, c)*/ );
  }
  //  ** initialize sampler for regression params **
  acceptbeta = 0;
  
  BetaParameters.n = NumIndividuals;
  BetaParameters.d = NumCovariates;
  BetaParameters.beta = beta;
  BetaSampler = new GaussianProposalMH( lr, dlr, ddlr);

  plotloglikelihood(1, Covariates);
}

void CoxRegression::InitializeOutputFile(const std::vector<std::string>& CovariateLabels, unsigned NumOutcomes)
{
  // Header line of paramfile
  for( unsigned i = 0; i < CovariateLabels.size(); i++ ){
    outputstream << CovariateLabels[i] << "\t";
  }
  if(NumOutcomes == RegNumber+1)outputstream << std::endl;
}

///reads start times, end times and event counts from file.
///Call this function before initialise.
void CoxRegression::ReadData(const DataMatrix& CoxData){
  //data should contain three columns: startTime, endTime and event(s)
   std::vector<double> endpoints;

  for(unsigned i = 0; i < CoxData.nRows(); ++i){
    endpoints.push_back(CoxData.get(i,0));
    endpoints.push_back(CoxData.get(i,1));
    BetaParameters.events.push_back((unsigned)CoxData.get(i,2));
  }

  //determine endpoints of intervals
  sort(endpoints.begin(), endpoints.end());
  vector<double>::iterator new_end = unique(endpoints.begin(), endpoints.end());//reorder, putting duplicates at end
  endpoints.erase(new_end, endpoints.end());//remove duplicates
  BetaParameters.NumIntervals = endpoints.size()-1;
  for(int t = 0; t < BetaParameters.NumIntervals; ++t){
    BetaParameters.IntervalLengths.push_back(endpoints[t+1]-endpoints[t]);
  }

  for(unsigned i = 0; i < CoxData.nRows(); ++i){
    for(int t = 0; t < BetaParameters.NumIntervals; ++t){
      bool atrisk = (CoxData.get(i,0) <= endpoints[t]) && (CoxData.get(i,1) >= endpoints[ t+1]);
      BetaParameters.atRisk.push_back(atrisk);
      sum_nr.push_back(0.0);
      if(atrisk) sum_nr[t]+= BetaParameters.events[i];
    }
  }
}

void CoxRegression::Update(bool sumbeta, const std::vector<double>& , double coolness ){
  // Sample for regression model parameters beta
  BetaParameters.Covariates = Covariates;//X;
  BetaParameters.coolness = coolness; 
  
  try{
    for( int j = 0; j < NumCovariates; j++ ){
      double oldbeta = beta[j];
      BetaParameters.beta0 = betamean[j];
      BetaParameters.priorprecision = betaprecision[j];
      BetaParameters.index = j;
      acceptbeta = BetaSampler->Sample( beta + j, &BetaParameters );
      if(acceptbeta){
	for(int i = 0; i < NumIndividuals; ++i){
	  ExpectedY[i] *= myexp(Covariates[i*(NumCovariates+1)+j+1] * (beta[j] - oldbeta)) ;  
	}
      }
    }
  }
  catch(string s){
    throw string("Error encountered while sampling Cox regression parameters: " + s);
  }
  
  
//   //sample hazard rates
//   for(unsigned t = 0; t < BetaParameters.HazardRates.size(); ++t){
//     double shape = mu * c + sum_nr[t];
//     double rate = 0.0;
//     for(int i = 0; i < NumIndividuals; ++i)if(BetaParameters.atRisk[i*BetaParameters.NumIntervals +t])rate += ExpectedY[i];
//     rate *= BetaParameters.IntervalLengths[t];
//     rate += c;
    
//     BetaParameters.HazardRates[t] = Rand::gengam(shape, rate);
//   }
  //TODO: broadcast hazard rates

  if(sumbeta){
    SumParameters();
    //TODO:accumulate sums of hazard rates
  }
}

//void CoxRegression::SetExpectedY(const double* const Covariates, const double* const beta, double* ExpectedY){
  //doesn't really set expected outcome. stores Xbeta for use in updates
//   for(unsigned i = 0; i < NumIndividuals; ++i){
//     ExpectedY[i] = 0.0;
//     for(unsigned j = 0; j < NumCovariates; ++j)
//       ExpectedY[i] += Covariates[i*(NumCovariates+1)+j+1]*beta[j]; 
//     ExpectedY[i] = myexp(ExpectedY[i]);
//   }
//}

void CoxRegression::OutputParams(ostream* out)const{
  Regression::OutputParams(out);
  //TODO: ?? output hazard rates
}

double CoxRegression::getDispersion()const{
  return 1.0;//TODO: what is the dispersion in a Cox regression?
}
///returns Derivative of Inverse Link Function for individual i
double CoxRegression::DerivativeInverseLinkFunction(unsigned)const{
  //TODO
  return 1.0;    
}

///given an array of regression parameters beta and covariates X, computes expected outcome EY = X * beta, 
///with index'th element of beta replaced with betaj.
void CoxRegression::getExpectedOutcome(const double* const beta, const double* const X, double* EY, int n, int dim, int index, double betaj){
  //NOTE: here, X has dimension n* (dim+1 ) as it includes the intercept
  for(int i = 0; i < n; ++i){
    EY[i] = 0.0;
    for(int j = 0; j < dim; ++j)
      if(j == index)
	EY[i] += X[i*(dim+1)+j+1]*betaj;
      else
	EY[i] += X[i*(dim+1)+j+1]*beta[j];
  }
}

///given an array of regression parameters beta and covariates X, computes expected outcome EY = X * beta
void CoxRegression::getExpectedOutcome(const double* const beta, const double* const X, double* EY, int n, int d){
  getExpectedOutcome(beta, X, EY, n, d, -1, 0.0);
}

double CoxRegression::lr(const double beta, const void* const vargs){
  const CoxBetaArgs* args = (const CoxBetaArgs*)vargs;
  const int n = args->n;
  const int T = args->NumIntervals;
  const int index = args->index ;
  const double beta0 = args->beta0;
  double f = 0.0; 


  //log likelihood
  if(args->coolness > 0.0){
    try{
      for(int i = 0; i < n; ++i){
	double x = args->Covariates[i*(args->d+1) +index+1];
	double xbetai =  x * beta;
	double exp_xbetai = myexp(x*(beta - args->beta[index]) ) ;  
	int sum_nr = 0;
	double sum_nalpha = 0.0;


	for(int t = 0; t < T; ++t){
	  if( args->atRisk[i*T +t] ){
	    int r = args->events[i];
	    //if event occurred in this interval
	    if(!args->atRisk[i*T +t+1] || t==T-1)sum_nr += r;
	    sum_nalpha += args->HazardRates[t] * args->IntervalLengths[t];
	  }
	}
	//cout << "i = " << i << "exp Xbeta = " << exp_Xbetai << endl;
	f += xbetai*(double)sum_nr  - EY[i] * exp_xbetai * sum_nalpha; 
      }
    }
    catch(string s){
      throw string("error in lr: " +s );
    }
  }
  //log prior on beta
  f -= 0.5 * args->priorprecision * (beta - beta0) * (beta - beta0);
  return f;
}

double CoxRegression::dlr(const double beta, const void* const vargs){
  const CoxBetaArgs* args = (const CoxBetaArgs*)vargs;
  const int n = args->n;
  const int T = args->NumIntervals;
  const int index = args->index ;
  const double beta0 = args->beta0;
  double f = 0.0;

  //std::vector<bool>::const_iterator atrisk = args->atRisk.begin();
  //log likelihood
  if(args->coolness > 0.0){
    try{
      
      for(int i = 0; i < n; ++i){
	int sum_nr = 0;
	double sum_nalpha = 0.0;
	double exp_Xbetai = EY[i] * myexp(args->Covariates[i*(args->d+1)+index+1] * (beta - args->beta[index])) ;  
	//std::vector<double>::const_iterator alpha = args->HazardRates.begin();
	for(int t = 0; t < T; ++t){
	  if( args->atRisk[i*T +t] ){//individual i was at risk during this interval
	    if(!args->atRisk[i*T +t+1] || t==T-1)sum_nr += args->events[i];
	    sum_nalpha += args->HazardRates[t]* args->IntervalLengths[t]; 
	  }
	}
	f += ((double)sum_nr - exp_Xbetai * sum_nalpha) * args->Covariates[i*(args->d+1) +index+1];
      }
    }
    catch(string s){
      throw string("error in dlr: " +s );
    }

  }
  f -= args->priorprecision * (beta - beta0);//log prior contribution
  return f;
}

double CoxRegression::ddlr(const double beta, const void* const vargs){
  const CoxBetaArgs* args = (const CoxBetaArgs*)vargs;
  const int n = args->n;
  const int T = args->NumIntervals;
  const int index = args->index ;
  double f = 0.0; 

  //std::vector<bool>::const_iterator atrisk = args->atRisk.begin();
  //log likelihood
  if(args->coolness > 0.0){
    try{
      for(int i = 0; i < n; ++i){
	double sum = 0.0;
	double exp_Xbetai = EY[i] * myexp(args->Covariates[i*(args->d+1)+index+1] * (beta - args->beta[index])) ;  
	//std::vector<double>::const_iterator alpha = args->HazardRates.begin();
	for(int t = 0; t < T; ++t){
	  if( args->atRisk[i*T +t] ){
	    sum += args->HazardRates[t]* args->IntervalLengths[t]; 
	  }
	}
	double x = args->Covariates[i*(args->d+1) +index+1];
	f -= sum * x * x *exp_Xbetai;
      }
    }
    catch(string s){
      throw string("error in ddlr: " +s );
    }

  }
  f -= args->priorprecision;//log prior contribution
  return f;
}

void CoxRegression::plotloglikelihood(int j, const double* Covariates){
  ofstream outfile("loglfile.txt");

  BetaParameters.beta0 = betamean[j];
  BetaParameters.priorprecision = betaprecision[j];
  BetaParameters.index = j;
  BetaParameters.Covariates = Covariates;//X;
  BetaParameters.coolness = 1.0; 
  for(double betaj = -10.0; betaj < 10.0; ++betaj){

    outfile<< betaj << "\t" << lr(betaj, &BetaParameters) << "\t" << dlr(betaj, &BetaParameters) << "\t" << ddlr(betaj, &BetaParameters) << endl;

  }
  outfile.close();
}

///returns log likelihood at current parameter values
double CoxRegression::getLogLikelihood(const std::vector<double>& Outcome)const{
  return getLogLikelihood(beta, BetaParameters.HazardRates, Outcome);
}
///returns loglikelihood at supplied parameter values
double CoxRegression::getLogLikelihood(const double* const _beta, const std::vector<double>& _HazardRates, 
				       const std::vector<double>& )const{
  double L = 0.0;
  double *Xbeta = new double[ NumIndividuals ];    
  Regression::getExpectedOutcome(_beta, Covariates, Xbeta, NumIndividuals, NumCovariates);
  
  for(int i = 0; i < NumIndividuals; ++i){
    double xbetai =  Xbeta[i];
    double exp_xbetai = myexp(xbetai ) ;  
    int sum_nr = 0;
    double sum_nalpha = 0.0;
    
    const int T = BetaParameters.NumIntervals;
    for(int t = 0; t < T; ++t){
      if( BetaParameters.atRisk[i*T +t] ){
	int r = BetaParameters.events[i];
	//if event occurred in this interval
	if(!BetaParameters.atRisk[i*T +t+1] || t==T-1)sum_nr += r;
	//sum_nalpha += BetaParameters.HazardRates[t] * BetaParameters.IntervalLengths[t];
	sum_nalpha += _HazardRates[t] * BetaParameters.IntervalLengths[t];
      }
    }
    //cout << "i = " << i << "exp Xbeta = " << exp_Xbetai << endl;
    L += xbetai*(double)sum_nr  - exp_xbetai * sum_nalpha; 
  }
  delete[] Xbeta;
  return L;
}
///returns LogLikelihood at posterior means of parameters
double CoxRegression::getLogLikelihoodAtPosteriorMeans(int iterations, const std::vector<double>& Outcome){
  double logL = 0.0;

  //set expected outcome at posterior means of regression parameters
  for(int i = 0; i < NumCovariates; ++i)SumBeta[i] /= (double)iterations; 
  //TODO:set Hazard Rates to posterior means

  //compute loglikelihood at posterior means
  logL = getLogLikelihood(SumBeta, BetaParameters.HazardRates/*replace this with PMs of HazardRates*/, Outcome);

  for(int i = 0; i < NumCovariates; ++i)SumBeta[i] *= (double)iterations; //restore sumbeta


  return logL;
}

