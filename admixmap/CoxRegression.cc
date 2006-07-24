#include "CoxRegression.h"
#include <algorithm>

using namespace::std;

CoxRegression::CoxRegression(){
  BetaSampler = 0;
  acceptbeta = 0;
  RegType = Cox;
}

CoxRegression::CoxRegression(const DataMatrix& CoxData){
  CoxRegression();
  ReadData(CoxData);
}

CoxRegression::~CoxRegression(){
  delete BetaSampler;
}

void CoxRegression::Initialise(unsigned Number, double priorPrecision, const IndividualCollection* const individuals, LogWriter &Log){
  Log.setDisplayMode(Quiet);
  //set regression number for this object
  RegNumber = Number;
  
  // ** Objects common to both regression types
  NumCovariates = individuals->GetNumCovariates();
  NumIndividuals = individuals->getSize();
  
  beta = new double[ NumCovariates ];
  SumBeta = new double[ NumCovariates ];
  fill(beta, beta + NumCovariates, 0.0);
  fill(SumBeta, SumBeta + NumCovariates, 0.0);
  
  betamean = new double[ NumCovariates ];
  fill(betamean, betamean + NumCovariates, 0.0);
  
  //initialise regression params at prior mean
  for(int j = 0; j < NumCovariates; ++j){
    beta[j] = betamean[j];
  }

  betaprecision = new double[NumCovariates];
  //TODO: ?? rescale precisions somehow 
  //double outcomeSampleVariance = individuals->getSampleVarianceOfOutcome(RegNumber);
  betaprecision[0] = priorPrecision;// / outcomeSampleVariance;
  Log << "\nGaussian priors on Cox regression parameters with zero means and precisions\n ("<< betaprecision[0];
  
  for(int j = 1; j < NumCovariates; ++j){
    betaprecision[j] = priorPrecision ;//* individuals->getSampleVarianceOfCovariate(j) / outcomeSampleVariance;
    Log << ", " << betaprecision[j];
  }
  Log << ")\n";
  
  //X = individuals->getCovariates();
  XtY = new double[NumCovariates];
  
  
  //initialisation specific to Cox Regression
  BetaParameters.c = 0.01;//TODO: decide sensible defaults for these, add options to allow user to set values
  BetaParameters.mu = 1;
  //num intervals = numindividuals
  //initialise hazard rates to their prior means
  for(int j = 0; j < BetaParameters.NumIntervals; ++j)
    HazardRates.push_back( 1.0);//BetaParameters.mu*intervalLength(j) );
  
  //  ** initialize sampler for regression params **
  acceptbeta = 0;
  
  BetaParameters.n = NumIndividuals;
  BetaParameters.d = NumCovariates;
  BetaParameters.startTimes = startTimes.begin();
  BetaParameters.endTimes = endTimes.begin();
  BetaParameters.events = events.begin();
  BetaParameters.endpoints = endpoints.begin();
  BetaParameters.HazardRates = HazardRates.begin();
  BetaParameters.beta = beta;
  BetaSampler = new GaussianProposalMH( lr, dlr, ddlr);
}

///reads start times, end times and event counts from file.
///Call this function before initialise.
void CoxRegression::ReadData(const DataMatrix& CoxData){
  //data should contain three columns: startTime, endTime and event(s)
  for(unsigned i = 1; i < CoxData.nRows(); ++i){
    startTimes.push_back((int)CoxData.get(i,0));
    endTimes.push_back((int)CoxData.get(i,1));
    events.push_back((unsigned)CoxData.get(i,2));
  }
  //determine endpoints of intervals
  endpoints.resize(startTimes.size() + endTimes.size());
  copy(startTimes.begin(), startTimes.end(), endpoints.begin());
  copy(endTimes.begin(), endTimes.end(), endpoints.begin()+startTimes.size());
  sort(endpoints.begin(), endpoints.end());
  vector<int>::iterator new_end = unique(endpoints.begin(), endpoints.end());//reorder, putting duplicates at end
  endpoints.erase(new_end, endpoints.end());//remove duplicates
  BetaParameters.NumIntervals = endpoints.size()-1;
  
}

bool CoxRegression::atRisk(unsigned ind, unsigned interval)const{
  return (bool) (startTimes[ind] <= endpoints[interval] && endTimes[ind] >= endpoints[ interval+1]);
}
unsigned CoxRegression::intervalLength(unsigned t)const{
  if(t > endpoints.size()) throw string("ERROR: invalid interval index in Cox Regression");
  return endpoints[t+1] - endpoints[t];
}
unsigned CoxRegression::numFailures(unsigned ind, unsigned interval)const{
  unsigned num = 0;
  if( atRisk(ind, interval)) num = events[ind]; 
  return num;
}
bool CoxRegression::atRisk( const vector<int>::const_iterator start, 
			    const vector<int>::const_iterator finish, const vector<int>::const_iterator endpts){
  return (bool) ( (*start <= *endpts) && (*finish >= *(endpts + 1) ));
}

void CoxRegression::Update(bool sumbeta, const std::vector<double>& , const double* const Covariates, double coolness
#ifdef PARALLEL
			   , MPI::Intracomm &Comm
#endif
			   ){
  // Sample for regression model parameters beta
  BetaParameters.Covariates = Covariates;//X;
  BetaParameters.coolness = coolness; 
  
  try{
    for( int j = 0; j < NumCovariates; j++ ){
      BetaParameters.beta0 = betamean[j];
      BetaParameters.priorprecision = betaprecision[j];
      BetaParameters.index = j;
      acceptbeta = BetaSampler->Sample( beta + j, &BetaParameters );
    }
  }
  catch(string s){
    throw string("Error encountered while sampling Cox regression parameters: " + s);
  }
  
#ifdef PARALLEL
  //broadcast parameters to workers
  Comm.Barrier();
  Comm.Bcast(beta, NumCovariates, MPI::DOUBLE, 0);
#endif
  
  //TODO: sample hazard rates and broadcast
//   for(vector<double>::forward_iterator t = HazardRates.begin(); t != HazardRates.end(); ++t){
//     *t = gengam(
//     }


  if(sumbeta){
    SumParameters();
    //TODO:accumulate sums of hazard rates
  }
  
}

void CoxRegression::OutputParams(ostream* out){
  Regression::OutputParams(out);
  //TODO: ?? output hazard rates
}

double CoxRegression::getDispersion()const{
  return 1.0;//TODO: what is the dispersion in a Cox regression?
}

double CoxRegression::lr(const double beta, const void* const vargs){
  const CoxBetaArgs* args = (const CoxBetaArgs*)vargs;
  const int n = args->n;
  const int T = args->NumIntervals;
  const int index = args->index ;
  const double beta0 = args->beta0;
  //const double c = args->c;
  double f = 0.0; 

  const vector<int>::const_iterator firstStart = args->startTimes;
  const vector<int>::const_iterator firstEnd = args->endTimes;
  const vector<int>::const_iterator firstEndpoint = args->endpoints;
  const vector<double>::const_iterator hazardRates = args->HazardRates;

  //log likelihood
  if(args->coolness > 0.0){
    double *Xbeta = new double[ n ];    
    try{
      Regression::getExpectedOutcome(args->beta, args->Covariates, Xbeta, n, args->d, index, beta);
      for(int i = 0; i < n; ++i)
	for(int t = 0; t < T; ++t){
	  if( atRisk(firstStart+i, firstEnd+i, firstEndpoint+t) ){
	    const double lambda = *(hazardRates + t)*myexp(Xbeta[i]);
	    int r = (*args->events +i);
	    f += (double)r * mylog(lambda) - lambda; 
	  }
	}
    }
    catch(string s){
      throw string("error in lr: " +s );
    }
    delete[] Xbeta;
  }
  //log prior on beta
  f -= 0.5 * args->priorprecision * (beta - beta0) * (beta - beta0);
  //log prior on hazard rates
  //   for(int t = 0; t < T; ++t){
  //     f += (args->mu * (*(firstEndpoint+t+1) - *(firstEndpoint+t)) * c - 1.0)*mylog(*(hazardRates +t)) - c *  *(hazardRates+t);
  //   }
  //TODO: loops in loglikelihood could be rearranged to include log prior on hazard rates
  return f;
}

double CoxRegression::dlr(const double beta, const void* const vargs){
  const CoxBetaArgs* args = (const CoxBetaArgs*)vargs;
  const int n = args->n;
  const int T = args->NumIntervals;
  const int index = args->index ;
  const double beta0 = args->beta0;
  //const double c = args->c;
  double f = 0.0;

  const vector<int>::const_iterator firstStart = args->startTimes;
  const vector<int>::const_iterator firstEnd = args->endTimes;
  const vector<int>::const_iterator firstEndpoint = args->endpoints;
  const vector<double>::const_iterator hazardRates = args->HazardRates;

  //log likelihood
  if(args->coolness > 0.0){
    double *Xbeta = new double[ n ];    
    try{
      Regression::getExpectedOutcome(args->beta, args->Covariates, Xbeta, n, args->d, index, beta);
      
      for(int i = 0; i < n; ++i){
	double sum = 0.0;
	for(int t = 0; t < T; ++t){
	  if( atRisk(firstStart+i, firstEnd+i, firstEndpoint+t) ){
	    const double lambda = *(args->HazardRates + t)*myexp(Xbeta[i]);
	    int r = (*args->events +i);
	    sum += (double)r - lambda; 
	  }
	}
	f += sum * args->Covariates[i*args->d +index];
      }
    }
    catch(string s){
      throw string("error in dlr: " +s );
    }

    delete[] Xbeta;
  }
  f -= args->priorprecision * (beta - beta0);//log prior contribution
//   for(int t = 0; t < T; ++t){
//     f += (args->mu * (*(firstEndpoint+t+1) - *(firstEndpoint+t)) * c - 1.0) / (*(hazardRates +t)) - c ;
//   }
  return f;
}

double CoxRegression::ddlr(const double beta, const void* const vargs){
  const CoxBetaArgs* args = (const CoxBetaArgs*)vargs;
  const int n = args->n;
  const int T = args->NumIntervals;
  const int index = args->index ;
  //const double c = args->c;
  double f = 0.0; 

  const vector<int>::const_iterator firstStart = args->startTimes;
  const vector<int>::const_iterator firstEnd = args->endTimes;
  const vector<int>::const_iterator firstEndpoint = args->endpoints;
  const vector<double>::const_iterator hazardRates = args->HazardRates;

  //log likelihood
  if(args->coolness > 0.0){
    double *Xbeta = new double[ n ]; 
    try{
      Regression::getExpectedOutcome(args->beta, args->Covariates, Xbeta, n, args->d, index, beta);
      
      for(int i = 0; i < n; ++i){
	double sum = 0.0;
	for(int t = 0; t < T; ++t){
	  if( atRisk(firstStart+i, firstEnd+i, firstEndpoint+t) ){
	    sum += myexp(Xbeta[i]); 
	  }
	}
	double x = args->Covariates[i*args->d +index];
	f -= sum * x * x;
      }
    }
    catch(string s){
      throw string("error in ddlr: " +s );
    }

    delete[] Xbeta;
  }
  f -= args->priorprecision;//log prior contribution
//   for(int t = 0; t < T; ++t){
//     f += (args->mu * (1.0 - *(firstEndpoint+t+1) - *(firstEndpoint+t)) * c ) / ( (*(hazardRates +t)) * (*(hazardRates +t)) );
//   }
  return f;
}

///returns log likelihood at current parameter values
double CoxRegression::getLogLikelihood(const IndividualCollection* const IC)const{
  return getLogLikelihood(beta, HazardRates, IC);
}
///returns loglikelihood at supplied parameter values
double CoxRegression::getLogLikelihood(const double* const _beta, const std::vector<double>& _HazardRates, 
				       const IndividualCollection* const IC)const{
  double L = 0.0;
  double *Xbeta = new double[ NumIndividuals ];    
  Regression::getExpectedOutcome(_beta, IC->getCovariates(), Xbeta, NumIndividuals, NumCovariates);
  
  for(int i = 0; i < NumIndividuals; ++i)
    for(int t = 0; t < BetaParameters.NumIntervals; ++t){
      if( atRisk(i, t) ){
	const double lambda = _HazardRates[t]*myexp(Xbeta[i]);
	int r = numFailures(i, t);
	L += r*mylog(lambda) - lambda; 
      }
    }
  delete[] Xbeta;
  return L;
}
///returns LogLikelihood at posterior means of parameters
double CoxRegression::getLogLikelihoodAtPosteriorMeans(IndividualCollection *IC, int iterations){
  //under construction!
  double logL = 0.0;

  //set expected outcome at posterior means of regression parameters
  for(int i = 0; i < NumCovariates; ++i)SumBeta[i] /= (double)iterations; 
  //TODO:set Hazard Rates to posterior means

  //compute loglikelihood at posterior means
  logL = getLogLikelihood(SumBeta, HazardRates/*replace this with PMs of HazardRates*/, IC);

  for(int i = 0; i < NumCovariates; ++i)SumBeta[i] *= (double)iterations; //restore sumbeta


  return logL;
}
