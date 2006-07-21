#include "LogisticRegression.h"

LogisticRegression::LogisticRegression(){
  BetaSampler = 0;
  acceptbeta = 0;
  RegType = Logistic;
}
///deallocates arrays specific to this class
LogisticRegression::~LogisticRegression(){
  delete BetaSampler;
}

void LogisticRegression::Initialise(unsigned Number, double priorPrecision, const IndividualCollection* const individuals, LogWriter &Log){
  Regression::Initialise(Number, priorPrecision, individuals, Log);
 
  // ** Initialise Logistic Regression objects
    //Log << "\nGaussian priors on logistic regression parameters with precision " << lambda << "\n";

  std::vector<double> v = individuals->getOutcome(RegNumber);
  double p = accumulate(v.begin(), v.end(), 0.0, std::plus<double>()) / (double)v.size();
  //check the outcomes are not all 0s or all 1s
  if(p==0.0 || p==1.0)throw string("Data Error: All binary outcomes are the same");

  betamean[0] = log( p / ( 1 - p ) );
  beta[0] = betamean[0];
  
  //  ** initialize sampler for logistic regression **
  acceptbeta = 0;
  
  BetaParameters.n = NumIndividuals;
  BetaParameters.d = NumCovariates;
  BetaParameters.beta = beta;
  BetaSampler = new GaussianProposalMH( lr, dlr, ddlr);
}

void LogisticRegression::Update(bool sumbeta, IndividualCollection* individuals, double coolness
#ifdef PARALLEL
			, MPI::Intracomm &Comm){
  if(Comm.Get_rank() == 0){
#else 
  ){
#endif

  // Sample for regression model parameters beta
  std::vector<double> Outcome = individuals->getOutcome(RegNumber);

  Y = &(Outcome[0]);
  X = individuals->getCovariates();

  BetaParameters.Covariates = X; 
  BetaParameters.coolness = coolness; 
  matrix_product(Y, X, XtY, 1, NumIndividuals, NumCovariates);//XtY = X' * Y
  
  for( int j = 0; j < NumCovariates; j++ ){
    BetaParameters.beta0 = betamean[j];
    BetaParameters.priorprecision = betaprecision[j];
    BetaParameters.index = j;
    BetaParameters.XtY = XtY[ j ];
    acceptbeta = BetaSampler->Sample( beta + j, &BetaParameters );
  }
  
#ifdef PARALLEL
  }
  //broadcast parameters to workers
  Comm.Barrier();
  Comm.Bcast(beta, NumCovariates, MPI::DOUBLE, 0);
  if(RegType == Linear)Comm.Bcast(&lambda, 1, MPI::DOUBLE, 0);
#endif
  
  individuals->SetExpectedY(RegNumber,beta);
  
  if(sumbeta){
    SumParameters();
  }
}//end Update

void LogisticRegression::OutputParams(ostream* out){
  for( int j = 0; j < NumCovariates; j++ ){
    out->width(9);
    (*out) << setprecision(6) << beta[j] << "\t";
  }
}

double LogisticRegression::getDispersion()const{
  return 1.0;
}

double LogisticRegression::getLogLikelihood(const IndividualCollection* const IC)const{
  int NumIndividuals = IC->getSize();
  double loglikelihood = 0.0;

  //loglikelihood is sum of logs of bernoulli probabilities, given by EY
  for(int i = 0; i < NumIndividuals; ++i){
    if(IC->getOutcome(RegNumber, i))loglikelihood += log( IC->getExpectedY(i, RegNumber) );
    else loglikelihood += log(1.0 - IC->getExpectedY(i, RegNumber));
  }
  return loglikelihood;
}
double LogisticRegression::getLogLikelihoodAtPosteriorMeans(IndividualCollection *IC, int iterations){
  double logL = 0.0;

  //set expected outcome at posterior means of regression parameters
  for(int i = 0; i < NumCovariates; ++i)SumBeta[i] /= (double)iterations; 
  IC->SetExpectedY(RegNumber,SumBeta);//computes X * BetaBar
  for(int i = 0; i < NumCovariates; ++i)SumBeta[i] *= (double)iterations; //restore sumbeta

  //set precision to posterior mean
  double temp = lambda;
  lambda = SumLambda/(double)iterations;

  logL = getLogLikelihood(IC);
  lambda = temp;//restore precision
  IC->SetExpectedY(RegNumber,beta);//restore Xbeta

  return logL;
}

//(unnormalised)log posterior for a single regression parameter
double LogisticRegression::lr( const double beta, const void* const vargs )
{
  const BetaArgs* args = (const BetaArgs*)vargs;

  int n = args->n;
  int index = args->index ;
  double beta0 = 0;
  if( index == 0 )
    beta0 = args->beta0;
  double f = 0.0; 

  if(args->coolness > 0.0){

    //log likelihood contribution
    f += args->XtY * beta;

    double *Xbeta = new double[ n ];    
    Regression::getExpectedOutcome(args->beta, args->Covariates, Xbeta, n, args->d, index, beta);
    
    for( int i = 0; i < n; i++ ){
      f -= log( 1.0 + exp( Xbeta[ i ] ) );}
    
    //anneal likelihood
    f *= args->coolness;
    
    delete[] Xbeta;
  }

  f -= 0.5 * args->priorprecision * (beta - beta0) * (beta - beta0); //log prior contribution
  return( f );
}

//first and second derivatives of log posterior for a single regression parameter; 
//used in Newton-Raphson algorithm in MH sampler
double LogisticRegression::dlr( const double beta, const void* const vargs )
{
  const BetaArgs* args = (const BetaArgs*)vargs;

  int n = args->n;
  int d = args->d;
  int index = args->index ;
  double beta0 = 0;
  if( index == 0 )
    beta0 = args->beta0;
  double f = 0.0;

  if(args->coolness > 0.0){
    f += args->XtY;
    
    double *Xbeta = new double[ n ];
    
    Regression::getExpectedOutcome(args->beta, args->Covariates, Xbeta, n, d, index, beta);
    for( int i = 0; i < n; i++ )
      {
	f -= args->Covariates[ i*d + index ] / ( 1.0 + myexp( -Xbeta[ i ] ) );
      }
    delete[] Xbeta;
    f *= args->coolness;
  }

  f -= args->priorprecision * (beta - beta0);//log prior contribution
  return( f );
}

double LogisticRegression::ddlr( const double beta, const void* const vargs )
{
  const BetaArgs* args = (const BetaArgs*)vargs;

  int n = args->n;
  int d = args->d;
  int index = args->index ;

  double f = 0.0;

  if(args->coolness > 0.0){
    double *Xbeta = new double[ n ];
    
    Regression::getExpectedOutcome(args->beta, args->Covariates, Xbeta, n, d, index, beta);
    
    for( int i = 0; i < n; i++ )
      {
	f -= args->Covariates[ i*d + index ] * args->Covariates[ i*d + index ] / ( 2.0 + exp( -Xbeta[ i ] ) + exp( Xbeta[ i ] ) );
      }
    delete[] Xbeta;
    f *= args->coolness;
  }

  f -= args->priorprecision;//log prior contribution
  return( f );
}
