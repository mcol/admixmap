/** 
 *   LogisticRegression.cc 
 *   Class to represent and update parameters of a logistic regression model
 *   Copyright (c) 2006-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "bclib/LogisticRegression.h"
#include "bclib/linalg.h"
#include "bclib/misc.h" //for myexp
#include <math.h>//for sqrt and exp
#include <numeric>//for accumulate

using namespace::std;
BEGIN_BCLIB_NAMESPACE

LogisticRegression::LogisticRegression(unsigned Number, double priorPrecision, const bclib::DataMatrix& Covars, const bclib::DataMatrix& Outcome, 
				       LogWriter &Log): Regression(Number, Logistic){
  BetaSampler = 0;
  acceptbeta = 0;
  //RegType = Logistic;
  Initialise(priorPrecision, Covars, Outcome, Log);
}
///deallocates arrays specific to this class
LogisticRegression::~LogisticRegression(){
  delete BetaSampler;
}

void LogisticRegression::Initialise(double priorPrecision, const bclib::DataMatrix& Covars, const bclib::DataMatrix& Outcome, 
				    LogWriter &Log){
  Regression::Initialise(Covars.nCols(), Covars.nRows(), Covars.getData());
  Log.setDisplayMode(Quiet);
  std::vector<double> v = Outcome.getCol(RegNumber);
  double p = accumulate(v.begin(), v.end(), 0.0, std::plus<double>()) / (double)v.size();
  //check the outcomes are not all 0s or all 1s
  if(p==0.0 || p==1.0)throw string("Data Error: All binary outcomes are the same");

  betamean[0] = log( p / ( 1 - p ) );
  beta[0] = betamean[0];

  //set prior precision
  betaprecision[0] = priorPrecision;
  Log << "\nGaussian priors on logistic regression parameters with zero means and precisions\n ("<< betaprecision[0];
  
  for(int j = 1; j < NumCovariates; ++j){
    //get sample variance of covariate

      double sum = 0.0, sumsq = 0.0, x = 0.0;
      for(int i = 0; i < NumIndividuals; ++i){
	x = Covariates[i*NumCovariates + j];
	sum += x;
	sumsq += x*x;
      }
      double svc = (sumsq - sum*sum / (double)NumIndividuals) / (double)NumIndividuals;
      if (svc < 0.0001) svc = 1.0;// incase covar is sampled and begins constant or near constant over individuals

    betaprecision[j] = priorPrecision * svc;
    Log << ", " << betaprecision[j];
  }
  Log << ")\n";
  
  SetExpectedY(beta);

  //  ** initialize sampler for logistic regression **
  acceptbeta = 0;
  
  BetaParameters.n = NumIndividuals;
  BetaParameters.d = NumCovariates;
  BetaParameters.beta = beta;
  BetaSampler = new GaussianProposalMH( lr, dlr, ddlr);
}

void LogisticRegression::Update(bool sumbeta, const std::vector<double>& Outcome, double coolness  ){
  // Sample for regression model parameters beta

  BetaParameters.Covariates = Covariates; 
  BetaParameters.coolness = coolness; 
  matrix_product(&(Outcome[0]), Covariates, XtY, 1, NumIndividuals, NumCovariates);//XtY = X' * Y
  
  for( int j = 0; j < NumCovariates; j++ ){
    BetaParameters.beta0 = betamean[j];
    BetaParameters.priorprecision = betaprecision[j];
    BetaParameters.index = j;
    BetaParameters.XtY = XtY[ j ];
    acceptbeta = BetaSampler->Sample( beta + j, &BetaParameters );
  }
  
  if(sumbeta){
    SumParameters();
  }
  SetExpectedY(beta);
}//end Update

void LogisticRegression::SetExpectedY(const double* const beta){
  //sets ExpectedY = X * Beta
  if(ExpectedY){
    matrix_product(Covariates, beta, ExpectedY, NumIndividuals, NumCovariates, 1);
    //for binary outcome sets EY as logit^-1(X*beta)
    for(int i = 0; i < NumIndividuals; i++ )
      ExpectedY[i] = 1 / ( 1 + exp( -ExpectedY[i] ) );
  }
}

double LogisticRegression::getDispersion()const{
  return 1.0;
}
///returns Derivative of Inverse Link Function for individual i
double LogisticRegression::DerivativeInverseLinkFunction(unsigned i)const{
  return ExpectedY[i] * (1.0 - ExpectedY[i]);
}
double LogisticRegression::getLogLikelihood(const std::vector<double>& Outcome)const{
  double loglikelihood = 0.0;

  //loglikelihood is sum of logs of bernoulli probabilities, given by EY
  for(int i = 0; i < NumIndividuals; ++i){
    if(Outcome[ i ])loglikelihood += log( ExpectedY[ i ] );
    else loglikelihood += log(1.0 - ExpectedY[ i ]);
  }
  return loglikelihood;
}
double LogisticRegression::getLogLikelihoodAtPosteriorMeans(int iterations, const std::vector<double>& Outcome){
  double logL = 0.0;

  //set expected outcome at posterior means of regression parameters
  for(int i = 0; i < NumCovariates; ++i)SumBeta[i] /= (double)iterations; 
  //IC->SetExpectedY(RegNumber,SumBeta);//computes X * BetaBar
  SetExpectedY(SumBeta);
  for(int i = 0; i < NumCovariates; ++i)SumBeta[i] *= (double)iterations; //restore sumbeta

  logL = getLogLikelihood(Outcome);

  //IC->SetExpectedY(RegNumber,beta);//restore Xbeta
  SetExpectedY(beta);
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
	f -= args->Covariates[ i*d + index ] / ( 1.0 + eh_exp( -Xbeta[ i ] ) );
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
END_BCLIB_NAMESPACE
