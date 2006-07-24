/** 
 *   ADMIXMAP
 *   Regression.cc 
 *   Class to represent and update parameters of a regression model as used in admixmap
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "Regression.h"
#include "IndividualCollection.h"
#include <numeric>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

std::ofstream Regression::outputstream;

Regression::Regression(){
  lambda = 0; 
  SumLambda = 0;
  beta = 0;
  betamean = 0;
  betaprecision = 0;
  SumBeta = 0;
  NumOutcomeVars = 0;
  NumCovariates = 0;
  //RegType = None;
  XtY = 0;
}

Regression::~Regression(){
  delete[] XtY;
  delete[] beta;
  delete[] SumBeta;

  delete[] betamean;
  delete[] betaprecision;
}

void Regression::OpenOutputFile(const unsigned NumOutcomes, const char* const filename, LogWriter &Log){
  //Open paramfile
  if ( NumOutcomes){ 
    if ( strlen( filename ) ){
      outputstream.open( filename, ios::out );
      if( !outputstream )
	{
	  Log.setDisplayMode(On);
	  Log << "ERROR: Couldn't open regparamfile\n";
	  exit( 1 );
	}
      else{
	Log.setDisplayMode(Quiet);
	Log << "Writing regression parameters to " << filename << "\n";
      }
    }
    else{
      Log << "No regparamfile given\n";
    }
  }
}

void Regression::InitializeOutputFile(const std::vector<std::string>& CovariateLabels, unsigned NumOutcomes)
{
  // Header line of paramfile
  outputstream << "intercept\t";
  for( unsigned i = 0; i < CovariateLabels.size(); i++ ){
    outputstream << CovariateLabels[i] << "\t";
  }
  if(NumOutcomes == RegNumber+1)outputstream << endl;
}

void Regression::Initialise(unsigned Number, const unsigned numCovariates){
  //set regression number for this object
  RegNumber = Number;
  
  // ** Objects common to all regression types
  NumCovariates = numCovariates;
  beta = new double[ NumCovariates ];
  lambda = 1.0; 
}

void Regression::Initialise(unsigned Number, double priorPrecision, const IndividualCollection* const individuals, LogWriter &Log){
  Log.setDisplayMode(Quiet);
  //set regression number for this object
  RegNumber = Number;
  
  // ** Objects common to all regression types
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
  double outcomeSampleVariance = individuals->getSampleVarianceOfOutcome(RegNumber);
  betaprecision[0] = priorPrecision / outcomeSampleVariance;
  Log << "\nGaussian priors on " << RegressionString[(int)RegType] 
      << " regression parameters with zero means and precisions\n ("<< betaprecision[0];
  
  for(int j = 1; j < NumCovariates; ++j){
    betaprecision[j] = priorPrecision * individuals->getSampleVarianceOfCovariate(j) / outcomeSampleVariance;
    Log << ", " << betaprecision[j];
  }
  Log << ")\n";
  
  //X = individuals->getCovariates();
  XtY = new double[NumCovariates];
  
  lambda = 1.0; 
  SumLambda = lambda;
  
}

///given an array of regression parameters beta and covariates X, computes expected outcome EY = X * beta, 
///with index'th element of beta replaced with betaj.
void Regression::getExpectedOutcome(const double* const beta, const double* const X, double* EY, int n, int dim, int index, double betaj){
  double* beta1 = new double[dim];

  for( int j = 0; j < dim; j++ )
    {
      if( j != index )
	beta1[ j ] = beta[j];
      else
	beta1[ j ] = betaj;//substitute supplied betaj for current value of beta[j]
    }
  //Xbeta = X * beta1;
  matrix_product(X, beta1, EY, n, dim, 1);
  delete[] beta1;
}

///given an array of regression parameters beta and covariates X, computes expected outcome EY = X * beta
void Regression::getExpectedOutcome(const double* const beta, const double* const X, double* EY, int n, int d){
  getExpectedOutcome(beta, X, EY, n, d, -1, 0.0);
}

void Regression::OutputParams(ostream* out)const{
  for( int j = 0; j < NumCovariates; j++ ){
    out->width(9);
    (*out) << setprecision(6) << beta[j] << "\t";
  }
}

void Regression::Output(const unsigned NumberOfOutcomes, bool toScreen, bool afterBurnIn){
  //output to screen
  if( toScreen )
    {
      if(NumberOfOutcomes>1)cout << "\nRegression " << RegNumber << "\t";
      OutputParams(&cout);
      cout << endl;
    }
  //Output to paramfile after BurnIn
  if( afterBurnIn ){
    OutputParams(&outputstream);
    if(RegNumber==NumberOfOutcomes-1) outputstream << endl;
    //output new line in paramfile when last regression model
  }
}

void Regression::OutputErgodicAvg(int samples, std::ofstream *avgstream)const{
 //output to ergodicaveragefile
  cout << "in base function" << endl;
  for( int j = 0; j < NumCovariates; j++ ){
    avgstream->width(9);
    *avgstream << setprecision(6) << SumBeta[j] / samples << "\t";
  }
}

void Regression::SumParameters(){
  // accumulate sum of parameters after burnin.
  if( NumCovariates > 0 )
      for(int j = 0; j < NumCovariates; ++j)
        SumBeta[j] += beta[j];
  SumLambda += lambda;
}
const double* Regression::getbeta()const{
  return beta;
}
double Regression::getlambda()const{
  return lambda;
}
int Regression::getNumCovariates()const{
  return NumCovariates;
}



