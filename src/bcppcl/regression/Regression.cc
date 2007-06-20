/** 
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
#include "bcppcl/Regression.h"
#include "bcppcl/linalg.h"
#include <iomanip>
#include <algorithm>
#include <vector>
#include <string>

using std::vector;
using std::string;

std::ofstream Regression::outputstream;
RObjectWriter Regression::EYStream;
int Regression::NumIndividuals;
int Regression::NumOutcomeVars;

Regression::Regression(){
  lambda = 0; 
  SumLambda = 0;
  beta = 0;
  betamean = 0;
  betaprecision = 0;
  SumBeta = 0;
  NumOutcomeVars = 0;
  NumCovariates = 0;
  Outcome = 0;
  Covariates = 0;
  ExpectedY = 0;
  //RegType = None;
  XtY = 0;
}

Regression::~Regression(){
  delete[] XtY;
  delete[] ExpectedY;
  delete[] beta;
  delete[] SumBeta;

  delete[] betamean;
  delete[] betaprecision;
}

void Regression::OpenOutputFile(const unsigned NumOutcomes, const char* const filename, LogWriter &Log){
  //Open paramfile
  NumOutcomeVars = NumOutcomes;
  if ( NumOutcomes){ 
    if ( strlen( filename ) ){
      outputstream.open( filename, std::ios::out );
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
  if(NumOutcomes == RegNumber+1)outputstream << std::endl;
}

void Regression::OpenExpectedYFile(const char* Filename, LogWriter & Log){
  try{
    EYStream.open(Filename);
  }
  catch(...){
    Log.setDisplayMode(On);
    Log<< "WARNING: Couldn't open expectedoutcomefile\n";
  }
  Log.setDisplayMode(Quiet);
  Log << "Writing expected values of outcome variable(s) to " << Filename << "\n";
}

void Regression::Initialise(unsigned Number, const unsigned numCovariates){
  //set regression number for this object
  RegNumber = Number;
  
  // ** Objects common to all regression types
  NumCovariates = numCovariates;
  beta = new double[ NumCovariates ];
  lambda = 1.0; 
}

void Regression::Initialise(unsigned Number, unsigned nCovariates, unsigned nIndivs, const double* const Covars){
  //set regression number for this object
  RegNumber = Number;
  
  // ** Objects common to all regression types
  NumCovariates = nCovariates;
  NumIndividuals = nIndivs;
  
  beta = new double[ NumCovariates ];
  SumBeta = new double[ NumCovariates ];
  std::fill(beta, beta + NumCovariates, 0.0);
  std::fill(SumBeta, SumBeta + NumCovariates, 0.0);
  
  betamean = new double[ NumCovariates ];
  std::fill(betamean, betamean + NumCovariates, 0.0);
  
  //initialise regression params at prior mean
  for(int j = 0; j < NumCovariates; ++j){
    beta[j] = betamean[j];
  }

  betaprecision = new double[NumCovariates];
  
  Covariates = Covars;
  XtY = new double[NumCovariates];
  ExpectedY = new double[NumIndividuals];
  
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

void Regression::OutputParams(std::ostream* out)const{
  for( int j = 0; j < NumCovariates; j++ ){
    out->width(9);
    (*out) << std::setprecision(6) << beta[j] << "\t";
  }
}

void Regression::Output(const unsigned NumberOfOutcomes, bool toScreen, bool afterBurnIn){
  //output to screen
  if( toScreen )
    {
      if(NumberOfOutcomes>1)std::cout << "\nRegression " << RegNumber << "\t";
      OutputParams(&std::cout);
      std::cout << std::endl;
    }
  //Output to paramfile after BurnIn
  if( afterBurnIn ){
    OutputParams(&outputstream);
    if(RegNumber==NumberOfOutcomes-1) outputstream << std::endl;
    //output new line in paramfile when last regression model
  }
}

void Regression::OutputExpectedY(){
  //output kth Expected Outcome to file
  if(EYStream.is_open()){
    for(int i = 0; i < NumIndividuals; ++i)
      EYStream << ExpectedY[i];
    EYStream << newline;
  }  
}

///finish writing expected outcome as R object
void Regression::FinishWritingEYAsRObject(unsigned NumIterations, const std::vector<std::string>& Labels){
  //dimensions are NumIndividuals, NumOutcomes, NumIterations
  vector<int> dims(3);
  dims[0] = NumIndividuals;
  dims[1] = Labels.size();
  dims[2] = NumIterations;

  vector<vector<string> > dimnames(3);
  for(unsigned j = 0; j < Labels.size(); ++j){
    dimnames[1].push_back( Labels[j]);
  }
  if(EYStream.is_open()){
    EYStream.close(dims, dimnames);
  }
}

void Regression::OutputErgodicAvg(int samples, std::ofstream *avgstream)const{
 //output to ergodicaveragefile
  for( int j = 0; j < NumCovariates; j++ ){
    avgstream->width(9);
    *avgstream << std::setprecision(6) << SumBeta[j] / samples << "\t";
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
const double* Regression::getExpectedOutcome()const{
  return ExpectedY;
}
double Regression::getExpectedOutcome(unsigned i)const{
  return ExpectedY[i];
}
RegressionType Regression::getRegressionType()const{
  return RegType;
}

