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
  X = 0;
  Y = 0;
  XtY = 0;
}

Regression::~Regression(){
  delete[] XtY;
  delete[] beta;
  delete[] SumBeta;

  delete[] betamean;
  delete[] betaprecision;
}

void Regression::OpenOutputFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
				const std::string* const PopulationLabels, LogWriter &Log){
  //Open paramfile
  if ( options->getNumberOfOutcomes()>0){ 
    if ( strlen( options->getRegressionOutputFilename() ) ){
      outputstream.open( options->getRegressionOutputFilename(), ios::out );
      if( !outputstream )
	{
	  Log.setDisplayMode(On);
	  Log << "ERROR: Couldn't open regparamfile\n";
	  exit( 1 );
	}
      else{
	Log.setDisplayMode(Quiet);
	Log << "Writing regression parameters to " << options->getRegressionOutputFilename() << "\n";
	InitializeOutputFile(options, individuals, PopulationLabels);
      }
    }
    else{
      Log << "No regparamfile given\n";
    }
  }
}

void Regression::InitializeOutputFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
				      const std::string* const PopulationLabels)
{
  // Header line of paramfile
  for( int kk = 0; kk < options->getNumberOfOutcomes(); kk++ ){
      outputstream << "intercept\t";
      const Vector_s& labels = individuals->getCovariateLabels();
      for( unsigned i = 0; i < labels.size(); i++ ){
	outputstream << labels[i] << "\t";
      }
      if( !options->getTestForAdmixtureAssociation() && !options->getHapMixModelIndicator() )
	for( int k = 1; k < options->getPopulations(); k++ ){
	  outputstream << "slope." << PopulationLabels[k] << "\t";
	}
      if( individuals->getOutcomeType(kk) == Continuous ){
	outputstream<< setprecision(6) << "precision\t";
      }
    }
  outputstream << endl;
}

void Regression::Initialise(unsigned Number, const IndividualCollection* const individuals){
  //set regression number for this object
  RegNumber = Number;
  
  // ** Objects common to both regression types
  NumCovariates = individuals->GetNumCovariates();
  beta = new double[ NumCovariates ];
  lambda = 1.0; 
}

void Regression::Initialise(unsigned Number, double priorPrecision, const IndividualCollection* const individuals, LogWriter &Log){
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
  
  //  std::vector<double> v = individuals->getOutcome(RegNumber);
  //double p = accumulate(v.begin(), v.end(), 0.0, std::plus<double>()) / (double)v.size();
  //     else if(RegType == Linear)
  //       betamean[0] = p;

  //initialise regression params at prior mean
  for(int j = 0; j < NumCovariates; ++j){
    beta[j] = betamean[j];
  }

  betaprecision = new double[NumCovariates];
  double outcomeSampleVariance = individuals->getSampleVarianceOfOutcome(RegNumber);
  betaprecision[0] = priorPrecision / outcomeSampleVariance;
  Log << "\nGaussian priors on " << (RegType==Linear? "Linear" : "Logistic") << " regression parameters with zero means and precisions\n ("<< betaprecision[0];
  
  for(int j = 1; j < NumCovariates; ++j){
    betaprecision[j] = priorPrecision * individuals->getSampleVarianceOfCovariate(j) / outcomeSampleVariance;
    Log << ", " << betaprecision[j];
  }
  Log << ")\n";
  
  X = individuals->getCovariates();
  XtY = new double[NumCovariates];
  
  lambda = 1.0; 
  SumLambda = lambda;
  
}

void Regression::SetExpectedY(IndividualCollection *IC)const{
  IC->SetExpectedY(RegNumber, beta);
}

void Regression::Output(int iteration, const AdmixOptions *options, LogWriter &Log){
  //output to logfile
  if( iteration == -1 )
    {
      if( RegType != None )
	{
	  Log.setDisplayMode(Off);
	  Log.setPrecision(6);
	  if(options->getNumberOfOutcomes()==2)Log <<"\nRegression " <<(int)RegNumber << ": ";
          for( int j = 0; j < NumCovariates; j++ )
	    {
	      Log << beta[j] << "\t";
	    }
          if( RegType == Linear )
	    {
	      Log << lambda;
	    }
	}
    }
  //output to screen
  if( options->getDisplayLevel()>2 )
    {
      if(options->getNumberOfOutcomes()==2)cout << "\nRegression " << RegNumber << "\t";
      OutputParams(&cout);
      cout << endl;
    }
  //Output to paramfile after BurnIn
  if( iteration > options->getBurnIn() ){
    OutputParams(&outputstream);
    if(options->getNumberOfOutcomes()< 2 || RegNumber==1) outputstream << endl;
    //output new line in paramfile when last regression model
  }
}

void Regression::OutputErgodicAvg(int samples, std::ofstream *avgstream)const{
 //output to ergodicaveragefile
  if( RegType != None ){
    for( int j = 0; j < NumCovariates; j++ ){
      avgstream->width(9);
      *avgstream << setprecision(6) << SumBeta[j] / samples << "\t";
    }
    avgstream->width(9);
    if( RegType == Linear )
      *avgstream << setprecision(6) << SumLambda / samples << "\t";
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
  if(beta)//in case beta not allocated (happens if no regression model); may be unnecessary if beta initialised to 0
    return beta;
  else return NULL;
}
double Regression::getlambda()const{
  return lambda;
}
int Regression::getNumCovariates()const{
  return NumCovariates;
}



