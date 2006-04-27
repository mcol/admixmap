// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   Regression.h
 *   Class to represent and update parameters of a regression model as used in admixmap
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef REGRESSION_H
#define REGRESSION_H 1

#include "common.h"
#include "Gaussian.h"
#include "AdmixOptions.h"
#include "IndividualCollection.h"
#include "GaussianProposalMH.h"
#include "LogWriter.h"

typedef struct{
  int n;//number of individuals
  int d;//numberof covariates
  int index;//index of current parameter
  double beta0;//prior mean
  double lambda;//prior precision 
  double XtY;
  const double* Covariates;
  const double* beta;//regression parameters
  double coolness;//for annealing

}BetaArgs;

class Regression{

public:
  Regression();
   ~Regression();
  void Initialise(unsigned RegNumber, double priorPrecision, const IndividualCollection* const, LogWriter &);
  void SetExpectedY(IndividualCollection* IC)const;
  void Update(bool sumbeta, IndividualCollection* individuals, double coolness);
  static void OpenOutputFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
			     const std::string *PopulationLabels, LogWriter &Log);  
  static void InitializeOutputFile(const AdmixOptions* const , const IndividualCollection* const individuals, 
				   const std::string* const PopulationLabels);
  void Output(int iteration, const AdmixOptions *, LogWriter &Log);
  void OutputParams(ostream* out);
  void OutputErgodicAvg(int iteration, std::ofstream *avgstream)const;
  const double* getbeta() const;
  double getlambda() const ;
  int getNumCovariates()const;
  double getDispersion()const;
  double getLogLikelihood(const IndividualCollection* const IC)const;
  double getLogLikelihoodAtPosteriorMeans(IndividualCollection *IC, int iterations);

private:
  int NumCovariates, NumOutcomeVars, NumIndividuals;
  RegressionType RegType;
  unsigned RegNumber;

  double *beta;//regression parameters
  double *betamean; //beta prior mean
  double *SumBeta;//running sums (for ergodic averages)
  double lambda; //precision parameter
  double SumLambda;
  const double* Y;

  // ** Linear Regression Objects
  double lambda0; //parameters of
  double lambda1; //Gamma prior for lambda
  double *betaprecision; //prior precision for beta
  double *R, *QY, *QX, *V, *betahat;
  Gaussian DrawBeta;//sampler

  void QRSolve(int dim1, int dim2, const double* a, const double* b, double* x);
  void SamplePrecision(double* lambda, const double* Y, const double* X, int NumIndivs, int NumCovars, double coolness);
  void SampleLinearRegressionParams(double* beta, const double* Y, const double* X, int NumIndivs, int NumCovars);
  void SampleLinearRegressionParametersWithAnnealing(const double* Y, const double* X, double* beta, double *lambda, 
							       double coolness);

  // ** Logistic Regression Objects
  GaussianProposalMH** BetaDrawArray;
  BetaArgs BetaParameters;
  int acceptbeta;
  double* XtY;
  const double *X;
  int *dims;
 
  static std::ofstream outputstream;//output to regparamfile

  void SumParameters();

  static double lr( const double beta, const void* const vargs );
  
  static double dlr( const double beta, const void* const vargs );
  
  static double ddlr( const double beta, const void* const vargs );
  

};



























#endif
