// *-*-C++-*-*
/** 
 *   LinearRegression.h 
 *   Class to represent and update parameters of a logistic regression model
 *   Copyright (c) 2006-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef LINEARREGRESSION_H
#define LINEARREGRESSION_H 1

#include "bclib/Regression.h"
#include "bclib/Gaussian.h"

BEGIN_BCLIB_NAMESPACE

///Class to sample the parameters of a linear regression
class LinearRegression : public bclib::Regression{
public:
  LinearRegression(unsigned Number, double priorPrecision, const bclib::DataMatrix& Covars, const bclib::DataMatrix& Outcome, 
		   LogWriter &Log);
  ~LinearRegression();
  void InitializeOutputFile(const std::vector<std::string>& CovariateLabels, unsigned NumOutcomes);

  double getDispersion()const;
  double DerivativeInverseLinkFunction(unsigned i)const;
  void OutputParams(bclib::Delimitedostream& out)const;
  void OutputErgodicAvg(int samples, std::ofstream& avgstream)const;
  void Update(bool sumbeta, const std::vector<double>& Outcome, double coolness);
  double getLogLikelihood(const std::vector<double>& Outcome)const;
  double getLogLikelihoodAtPosteriorMeans(int iterations, const std::vector<double>& Outcome);
private:
  double lambda0; //parameters of
  double lambda1; //Gamma prior for lambda
  double *R, *QY, *QX, *V, *betahat;
  Gaussian DrawBeta;//sampler

  void SetExpectedY(const double* const _beta);
  void QRSolve(int dim1, int dim2, const double* a, const double* b, double* x);
  void SamplePrecision(double* lambda, const double* Y, const double* X, int NumIndivs, int NumCovars, double coolness);
  void SampleLinearRegressionParams(double* beta, const double* Y, const double* X, int NumIndivs, int NumCovars);
  void SampleLinearRegressionParametersWithAnnealing(const double* Y, const double* X, double* beta, double *lambda, 
							       double coolness);

  LinearRegression();
  void Initialise(double priorPrecision, const bclib::DataMatrix& Covars, const bclib::DataMatrix& Outcome, 
		  LogWriter &Log);

};
END_BCLIB_NAMESPACE

#endif
