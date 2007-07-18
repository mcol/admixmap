// *-*-C++-*-*
/** 
 *   LogisticRegression.h
 *   Class to represent and update parameters of a logistic regression model
 *   Copyright (c) 2006-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef LOGISTICREGRESSION_H
#define LOGISTICREGRESSION_H 1

#include "bclib/Regression.h"
#include "bclib/GaussianProposalMH.h"

BEGIN_BCLIB_NAMESPACE

///Struct to hold arguments for sampling logistic regression parameters
typedef struct{
  ///number of individuals
  int n;
  ///numberof covariates
  int d;
  ///index of current parameter
  int index;
  ///prior mean
  double beta0;
  ///prior precision 
  double priorprecision;
  double XtY;
  const double* Covariates;
  ///regression parameters
  const double* beta;
  ///for annealing
  double coolness;

}BetaArgs;

///class to sample the parameters of a logistic regression
class LogisticRegression : public bclib::Regression{
public:
  LogisticRegression(unsigned Number, double priorPrecision, const bclib::DataMatrix& Covars, const bclib::DataMatrix& Outcome, 
		     LogWriter &Log);
  ~LogisticRegression();

  double getDispersion()const;
  double DerivativeInverseLinkFunction(unsigned i)const;
  void Update(bool sumbeta, const std::vector<double>& Outcome, double coolness);
  double getLogLikelihood(const std::vector<double>& Outcome)const;
  double getLogLikelihoodAtPosteriorMeans(int iterations, const std::vector<double>& Outcome);
private:
  // ** Logistic Regression Objects
  GaussianProposalMH* BetaSampler;
  BetaArgs BetaParameters;
  int acceptbeta;

  void SetExpectedY(const double* const beta) ;

  static double lr( const double beta, const void* const vargs );
  
  static double dlr( const double beta, const void* const vargs );
  
  static double ddlr( const double beta, const void* const vargs );

  LogisticRegression();
  void Initialise(double priorPrecision, const bclib::DataMatrix& Covars, const bclib::DataMatrix& Outcome, 
		  LogWriter &Log);

};
END_BCLIB_NAMESPACE
#endif
