// *-*-C++-*-*
/** 
 *   Regression.h
 *   Abstract base class for a Bayesian regression
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef REGRESSIONBASE_H
#define REGRESSIONBASE_H 1

#include "bclib/bclib.h"
#include "bclib/common.h"
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "bclib/LogWriter.h"
#include "bclib/DataMatrix.h"
#include "bclib/RObjectWriter.h"

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */



///Abstract Base Class for a generic Bayesian Regression
class Regression{
public:
  Regression(unsigned Number, RegressionType RT);
  virtual ~Regression();
  virtual void Initialise(double priorPrecision, 
			  const bclib::DataMatrix& Covars, const bclib::DataMatrix& Outcome, LogWriter &) = 0;
  //void Initialise(unsigned numCovariates);
  virtual void Update(bool sumbeta, const std::vector<double>& Outcome, double coolness) = 0;
  virtual double getLogLikelihood(const std::vector<double>& Outcome)const = 0;
  virtual double getLogLikelihoodAtPosteriorMeans(int iterations, const std::vector<double>& Outcome) = 0;
  static void OpenOutputFile(const unsigned NumOutcomes, const char* const filename, LogWriter &Log); 
  static void OpenExpectedYFile(const char* Filename, LogWriter & Log); 
  virtual void InitializeOutputFile(const std::vector<std::string>& CovariateLabels, unsigned NumOutcomes);
  virtual void OutputParams(unsigned NumOutcomes);
  virtual void OutputParams(bclib::Delimitedostream& out)const;
  virtual void OutputErgodicAvg(int iteration, std::ofstream& avgstream)const;
  virtual void OutputExpectedY();
  static void FinishWritingEYAsRObject(unsigned NumIterations, const std::vector<std::string>& Labels);
  const double* getbeta() const;
  double getlambda() const ;
  int getNumCovariates()const;
  virtual double getDispersion()const = 0;
  virtual double DerivativeInverseLinkFunction(unsigned i)const = 0;
  const double* getExpectedOutcome()const;
  double getExpectedOutcome(unsigned i)const;
  RegressionType getRegressionType()const;
protected:
  int NumCovariates;
  static int NumOutcomeVars;
  static int NumIndividuals;
  const RegressionType RegType;
  const unsigned RegNumber;

  double *beta;//regression parameters
  double *betamean; //beta prior mean
  double *betaprecision; //prior precision for beta
  double *SumBeta;//running sums (for ergodic averages)
  double lambda; //precision parameter
  double SumLambda;
  const double* Outcome;
  const double* Covariates;
  double* ExpectedY;
  double* XtY;

  static bclib::DelimitedFileWriter outputstream;///< stream for regression parameters
  static RObjectWriter EYStream;///< stream for expected outcomes

  void Initialise(unsigned nCovariates, unsigned nIndivs, const double* const Covars);
  void SumParameters();
  static void getExpectedOutcome(const double* const beta, const double* const X, double* EY, int n, int d);
  static void getExpectedOutcome(const double* const beta, const double* const X, double* Y, int n, int dim, int index, double betaj);

private:
  Regression();
};


/** @} */

END_BCLIB_NAMESPACE

#endif
