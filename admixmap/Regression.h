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
#ifndef REGRESSIONBASE_H
#define REGRESSIONBASE_H 1

#include "common.h"
#include <fstream>
#include "LogWriter.h"

class IndividualCollection;//to avoid circular includes

///Abstract Base Class for a generic Regression
class Regression{
public:
  Regression();
  virtual ~Regression();
  virtual void Initialise(unsigned RegNumber, double priorPrecision, const IndividualCollection* const, LogWriter &) = 0;
  void Initialise(unsigned RegNumber, unsigned numCovariates);
  virtual void Update(bool sumbeta, const std::vector<double>& Outcome, const double* const Covariates, double coolness
#ifdef PARALLEL
			, MPI::Intracomm &Comm
#endif
	      ) = 0;
  virtual double getLogLikelihood(const std::vector<double>& Outcome)const = 0;
  virtual double getLogLikelihoodAtPosteriorMeans(int iterations, const std::vector<double>& Outcome) = 0;
  static void OpenOutputFile(const unsigned NumOutcomes, const char* const filename, LogWriter &Log);  
  virtual void InitializeOutputFile(const std::vector<std::string>& CovariateLabels, unsigned NumOutcomes);
  virtual void Output(const unsigned NumberOfOutcomes, bool toScreen, bool afterBurnIn);
  virtual void OutputParams(std::ostream* out)const;
  virtual void OutputErgodicAvg(int iteration, std::ofstream *avgstream)const;
  const double* getbeta() const;
  double getlambda() const ;
  int getNumCovariates()const;
  virtual double getDispersion()const = 0;
  virtual double DerivativeInverseLinkFunction(unsigned i)const = 0;
  const double* getExpectedOutcome()const;

protected:
  int NumCovariates, NumOutcomeVars, NumIndividuals;
  RegressionType RegType;
  unsigned RegNumber;

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

  static std::ofstream outputstream;
  void Initialise(unsigned Number, unsigned nCovariates, unsigned nIndivs, const double* const Covars);
  void SumParameters();
  static void getExpectedOutcome(const double* const beta, const double* const X, double* EY, int n, int d);
  static void getExpectedOutcome(const double* const beta, const double* const X, double* Y, int n, int dim, int index, double betaj);
};


#endif
