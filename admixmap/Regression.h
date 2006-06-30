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
#include "AdmixOptions.h"
//#include "IndividualCollection.h"
#include "LogWriter.h"

class IndividualCollection;//to avoid circular includes

///Abstract Base Class for a generic Regression
class Regression{
public:
  Regression();
  virtual ~Regression();
  virtual void Initialise(unsigned RegNumber, double priorPrecision, const IndividualCollection* const, LogWriter &) = 0;
  void Initialise(unsigned Number, const IndividualCollection* const individuals);
  void SetExpectedY(IndividualCollection* IC)const;
  virtual void Update(bool sumbeta, IndividualCollection* individuals, double coolness
#ifdef PARALLEL
			, MPI::Intracomm &Comm
#endif
	      ) = 0;
  virtual double getLogLikelihood(const IndividualCollection* const IC)const = 0;
  virtual double getLogLikelihoodAtPosteriorMeans(IndividualCollection *IC, int iterations) = 0;
  static void OpenOutputFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
  			     const std::string *PopulationLabels, LogWriter &Log);  
  static void InitializeOutputFile(const AdmixOptions* const , const IndividualCollection* const individuals, 
  				   const std::string* const PopulationLabels);
  void Output(int iteration, const AdmixOptions *, LogWriter &Log);
  virtual void OutputParams(ostream* out) = 0;
  void OutputErgodicAvg(int iteration, std::ofstream *avgstream)const;
  const double* getbeta() const;
  double getlambda() const ;
  int getNumCovariates()const;
  virtual double getDispersion()const = 0;

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
  const double* Y;
  double* XtY;
  const double *X;

  static std::ofstream outputstream;
  void SumParameters();
};


#endif
