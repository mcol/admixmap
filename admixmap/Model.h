// *-*-C++-*-*
/* 
 *   Model.h
 *   
 *   Copyright (c) 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef MODEL_TOP_H
#define MODEL_TOP_H 1

#include <stdlib.h>    /* for exit, strtol */

#include <string.h>    /* for strcmp, strcpy */
#include <string>
#include "common.h"
#include "regression/LinearRegression.h"
#include "regression/LogisticRegression.h"
#include "regression/CoxRegression.h"
#include "AlleleFreqs.h"
#include "AdmixOptions.h"
#include "InputData.h"
#include "IndividualCollection.h"
#include "Chromosome.h"
#include "ScoreTests.h"
#include "HWTest.h"
#include "Comms.h"

void MakeResultsDir(const char* dirname, bool verbose, bool DeleteExistingFiles=true);

void WriteIterationNumber(const int iteration, const int width, int displayLevel);

void PrintCopyrightNotice(LogWriter & Log);

void PrintOptionsMessage();
void ThrowException(const string& msg, LogWriter & Log);

///Abstract base class for model used in ADMIXMAP. 
///Currently, there are two derived classes: AdmixMapModel (admixture mapping) and HapMixModel (haplotype mixture model)
class Model{
public:
  
  Model();
  virtual ~Model();
  void Run(AdmixOptions& options, InputData& data, LogWriter& Log);
  void Iterate(const int & samples, const int & burnin, const double* Coolnesses, unsigned coolness,
	       AdmixOptions & options, InputData & data, LogWriter& Log, 
	       double & SumEnergy, double & SumEnergySq, 
	       double& logz, bool AnnealedRun, ofstream & loglikelihoodfile);
  
  virtual void SubIterate(int iteration, const int& burnin, AdmixOptions & options, InputData & data, 
			  LogWriter& Log, double & SumEnergy, double & SumEnergySq, 
			  bool AnnealedRun) = 0;
  
  virtual void Initialise(AdmixOptions & options, InputData& data,  LogWriter& Log) = 0;
  void InitialiseRegressionObjects(AdmixOptions & options, InputData& data,  LogWriter& Log) ;
  virtual void InitialiseTests(AdmixOptions& options, const InputData& data, LogWriter& Log) = 0;
  virtual void InitializeErgodicAvgFile(const AdmixOptions* const options, LogWriter &Log,  
					const Vector_s& PopLabels, const Vector_s& CovariateLabels) = 0;
  virtual void ResetStepSizeApproximators(int resetk);
  virtual void PrintAcceptanceRates(const AdmixOptions& options, LogWriter& Log) = 0;
  
  std::vector<Regression*>& getRegression(){return R;};
  virtual unsigned getNumIndividuals()const = 0; 
  virtual double* getSumEnergy()const = 0;
  virtual double* getSumEnergySq()const = 0; 
  virtual double getDevianceAtPosteriorMean(const AdmixOptions* const options, LogWriter& Log) = 0;
  virtual void Finalize(const AdmixOptions& options, LogWriter& Log, const InputData& data)=0 ;
  
  //this function is used only in admixmap model
  virtual void getOnePopOneIndLogLikelihood(LogWriter& , const std::vector<std::string>& ){};
  
protected:
  void InitialiseLoci(const Options& options, InputData& data, LogWriter& Log);
  virtual void UpdateParameters(int iteration, const AdmixOptions *options, 
				LogWriter& Log, const Vector_s& PopulationLabels, double coolness, bool anneal) = 0;
  virtual void OutputParameters(int iteration, const AdmixOptions *options, LogWriter& Log) = 0;

  void OutputErgodicAvgDeviance(int samples, double & SumEnergy, double & SumEnergySq);

  AlleleFreqs* pA;
  Genome Loci;
  IndividualCollection *IC;
  vector<Regression*> R;//vector of regression pointers

  std::ofstream avgstream; //output to ErgodicAverageFile

  HWTest HWtest;
  ScoreTests Scoretests;
private:

};

#endif /* MODEL_TOP_H */
