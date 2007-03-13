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

#include <cstdlib>    /* for exit, strtol */
#include <cstring>    /* for strcmp, strcpy */
#include <string>
#include "common.h"
#include "regression/LinearRegression.h"
#include "regression/LogisticRegression.h"
#include "regression/CoxRegression.h"
#include "AlleleFreqs.h"
#include "Options.h"
#include "InputData.h"
#include "IndividualCollection.h"
#include "Chromosome.h"
#include "ScoreTests.h"
#include "HWTest.h"
#include "Comms.h"
#include "Annealer.h"

void MakeResultsDir(const char* dirname, bool verbose, bool DeleteExistingFiles=true);

void WriteIterationNumber(const int iteration, const int width, int displayLevel);

void PrintCopyrightNotice(LogWriter & Log);

void PrintOptionsMessage();

void ThrowException(const string& msg, LogWriter & Log);

void PrintBuildInfo(LogWriter& Log);

///Abstract base class for top-level Model classes
///Currently, there are two derived classes: AdmixMapModel (admixture mapping) and HapMixModel (haplotype mixture model)
class Model{
public:
  
  Model();
  virtual ~Model();

  void Run(Options& options, InputData& data, LogWriter& Log, 
		 int NumAnnealedRuns);
  void TestIndivRun(Options& options, InputData& data, LogWriter& Log,  
		    int NumAnnealedRuns);

  void Iterate(const int & samples, const int & burnin, const double* Coolnesses, unsigned coolness,
	       Options & options, InputData & data, LogWriter& Log, 
	       double & SumEnergy, double & SumEnergySq, 
	       bool AnnealedRun);
  
  virtual void SubIterate(int iteration, const int& burnin, Options & options, InputData & data, 
                          LogWriter& Log, double & SumEnergy, double & SumEnergySq, 
                          bool AnnealedRun) = 0;
  
  //virtual void Initialise(Options & options, InputData& data,  LogWriter& Log) = 0;
  void InitialiseRegressionObjects(Options & options, InputData& data,  LogWriter& Log) ;
  virtual void InitialiseTests(Options& options, const InputData& data, LogWriter& Log) = 0;
  //virtual void InitializeErgodicAvgFile(const Options* const options, LogWriter &Log,  
  //				const Vector_s& PopLabels, const Vector_s& CovariateLabels) = 0;
  virtual void ResetStepSizeApproximators(int resetk);
  virtual void PrintAcceptanceRates(const Options& options, LogWriter& Log) = 0;
  
  std::vector<Regression*>& getRegression(){return R;};
  virtual unsigned getNumIndividuals()const = 0; 
  virtual double* getSumEnergy()const = 0;
  virtual double* getSumEnergySq()const = 0; 
  virtual double getDevianceAtPosteriorMean(const Options* const options, LogWriter& Log) = 0;
  virtual void Finalize(const Options& options, LogWriter& Log, const InputData& data)=0 ;
  
  //this function is used only in admixmap model
  virtual void getOnePopOneIndLogLikelihood(LogWriter& , const std::vector<std::string>& ){};
  
protected:
  void InitialiseLoci(const Options& options, InputData& data, LogWriter& Log);
  virtual void UpdateParameters(int iteration, const Options *options, 
                                LogWriter& Log, const Vector_s& PopulationLabels, const double* Coolnesses, double coolness, bool anneal) = 0;
  //virtual void OutputParameters(int iteration, const Options *options, LogWriter& Log) = 0;

  void OutputErgodicAvgDeviance(int samples, double & SumEnergy, double & SumEnergySq);

  AlleleFreqs* pA;
  Genome Loci;
  IndividualCollection *IC;
  vector<Regression*> R;//vector of regression pointers

  std::ofstream avgstream; //output to ErgodicAverageFile

  HWTest HWtest;
  ScoreTests Scoretests;
private:
  std::ofstream loglikelihoodfile;
  Annealer A;
  double AISsumlogz; //for computing marginal likelihood by Annealed Importance Sampling

  void Start(Options& options, InputData& data, LogWriter& Log, int NumAnnealedRuns);
  void Finish(Options& options, InputData& data, LogWriter& Log);
};

#endif /* MODEL_TOP_H */
