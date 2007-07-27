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
#include "bclib/LinearRegression.h"
#include "bclib/LogisticRegression.h"
#include "bclib/CoxRegression.h"
#include "AlleleFreqs.h"
#include "Options.h"
#include "InputData.h"
#include "IndividualCollection.h"
#include "Chromosome.h"
#include "ResidualLDTest.h"
#include "HWTest.h"
#include "Annealer.h"

//prototypes for miscellaneous functions
void CreateDirectory(const char* dirname, bool DeleteExistingFiles=true);

void WriteIterationNumber(const int iteration, const int width, int displayLevel);

void PrintCopyrightNotice(bclib::LogWriter & Log);

void PrintUsage(const char* ProgName);

void ThrowException(const string& msg, bclib::LogWriter & Log);

void PrintBuildInfo(bclib::LogWriter& Log);

///Abstract base class for top-level Model classes
///Currently, there are two derived classes: AdmixMapModel (admixture mapping) and HapMixModel (haplotype mixture model)
class Model{
public:
  
  Model();
  virtual ~Model();

  void Run(Options& options, InputData& data, bclib::LogWriter& Log);

  virtual void Iterate(const int & samples, const int & burnin, const double* Coolnesses, unsigned coolness,
		       Options & options, InputData & data, bclib::LogWriter& Log, 
		       double & SumEnergy, double & SumEnergySq, 
		       bool AnnealedRun) = 0;
  
  void InitialiseRegressionObjects(Options & options, InputData& data,  bclib::LogWriter& Log) ;
  virtual void InitialiseTests(Options& options, const InputData& data, bclib::LogWriter& Log) = 0;
  virtual void InitializeErgodicAvgFile(const Options& options, bclib::LogWriter &Log,  
					const Vector_s& PopLabels, const Vector_s& CovariateLabels) = 0;

  virtual void ResetStepSizeApproximators(int resetk);
  virtual void PrintAcceptanceRates(const Options& options, bclib::LogWriter& Log) = 0;
  
  std::vector<bclib::Regression*>& getRegression(){return R;};
  virtual unsigned getNumIndividuals()const = 0; 
  virtual double* getSumEnergy()const = 0;
  virtual double* getSumEnergySq()const = 0; 
  virtual double getDevianceAtPosteriorMean(const Options& options, bclib::LogWriter& Log) = 0;
  virtual void Finalize(const Options& options, bclib::LogWriter& Log, const InputData& data)=0 ;
  
  //this function is used only in admixmap model
  virtual void getOnePopOneIndLogLikelihood(bclib::LogWriter& , const std::vector<std::string>& ){};
  
protected:
  void InitialiseGenome(Genome& G, const Options& options, InputData& data, bclib::LogWriter& Log);
  //virtual void OutputParameters(int iteration, const Options *options, bclib::LogWriter& Log) = 0;

  void OutputErgodicAvgDeviance(int samples, double & SumEnergy, double & SumEnergySq);

  void AccumulateEnergy(const double* Coolnesses, unsigned coolness, const Options & options, double & SumEnergy, double & SumEnergySq, 
			double& AISz, bool AnnealedRun, int iteration);
  AlleleFreqs* pA;
  IndividualCollection *IC;
  vector<bclib::Regression*> R;//vector of regression pointers

  std::ofstream avgstream; //output to ErgodicAverageFile
  bclib::RObjectWriter paramstream;

  HWTest HWtest;
  ResidualLDTest ResidualAllelicAssocScoreTest;
  double AISsumlogz; //for computing marginal likelihood by Annealed Importance Sampling
  std::ofstream loglikelihoodfile;
  Annealer _Annealer;

  void Start(Options& options, InputData& data, bclib::LogWriter& Log);
  void Finish(Options& options, InputData& data, bclib::LogWriter& Log);
};

#endif /* MODEL_TOP_H */
