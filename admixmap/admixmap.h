// *-*-C++-*-*
/* 
 *   ADMIXMAP
 *   admixmap.h 
 *   
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef ADMIXMAP_H
#define ADMIXMAP_H 1

#include <stdlib.h>    /* for exit, strtol */

#include <string.h>    /* for strcmp, strcpy */
#include <string>
#include "common.h"
#include "PopAdmix.h"
#include "PopHapMix.h"
#include "regression/LinearRegression.h"
#include "regression/LogisticRegression.h"
#include "regression/CoxRegression.h"
#include "AlleleFreqs.h"
#include "AdmixOptions.h"
#include "InputData.h"
#include "IndividualCollection.h"
#include "Chromosome.h"
#include "StratificationTest.h"
#include "DispersionTest.h"
#include "ScoreTests.h"
#include "MisSpecAlleleFreqTest.h"
#include "HWTest.h"
#include "Comms.h"

void MakeResultsDir(const char* dirname, bool verbose);

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
  void Iterate(const int & samples, const int & burnin, const double* Coolnesses, unsigned coolness,
		       AdmixOptions & options, InputData & data, const Genome & Loci, LogWriter& Log, 
		       double & SumEnergy, double & SumEnergySq, 
		       double& logz, bool AnnealedRun, ofstream & loglikelihoodfile);

  virtual void SubIterate(int iteration, const int& burnin, AdmixOptions & options, InputData & data, const Genome & Loci, 
			  LogWriter& Log, double & SumEnergy, double & SumEnergySq, 
			   bool AnnealedRun) = 0;

  virtual void Initialise(Genome& Loci, AdmixOptions & options, InputData& data,  LogWriter& Log) = 0;
  virtual void InitialiseTests(AdmixOptions& options, const InputData& data, const Genome& Loci, 
		       LogWriter& Log) = 0;
  virtual void InitializeErgodicAvgFile(const AdmixOptions* const options, LogWriter &Log,  
				const Vector_s& PopLabels, const Vector_s& CovariateLabels) = 0;
  virtual void ResetStepSizeApproximators(int resetk);
  virtual void PrintAcceptanceRates(const AdmixOptions& options, const Genome& Loci,LogWriter& Log) = 0;

  std::vector<Regression*>& getRegression(){return R;};
  unsigned getNumIndividuals()const{return IC->getSize();}; 
  double* getSumEnergy()const{return IC->getSumEnergy();};
  double* getSumEnergySq()const{return IC->getSumEnergySq();}; 
  virtual double getDevianceAtPosteriorMean(const AdmixOptions* const options, Genome* Loci, LogWriter& Log) = 0;
  void getOnePopOneIndLogLikelihood(LogWriter& Log, const std::vector<std::string>& PopLabels){IC->getOnePopOneIndLogLikelihood(Log, PopLabels);};
  virtual void Finalize(const AdmixOptions& options, LogWriter& Log, const InputData& data, const Genome& Loci)=0 ;

 protected:

  virtual void UpdateParameters(int iteration, const AdmixOptions *options, 
			const Genome *Loci, LogWriter& Log, const Vector_s& PopulationLabels, double coolness, bool anneal) = 0;
  virtual void OutputParameters(int iteration, const AdmixOptions *options, LogWriter& Log) = 0;

  void OutputErgodicAvgDeviance(int samples, double & SumEnergy, double & SumEnergySq);

  IndividualCollection *IC;
  AlleleFreqs A;
  vector<Regression*> R;//vector of regression pointers

  std::ofstream avgstream; //output to ErgodicAverageFile
  StratificationTest StratTest;
  DispersionTest DispTest;
  MisSpecAlleleFreqTest AlleleFreqTest;
  HWTest HWtest;
  ScoreTests Scoretests;
 private:

};

class AdmixMapModel : public Model{
public:
    ~AdmixMapModel();
  void Initialise(Genome& Loci, AdmixOptions & options, InputData& data,  LogWriter& Log);
  void InitialiseTests(AdmixOptions& options, const InputData& data, const Genome& Loci, 
		       LogWriter& Log);
    void SubIterate(int iteration, const int & burnin,
		    AdmixOptions & options, InputData & data, const Genome & Loci, LogWriter& Log, 
		    double & SumEnergy, double & SumEnergySq, 
		    bool AnnealedRun);

  void PrintAcceptanceRates(const AdmixOptions& options, const Genome& Loci,LogWriter& Log);
  void Finalize(const AdmixOptions& options, LogWriter& Log, const InputData& data, const Genome& Loci) ;
  void ResetStepSizeApproximators(int resetk);
    double getDevianceAtPosteriorMean(const AdmixOptions* const options, Genome* Loci, LogWriter& Log);
private:
  PopAdmix* L;
  void UpdateParameters(int iteration, const AdmixOptions *options, 
			const Genome *Loci, LogWriter& Log, const Vector_s& PopulationLabels, double coolness, bool anneal);
  void OutputParameters(int iteration, const AdmixOptions *options, LogWriter& Log);
  void InitializeErgodicAvgFile(const AdmixOptions* const options, LogWriter &Log,  
				const Vector_s& PopLabels, const Vector_s& CovariateLabels);

};
class HapMixModel : public Model{
public:
    ~HapMixModel();
  void Initialise(Genome& Loci, AdmixOptions & options, InputData& data,  LogWriter& Log);
  void InitialiseTests(AdmixOptions& options, const InputData& data, const Genome& Loci, 
		       LogWriter& Log);
  void SubIterate(int iteration, const int & burnin, AdmixOptions & options, InputData & data, const Genome & Loci, 
		  LogWriter& Log, double & SumEnergy, double & SumEnergySq, 
		  bool AnnealedRun);

  void PrintAcceptanceRates(const AdmixOptions& options, const Genome& Loci,LogWriter& Log);
  void Finalize(const AdmixOptions& options, LogWriter& Log, const InputData& data, const Genome& Loci) ;
  double getDevianceAtPosteriorMean(const AdmixOptions* const options, Genome* Loci, LogWriter& Log);
private:
  PopHapMix* L;
  void UpdateParameters(int iteration, const AdmixOptions *options, 
			const Genome *Loci, LogWriter& Log, const Vector_s& PopulationLabels, double coolness, bool anneal);
  void OutputParameters(int iteration, const AdmixOptions *options, LogWriter& Log);
  void InitializeErgodicAvgFile(const AdmixOptions* const options, LogWriter &Log,  
				const Vector_s& PopLabels, const Vector_s& CovariateLabels);
};


#endif /* ADMIXMAP_H */
