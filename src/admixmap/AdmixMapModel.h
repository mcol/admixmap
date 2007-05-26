// *-*-C++-*-*
/* 
 *   ADMIXMAP
 *   AdmixMapModel.h
 *   
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef ADMIXMAPMODEL_H
#define ADMIXMAPMODEL_H 1

#include "Model.h"
#include "Genome.h"
#include "AdmixOptions.h"
#include "InputAdmixdata.h"
#include "PopAdmix.h"
#include "StratificationTest.h"
#include "DispersionTest.h"
#include "MisSpecAlleleFreqTest.h"
#include "AdmixIndividualCollection.h"
#include "DispersionFreqs.h"
#include "ScoreTests.h"

class InputAdmixData;

class AdmixMapModel : public Model{
public:
    ~AdmixMapModel();
  void Initialise(AdmixOptions & options, InputAdmixData& data,  LogWriter& Log);
  void InitialiseTests(Options& options, const InputData& data, LogWriter& Log);
  void Iterate(const int & samples, const int & burnin, const double* Coolnesses, unsigned coolness,
	       Options & options, InputData & data, LogWriter& Log, 
	       double & SumEnergy, double & SumEnergySq, 
	       bool AnnealedRun);
  void SubIterate(int iteration, const int & burnin,
                  Options & options, InputData & data, LogWriter& Log, 
                  double & SumEnergy, double & SumEnergySq, 
                  bool AnnealedRun);

  void PrintAcceptanceRates(const Options& options,LogWriter& Log);
  void Finalize(const Options& options, LogWriter& Log, const InputData& data) ;
  void ResetStepSizeApproximators(int resetk);
  double getDevianceAtPosteriorMean(const Options* const options, LogWriter& Log);
  void getOnePopOneIndLogLikelihood(LogWriter& Log, const std::vector<std::string>& PopLabels){IC->getOnePopOneIndLogLikelihood(Log, PopLabels);};
  unsigned getNumIndividuals()const{return IC->getSize();}; 
  double* getSumEnergy()const;
  double* getSumEnergySq()const;
private:
  Genome Loci;
  PopAdmix* L;
  DispersionFreqs A;
  AdmixIndividualCollection* AdmixedIndividuals;
  StratificationTest StratTest;
  DispersionTest DispTest;
  MisSpecAlleleFreqTest AlleleFreqTest;
  ScoreTests Scoretests;

  void UpdateParameters(int iteration, const Options *options, 
			LogWriter& Log, const Vector_s& PopulationLabels, const double* Coolnesses, double coolness, bool anneal);
  void OutputParameters(int iteration, const AdmixOptions *options, LogWriter& Log);
  void InitializeErgodicAvgFile(const AdmixOptions* const options, LogWriter &Log,  
				const Vector_s& PopLabels, const Vector_s& CovariateLabels);

};
#endif /* MODEL_TOP_H */
