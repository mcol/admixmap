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
#include "PopAdmix.h"
#include "StratificationTest.h"
#include "DispersionTest.h"
#include "MisSpecAlleleFreqTest.h"
#include "AdmixIndividualCollection.h"

class AdmixMapModel : public Model{
public:
    ~AdmixMapModel();
  void Initialise(AdmixOptions & options, InputData& data,  LogWriter& Log);
  void InitialiseTests(AdmixOptions& options, const InputData& data, LogWriter& Log);
    void SubIterate(int iteration, const int & burnin,
		    AdmixOptions & options, InputData & data, LogWriter& Log, 
		    double & SumEnergy, double & SumEnergySq, 
		    bool AnnealedRun);

  void PrintAcceptanceRates(const AdmixOptions& options,LogWriter& Log);
  void Finalize(const AdmixOptions& options, LogWriter& Log, const InputData& data) ;
  void ResetStepSizeApproximators(int resetk);
  double getDevianceAtPosteriorMean(const AdmixOptions* const options, LogWriter& Log);
  void getOnePopOneIndLogLikelihood(LogWriter& Log, const std::vector<std::string>& PopLabels){IC->getOnePopOneIndLogLikelihood(Log, PopLabels);};
  unsigned getNumIndividuals()const{return IC->getSize();}; 
  double* getSumEnergy()const;
  double* getSumEnergySq()const;
private:
  PopAdmix* L;
  AlleleFreqs A;
  AdmixIndividualCollection* AdmixedIndividuals;
  StratificationTest StratTest;
  DispersionTest DispTest;
  MisSpecAlleleFreqTest AlleleFreqTest;

  void UpdateParameters(int iteration, const AdmixOptions *options, 
			LogWriter& Log, const Vector_s& PopulationLabels, double coolness, bool anneal);
  void OutputParameters(int iteration, const AdmixOptions *options, LogWriter& Log);
  void InitializeErgodicAvgFile(const AdmixOptions* const options, LogWriter &Log,  
				const Vector_s& PopLabels, const Vector_s& CovariateLabels);

};
#endif /* MODEL_TOP_H */
