// *-*-C++-*-*
/* 
 *   HAPMIXMAP
 *   HapMixModel.h
 *   
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef HAPMIXMODEL_H
#define HAPMIXMODEL_H 1

#include "Model.h"
#include "PopHapMix.h"
#include "HapMixFreqs.h"
#include "HapMixIndividualCollection.h"
#include "MantelHaenszelTest.h"

class HapMixModel : public Model{
public:
    ~HapMixModel();
  void Initialise(AdmixOptions & options, InputData& data,  LogWriter& Log);
  void InitialiseTests(AdmixOptions& options, const InputData& data, LogWriter& Log);
  void SubIterate(int iteration, const int & burnin, AdmixOptions & options, InputData & data, 
		  LogWriter& Log, double & SumEnergy, double & SumEnergySq, 
		  bool AnnealedRun);

  void PrintAcceptanceRates(const AdmixOptions& options,LogWriter& Log);
  void Finalize(const AdmixOptions& options, LogWriter& Log, const InputData& data) ;
  double getDevianceAtPosteriorMean(const AdmixOptions* const options, LogWriter& Log);
  unsigned getNumIndividuals()const{return IC->getSize();}; 
  double* getSumEnergy()const{return 0;};
  double* getSumEnergySq()const{return 0;};
private:
  PopHapMix* L;
  HapMixFreqs A;
  MantelHaenszelTest MHTest;
  void UpdateParameters(int iteration, const AdmixOptions *options, 
			LogWriter& Log, const Vector_s& PopulationLabels, double coolness, bool anneal);
  void OutputParameters(int iteration, const AdmixOptions *options, LogWriter& Log);
  void InitializeErgodicAvgFile(const AdmixOptions* const options, LogWriter &Log,  
				const Vector_s& PopLabels, const Vector_s& CovariateLabels);
};


#endif /* HAPMIXMODEL_H */
