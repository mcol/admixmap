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
#include "InputHapMixData.h"
#include "HapMixGenome.hh"
#include "HapMixOptions.h"
#include "PopHapMix.h"
#include "HapMixFreqs.h"
#include "HapMixIndividualCollection.h"
#include "MantelHaenszelTest.h"
#include "HapMixAllelicAssocTest.h"

class HapMixModel : public Model{
public:
  HapMixModel();
  ~HapMixModel();
  void Initialise(HapMixOptions & options, InputHapMixData& data,  bclib::LogWriter& Log);
  void Iterate(const int & samples, const int & burnin, const double* Coolnesses, unsigned coolness,
	       Options & options, InputData & data, bclib::LogWriter& Log, 
	       double & SumEnergy, double & SumEnergySq, 
	       bool AnnealedRun);

  void PrintAcceptanceRates(const Options& options,bclib::LogWriter& Log);
  void Finalize(const Options& options, bclib::LogWriter& Log, const InputData& data) ;
  double getDevianceAtPosteriorMean(const Options& options, bclib::LogWriter& Log);
  unsigned getNumIndividuals()const{return IC->getSize();}; 
  double* getSumEnergy()const{return 0;};
  double* getSumEnergySq()const{return 0;};

private:
  HapMixGenome Loci;
  PopHapMix* L;
  HapMixFreqs A;
  MantelHaenszelTest MHTest;
  HapMixAllelicAssocTest AllelicAssocTest;
  HapMixIndividualCollection* HMIC;

  void ReadInitialValuesFromFile(unsigned startnum, const HapMixOptions& options, bclib::LogWriter& Log);
  void UpdateParameters(int iteration, const Options& _options, bclib::LogWriter&, 
			const InputData & data, const double* Coolnesses, unsigned coolness_index, bool anneal, 
			double & SumEnergy, double & SumEnergySq, double& AISz);
  void OutputParameters(int iteration, const Options& options, bclib::LogWriter& Log);
  void OutputErgodicAverages(int samples, double & SumEnergy, double & SumEnergySq);
  void OutputTests(HapMixOptions& options, InputData & data, bclib::LogWriter& Log  );
  void InitialiseTests(Options& options, const InputData& data, bclib::LogWriter& Log);
  void InitializeErgodicAvgFile(const Options& options, bclib::LogWriter &Log,  
				const Vector_s& PopLabels, const Vector_s& CovariateLabels);
  void WriteParamsAsRObjectDimensions(const HapMixOptions& options, const InputData& data);
};


#endif /* HAPMIXMODEL_H */
