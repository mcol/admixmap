//=============================================================================
//
// Copyright (C) 2002-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file AdmixMapModel.h
/// Definition of the AdmixMapModel class.
//=============================================================================

#ifndef ADMIXMAPMODEL_H
#define ADMIXMAPMODEL_H 1

#include "Model.h"
#include "Genome.h"
#include "InputAdmixData.h"
#include "PopAdmix.h"
#include "StratificationTest.h"
#include "DispersionTest.h"
#include "MisSpecAlleleFreqTest.h"
#include "AdmixIndividualCollection.h"
#include "ScoreTests.h"

class AdmixFreqs;
class AdmixOptions;


/** \addtogroup admixmap
 * @{ */


class AdmixMapModel : public Model{
public:
  AdmixMapModel();
    ~AdmixMapModel();
  void Initialise(AdmixOptions & options, InputAdmixData& data,  bclib::LogWriter& Log);
  void InitialiseTests(Options& options, const InputData& data, bclib::LogWriter& Log);
  void TestIndivRun(Options& options, InputData& data, bclib::LogWriter& Log);
  void Iterate(const double* Coolnesses, unsigned coolness,
	       Options & options, InputData & data, bclib::LogWriter& Log, 
	       double & SumEnergy, double & SumEnergySq, 
	       bool AnnealedRun);
  void SubIterate(int iteration,
                  Options & options, InputData & data, bclib::LogWriter& Log, 
                  double & SumEnergy, double & SumEnergySq, 
                  bool AnnealedRun);

  void PrintAcceptanceRates(const Options& options, bclib::LogWriter& Log,
                            const Vector_s& PopulationLabels);
  void Finalize(const Options& options, bclib::LogWriter& Log, const InputData& data) ;
  void ResetStepSizeApproximators(int resetk);
  double getDevianceAtPosteriorMean(const Options& options, bclib::LogWriter& Log);
  void getOnePopOneIndLogLikelihood(bclib::LogWriter& Log, const std::vector<std::string>& PopLabels){IC->getOnePopOneIndLogLikelihood(Log, PopLabels);};
  unsigned getNumIndividuals()const{return IC->getSize();}; 
  double* getSumEnergy()const;
  double* getSumEnergySq()const;
private:
  Genome Loci;
  PopAdmix* L;
  AdmixFreqs* A;
  AdmixIndividualCollection* AdmixedIndividuals;
  StratificationTest StratTest;
  DispersionTest DispTest;
  MisSpecAlleleFreqTest AlleleFreqTest;
  ScoreTests Scoretests;

  void UpdateParameters(int iteration, const Options& options, 
			bclib::LogWriter& Log, const Vector_s& PopulationLabels, const double* Coolnesses, double coolness, bool anneal);
  void OutputParameters(int iteration, const AdmixOptions *options, bclib::LogWriter& Log);
  void InitializeErgodicAvgFile(const Options& options, bclib::LogWriter &Log,  
				const Vector_s& PopLabels, const Vector_s& CovariateLabels);
  void WriteParamsAsRObjectDimensions(const AdmixOptions& options, const InputData& data);
};


/** @} */


#endif /* MODEL_TOP_H */
