// *-*-C++-*-*
/** 
 *   DispersionFreqs.h
 *   Class to hold and update allele frequencies in a correlated frequencies model
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef CORRELATEDFREQS_H
#define CORRELATEDFREQS_H

#include "AdmixFreqs.h"
#include "bclib/DispersionSampler.h"

#define ETASAMPLER 1 //1 = RANDOM WALK METROPOLIS
                     //2 = HAMILTONIAN

///class to implement a correlated allele frequency model
class CorrelatedFreqs : public AdmixFreqs{
public:
  CorrelatedFreqs();
  void Initialise(Options* const a_options, InputData* const a_data, 
		  Genome *pLoci, bclib::LogWriter &Log, bool MAP);
  ///print parameters of prior on frequencies and their Dirichlet parameters
  void PrintPrior(const Vector_s& PopLabels, bclib::LogWriter& Log)const;
  void Update(IndividualCollection*IC , bool afterBurnIn, double coolness);
  ///write ergodic average of dispersion to file
  void OutputErgodicAvg( int samples, std::ofstream *avgstream)const;
  ///write dispersion parameter to file
  void OutputParams();
  ///write dispersion parameter to stream
  void OutputParams(bclib::Delimitedostream& os)const;
  void PrintAcceptanceRates(bclib::LogWriter& Log);

private:
  double eta; ///<dispersion parameter
  double SumEta;
  double psi,tau;// eta has Gamma prior with shape and scale parameters psi and tau
  double psi0;
  bclib::MuSampler *muSampler;
  //new sampler for eta
  bclib::DispersionSampler EtaSampler;
 
  ///adaptive RW sampler for eta
  bclib::StepSizeTuner TuneEtaSampler;
  int NumberOfEtaUpdates;
  int NumberAccepted; 
  double etastep;///< step size in random walk sampler
  double etastep0;///< initial step size
  int w;//The eta sampler is tuned every w updates.

  // sampler for Dirichlet proportion parameters 
  std::vector<bclib::StepSizeTuner> MuProposal;
  
  bclib::DelimitedFileWriter outputstream;//outputs eta to paramfile

  void LoadAlleleFreqs(const Matrix_s& NewFreqs, int i, unsigned row0, bool);
  void InitializeEtaOutputFile(const char* filename, const Vector_s& PopulationLabels, bclib::LogWriter &Log);
  void SetDefaultPriorParams(int i, double defaultpriorparams);
  void SampleDirichletParams();
  void SampleEtaWithRandomWalk(bool updateSumEta);
  const double* GetStatsForEta( int locus )const;
  void UpdatePriorAlleleFreqs(const vector<vector<double> >& mu);
  void SampleAlleleFreqs(int i, double coolness);
  float getEtaSamplerAcceptanceRate()const;
  float getEtaSamplerStepsize()const;
};

#endif
