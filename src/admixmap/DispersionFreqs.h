// *-*-C++-*-*
/** 
 *   DispersionFreqs.h
 *   header file for DispersionFreqs class
 *   Copyright (c) 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef DISPFREQS_H
#define DISPFREQS_H 1

#define ETASAMPLER 1 //1 = RANDOM WALK METROPOLIS
                     //2 = HAMILTONIAN

#include "AlleleFreqs.h"
#include "bclib/MuSampler.h"
#include "bclib/DispersionSampler.h"
#include "bclib/StepSizeTuner.h"
#include "bclib/AdaptiveRejection.h"
class AdmixOptions;
class InputAdmixData;

/// Class to hold allele/haplotype frequencies and their priors in a dispersion model.
class DispersionFreqs : public AlleleFreqs{

public:
  DispersionFreqs();
  ~DispersionFreqs();
  void Initialise(AdmixOptions* const options, InputAdmixData* const Data, Genome *pLoci, bclib::LogWriter &Log, bool MAP=false);
  void Update(IndividualCollection*IC , bool afterBurnIn, double coolness);
  void PrintPrior(const Vector_s&, bclib::LogWriter& Log)const;
  ///initialize output file for samples of dispersion parameters
  void InitializeEtaOutputFile(const AdmixOptions* const options, const Vector_s& PopulationLabels, bclib::LogWriter &Log);

  ///outputs ergodic averages of dispersion parameters (SumEta)  to ErgodicAverageFile
  virtual void OutputErgodicAvg( int iteration, std::ofstream *avgstream)const;
  ///output samples of dispersion parameters (eta) to dispparamfile
  void OutputEta(int iteration, const AdmixOptions *options, bclib::LogWriter &Log);

  void OutputFST();

  void UpdateFst();
  const double *GetStatsForEta( int , int locus)const;

  float getAlphaSamplerAcceptanceRate(int)const;
  float getAlphaSamplerStepsize(int)const;
  float getEtaSamplerAcceptanceRate(int)const;
  float getEtaSamplerStepsize(int)const;

private:
  bool IsHistoricAlleleFreq;//indicator for dispersion model
  bool CorrelatedAlleleFreqs;
  double *eta; //dispersion parameter
  double *SumEta;
  double *psi,*tau;// eta has Gamma prior with shape and scale parameters psi and tau
  double psi0;
  bclib::MuSampler *muSampler;
  //new sampler for eta
  bclib::DispersionSampler *EtaSampler;
 
  double **HistoricAlleleFreqs;
  double **HistoricAlleleCounts;
  bool calculateFST;
  double** Fst;
  double** SumFst;

  ///adaptive RW sampler for eta
  bclib::StepSizeTuner *TuneEtaSampler;
  int NumberOfEtaUpdates;
  int *NumberAccepted; 
  double *etastep;
  double etastep0;
  int w;//The eta sampler is tuned every w updates.

  // sampler for Dirichlet proportion parameters 
  std::vector<bclib::StepSizeTuner> *MuProposal;
  
  std::ofstream outputstream;//outputs eta to paramfile
  std::ofstream fstoutputstream;

  void OpenFSTFile(const std::string& ResultsDir, bclib::LogWriter &Log); 

  void LoadAlleleFreqs(const Matrix_s& NewFreqs, int i, unsigned row0, bool);
  void LoadAlleleFreqs(AdmixOptions* const options, InputAdmixData* const data, bclib::LogWriter &Log);
  void SetDefaultPriorParams(int i, double defaultpriorparams);
  void SampleAlleleFreqs(int, const double coolness);
  void SampleDirichletParams1D( int );
  void SampleDirichletParamsMultiDim( int);
  void SampleDirichletParams();
  void UpdatePriorAlleleFreqs( int, const std::vector<std::vector<double> >& );
  void SampleEtaWithRandomWalk(int k, bool updateSumEta);

  static double muEnergyFunction(unsigned K, const double * const alpha, const double* const *args);
  static void muGradient(unsigned K, const double * const alpha, const double* const *args, double *g);

};
// functions required to update proportion vector Mu with adaptive rejection sampler
// likelihood, 1st and 2nd derivatives of log-likelihood
//Note that these are not part of AlleleFreqs class
double fMu( double alpha, const void* const args );
double dfMu( double alpha, const void* const args );
double ddfMu( double alpha, const void* const args );

#endif
