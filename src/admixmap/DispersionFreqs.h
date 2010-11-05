//=============================================================================
//
// Copyright (C) 2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file DispersionFreqs.h
/// Definition of the DispersionFreqs class.
//=============================================================================

#ifndef DISPFREQS_H
#define DISPFREQS_H 1

#define ETASAMPLER 1 //1 = RANDOM WALK METROPOLIS
                     //2 = HAMILTONIAN

#include "AdmixFreqs.h"
#include "bclib/MuSampler.h"
#include "bclib/DispersionSampler.h"


/** \addtogroup admixmap
 * @{ */


///class to implement a dispersion model for allele frequencies
class DispersionFreqs : public AdmixFreqs{
public:
  DispersionFreqs();
  ~DispersionFreqs();

  void Initialise(Options* const options, InputData* const Data, Genome *pLoci, bclib::LogWriter &Log, bool MAP=false);
  void Update(IndividualCollection*IC , bool afterBurnIn, double coolness);
  void PrintPrior(const Vector_s&, bclib::LogWriter& Log)const;

  ///initialize output file for samples of dispersion parameters
  void InitializeEtaOutputFile(const AdmixOptions* const options, const Vector_s& PopulationLabels, bclib::LogWriter &Log);

  ///outputs ergodic averages of dispersion parameters (SumEta)  to ErgodicAverageFile
  virtual void OutputErgodicAvg( int iteration, std::ofstream *avgstream)const;
  ///write sampled values of dispersion parameters to file
  void OutputParams();
  ///write sampled values of dispersion parameters to stream
  void OutputParams(bclib::Delimitedostream& os)const;

  void OutputFST();
  void UpdateFst();

  void PrintAcceptanceRates(bclib::LogWriter& Log);

private:
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
  
  bclib::DelimitedFileWriter outputstream;//outputs eta to paramfile
  std::ofstream fstoutputstream;

  void OpenFSTFile(const std::string& ResultsDir, bclib::LogWriter &Log); 

  void LoadAlleleFreqs(Options* const options, InputData* const data, bclib::LogWriter &Log);
  void LoadAlleleFreqs(const Matrix_s& New, int i, unsigned row0, bool );
  void SampleAlleleFreqs(int i, double coolness);
  void SampleDirichletParams1D( int );
  void SampleDirichletParamsMultiDim( int);
  void SampleDirichletParams();
  void UpdatePriorAlleleFreqs( int, const std::vector<std::vector<double> >& );
  void SampleEtaWithRandomWalk(int k, bool updateSumEta);
  const double *GetStatsForEta( int , int locus)const;

  static double muEnergyFunction(unsigned K, const double * const alpha, const double* const *args);
  static void muGradient(unsigned K, const double * const alpha, const double* const *args, double *g);

  float getAlphaSamplerAcceptanceRate(int)const;
  float getAlphaSamplerStepsize(int)const;
  float getEtaSamplerAcceptanceRate(int)const;
  float getEtaSamplerStepsize(int)const;

};

// functions required to update proportion vector Mu with adaptive rejection sampler
// likelihood, 1st and 2nd derivatives of log-likelihood
//Note that these are not part of AlleleFreqs class
double fMu( double alpha, const void* const args );
double dfMu( double alpha, const void* const args );
double ddfMu( double alpha, const void* const args );


/** @} */


#endif
