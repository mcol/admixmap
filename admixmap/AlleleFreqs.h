// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   AlleleFreqs.h 
 *   header file for AlleleFreqs class
 *   Copyright (c) 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#ifndef ALLELEFREQS_H
#define ALLELEFREQS_H 1
#include "InputData.h"
#include "Genome.h"
#include "AdmixOptions.h"
#include "LogWriter.h"

class AlleleFreqs{

public:
  AlleleFreqs(Genome *pLoci);
  ~AlleleFreqs();
  void Initialise(AdmixOptions *options, const Matrix_d& etaprior,LogWriter *Log, std::string *);
  void Update(int iteration,int);

  //initialize output file for samples of dispersion parameters
  void InitializeEtaOutputFile(AdmixOptions *options, std::string *PopulationLabels, LogWriter *Log);

  //outputs ergodic averages of dispersion parameters (SumEta)  to ErgodicAverageFile
  void OutputErgodicAvg( int iteration,AdmixOptions *options, std::ofstream *avgstream);
  //output samples of dispersion parameters (eta) to dispparamfile
  void OutputEta(int iteration, AdmixOptions *options, std::ofstream *LogFileStreamPtr);

  Genome *getLoci();
  CompositeLocus *getLocus(int);

  int GetNumberOfCompositeLoci();

  void OutputAlleleFreqs();
  void CloseOutputFile(int iterations, string* PopulationLabels);

  void OutputFST(bool IsPedFile);

  void LoadAlleleFreqs(AdmixOptions *options, Chromosome ***chrm,LogWriter *Log, InputData *data);

  void ResetAlleleCounts();//resets Allelecounts to zero at start of iteration
  bool IsRandom();
  void UpdateFst();
  double *GetStatsForEta( int , int locus);
  double GetAlleleProbsMAP( int x, int ancestry , int locus);
  Vector_d GetPriorAlleleFreqs( int locus, int population );
  Vector_i GetAlleleCounts( int locus, int population );
  Vector_d getAlleleFreqsMAP( int locus, int population );
  Matrix_d &GetAlleleFreqs(int locus);
  Matrix_d *GetAlleleFreqs();
  Matrix_i &GetAlleleCounts(int locus);
  
  Matrix_d AlleleFreqs::GetSumAlleleFreqs(int locus);//is this used?

  void UpdateAlleleCounts(int locus, int h[2], Vector_i ancestry );
  void UpdateAlleleCounts_HaploidData(int locus, unsigned short **genotype, int ancestry );
  void ResetSumAlleleFreqs();
  void setAlleleFreqsMAP();
 
 // function to merge rare haplotypes for construction of score tests
  void SetMergedHaplotypes(Vector_d *alpha0, std::ofstream *LogFileStreamPtr, bool IsPedFile);

private:
  int Populations, NumberOfCompositeLoci;
  double *eta; //dispersion parameter
  double *SumEta;
  double *psi,*tau;// eta has Gamma prior with shape and scale parameters psi and tau
  double psi0;
 
  Matrix_d *Freqs;// allele frequencies except for last allele
  Matrix_d *AlleleFreqsMAP; // posterior mode of allele freqs
  Matrix_d *HistoricAlleleFreqs;
  Matrix_i *AlleleCounts;
  Matrix_d *HistoricLikelihoodAlleleFreqs;
  Matrix_d *PriorAlleleFreqs;

  Matrix_d *SumAlleleFreqs;// used to compute ergodic average

  dmatrix Fst;
  dmatrix SumFst;
  bool IsHistoricAlleleFreq;//indicator for dispersion model
  bool RandomAlleleFreqs;//indicator for whether allele freqs are fixed or random

  Genome *Loci;//pointer to Loci object

  TuneRW *TuneEtaSampler;
  int Number,w; // Number is the number of updates of eta. The eta sampler is tuned every w updates. 

  double *etastep;
  double etastep0;

  int *NumberAccepted;
  double *SumAcceptanceProb;

  double *pp;//used to set merged haplotypes, which are used in the allelic association test

//    DARS SampleMu;
   std::vector<TuneRW> *MuProposal;

  std::ofstream allelefreqoutput;// object to output allele frequencies
  std::ofstream outputstream;//outputs eta to paramfile
  std::ofstream fstoutputstream;

  void OpenFSTFile(AdmixOptions *options,LogWriter *Log); 

  // we have three different functions to initialize allele freqs
  // would be simpler to have one method
  // these functions are called by method LoadAlleleFreqs  
  void InitialiseAlleleFreqs(Matrix_d NewAlleleFreqs, int i, int Populations);
  void InitialisePriorAlleleFreqs(Matrix_d NewFreqs, int i, bool fixed, bool Historic);
  void SetDefaultAlleleFreqs(int Pops);

  void SamplePriorAlleleFreqs1D( int );
  void SamplePriorAlleleFreqsMultiDim( int);
  void SampleAlleleFreqs(int, int);
  void UpdatePriorAlleleFreqs( int, const std::vector<Vector_d>& );

  void OpenOutputFile(AdmixOptions *options);

};
// functions required to update proportion vector Mu with adaptive rejection sampler
// likelihood, 1st and 2nd derivatives of log-likelihood
//Note that these are not part of AlleleFreqs class
double fMu( Vector_d &, Matrix_i &, Matrix_d &, double );
double dfMu( Vector_d &, Matrix_i &, Matrix_d &, double );
double ddfMu( Vector_d &, Matrix_i &, Matrix_d &, double );





#endif
