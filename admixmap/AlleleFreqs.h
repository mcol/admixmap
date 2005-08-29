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

#define ETASAMPLER 1 //1 = RANDOM WALK METROPOLIS
                     //2 = HAMILTONIAN
#include "InputData.h"
#include "Genome.h"
#include "AdmixOptions.h"
#include "LogWriter.h"


#include "AdaptiveRandomWalkMH.h"
class AdaptiveRandomWalkMH;
#if ETASAMPLER ==2
#include "HamiltonianMonteCarlo.h"
#endif


class AlleleFreqs{

public:
  AlleleFreqs(Genome *pLoci);
  ~AlleleFreqs();
  void Initialise(AdmixOptions *options, InputData *Data, LogWriter *Log);
  void Update(int iteration,int);

  //initialize output file for samples of dispersion parameters
  void InitializeEtaOutputFile(AdmixOptions *options, std::string *PopulationLabels, LogWriter *Log);

  //outputs ergodic averages of dispersion parameters (SumEta)  to ErgodicAverageFile
  void OutputErgodicAvg( int iteration,AdmixOptions *options, std::ofstream *avgstream);
  //output samples of dispersion parameters (eta) to dispparamfile
  void OutputEta(int iteration, AdmixOptions *options, LogWriter *Log);

  Genome *getLoci();
  CompositeLocus *getLocus(int);

  int GetNumberOfCompositeLoci();

  void OutputAlleleFreqs();
  void CloseOutputFile(int iterations, string* PopulationLabels);

  void OutputFST();

  void LoadAlleleFreqs(AdmixOptions *options, InputData *data);

  void ResetAlleleCounts();//resets Allelecounts to zero at start of iteration
  bool IsRandom();
  void UpdateFst();
  double *GetStatsForEta( int , int locus);
  double GetAlleleProbsMAP( int x, int ancestry , int locus);
  Vector_d GetPriorAlleleFreqs( int locus, int population );
  Vector_i GetAlleleCounts( int locus, int population );
  Vector_d getAlleleFreqsMAP( int locus, int population );
  Vector_d GetAlleleFreqs( int locus, int population );
  double *GetAlleleFreqs(int locus);
  double **GetAlleleFreqs();
  int *GetAlleleCounts(int locus);
  
  //Matrix_d AlleleFreqs::GetSumAlleleFreqs(int locus);//is this used?

  void UpdateAlleleCounts(int locus, int h[2], int ancestry[2], bool diploid );
  void ResetSumAlleleFreqs();
  void setAlleleFreqsMAP();
 
private:
  int Populations, NumberOfCompositeLoci;
  double *eta; //dispersion parameter
  double *SumEta;
  double *psi,*tau;// eta has Gamma prior with shape and scale parameters psi and tau
  double psi0;
 
  double **Freqs;// allele frequencies except for last allele
  double **AlleleFreqsMAP; // posterior mode of allele freqs
  double **HistoricAlleleFreqs;
  int **AlleleCounts;
  double **HistoricAlleleCounts;
  Matrix_d *PriorAlleleFreqs;

  //double **SumAlleleFreqs;// used to compute ergodic average

  double** Fst;
  double** SumFst;
  bool IsHistoricAlleleFreq;//indicator for dispersion model
  bool RandomAlleleFreqs;//indicator for whether allele freqs are fixed or random

  Genome *Loci;//pointer to Loci object

#if ETASAMPLER == 1
  //adaptive RW sampler for eta
  AdaptiveRandomWalkMH *TuneEtaSampler;
  int Number,w; // Number is the number of updates of eta. The eta sampler is tuned every w updates. 
  double *etastep;
  double etastep0;
  int *NumberAccepted;
  double *SumAcceptanceProb;
#elif ETASAMPLER == 2
  //Hamiltonian MC sampler for eta
  double **EtaArgs;
  double *logeta;
  HamiltonianMonteCarlo EtaSampler;
  double initialEtaStepsize;
  float targetEtaAcceptRate;
#endif
  
  //sampler for univariate mu (beta proportion parameters)
  //    DARS SampleMu;
  // DARS SampleMu is now initialized each time the allele freqs at a locus in a population are sampled

  // sampler for Dirichlet proportion parameters 
  std::vector<AdaptiveRandomWalkMH> *MuProposal;
  
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
#if ETASAMPLER == 2
  static double etaEnergyFunction(unsigned dim, const double* const theta, const double* const*args);
  static void etaGradient(unsigned dim,const double* const theta, const double* const* args, double *g);
#endif
};
// functions required to update proportion vector Mu with adaptive rejection sampler
// likelihood, 1st and 2nd derivatives of log-likelihood
//Note that these are not part of AlleleFreqs class
double fMu( const double*, const int *, const double*, double );
double dfMu( const double*, const int *, const double*, double );
double ddfMu( const double*, const int *, const double*, double );

#endif
