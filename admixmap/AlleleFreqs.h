// *-*-C++-*-*
/** 
 *   AlleleFreqs.h 
 *   header file for AlleleFreqs class
 *   Copyright (c) 2005, 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef ALLELEFREQS_H
#define ALLELEFREQS_H 1

#define ETASAMPLER 1 //1 = RANDOM WALK METROPOLIS
                     //2 = HAMILTONIAN
#define FREQ_CONJUGATE_SAMPLER 1
#define FREQ_HAMILTONIAN_SAMPLER 2

#include "InputData.h"
#include "Genome.h"
#include "AdmixOptions.h"
#include "utils/LogWriter.h"
#include "AlleleFreqSampler.h"
#include "samplers/MuSampler.h"
#include "samplers/DispersionSampler.h"
#include "samplers/StepSizeTuner.h"
#include "common.h"
#include "samplers/AdaptiveRejection.h"

#ifndef PARALLEL
#define ARRAY2D
#endif
 
/**
   struct to hold allelecounts in either a 1d (where Number of alleles is fixed) or 2d array. 
   usage, 1d array:
   array_of_allelecounts AlleleCounts;
   AlleleCounts.stride = K*2;
   AlleleCounts.array = new int[ NumberOfLoci*K*2];
   
   usage, 2darray:
   AlleleCounts.array = new int*[NumberOfLoci ];
   for(unsigned i = 0; i < NumberOfLoci; ++i) AlleleCounts.array[i] = new int[K*NumberOfStates[i]];
   
   AlleleCounts[i]; //accesses counts for ith locus
   
   AlleleCounts[i][k*2 +a]; //accesses count of ath allele in kth pop at ith locus
   */
#ifdef ARRAY2D
typedef struct{

  int **array;

  int* operator[](unsigned i){//for reading or writing element
    return array[i];
  };
  const int* operator[](unsigned i)const{//for read-only
    return array[i];
  };
  void dealloc(int L){
    if(array){
      for(int i = L-1; i >=0 ; --i)
	if(array[i]){
	  delete[] array[i];
	  array[i] = 0;
	}
      delete[]array;
      array = 0;
    }
  };
}array_of_allelecounts;

/**
   struct to hold allelefreqs in either a 1d (where Number of alleles is fixed) or 2d array. 
   See array_of_allelecounts for details.
*/
typedef struct{

  double **array;

  double* operator[](unsigned i){//for reading or writing element
    return array[i];
  };
  const double* operator[](unsigned i)const{//for read-only
    return array[i];
  };
  void dealloc(int L){
    if(array){
      for(int i = L-1; i >=0 ; --i)
	if(array[i]){
	  delete[] array[i];
	  array[i] = 0;
	}
      delete[] array;
      array = 0;
    }
  };
}array_of_allelefreqs;
#else
typedef struct{
  int* array;
  unsigned stride;

  int* operator[](unsigned i){
    return array + i*stride;
  };
  const int* operator[](unsigned i)const{
    return array + i*stride;
  };
  void dealloc(int ){
    if(array){
      delete[] array;
      array = 0;
    }
  };
}array_of_allelecounts;
/**
   struct to hold allelefreqs in either a 1d (where Number of alleles is fixed) or 2d array. 
   See array_of_allelecounts for details.
*/

typedef struct{
  double* array;
  unsigned stride;

  double* operator[](unsigned i){
    return array + i*stride;
  };
  const double* operator[](unsigned i)const{
    return array + i*stride;
  };
  void dealloc(int ){
    if(array){
      delete[] array;
      array = 0;
    }
  };
}array_of_allelefreqs;
#endif

/// Class to hold allele/haplotype frequencies and their priors.
class AlleleFreqs{

public:
  AlleleFreqs();
  virtual ~AlleleFreqs();
  virtual void Initialise(AdmixOptions* const options, InputData* const Data, Genome *pLoci, LogWriter &Log);
  void AllocateAlleleCountArrays(unsigned K);
  void PrintPrior(const Vector_s&, LogWriter& Log)const;
  virtual void Update(IndividualCollection*IC , bool afterBurnIn, double coolness);

  ///initialize output file for samples of dispersion parameters
  void InitializeEtaOutputFile(const AdmixOptions* const options, const Vector_s& PopulationLabels, LogWriter &Log);

  ///outputs ergodic averages of dispersion parameters (SumEta)  to ErgodicAverageFile
  virtual void OutputErgodicAvg( int iteration, std::ofstream *avgstream)const;
  ///output samples of dispersion parameters (eta) to dispparamfile
  void OutputEta(int iteration, const AdmixOptions *options, LogWriter &Log);

  void OutputAlleleFreqs();
  void OutputAlleleFreqs(const char* filename);
  void CloseOutputFile(int iterations, const Vector_s& PopulationLabels);

  void OutputFST();

  ///resets Allelecounts to zero at start of iteration
  void ResetAlleleCounts(unsigned K);
  bool IsRandom()const;
  void UpdateFst();
  const double *GetStatsForEta( int , int locus)const;
  double GetAlleleProbsMAP( int x, int ancestry , int locus)const;
  std::vector<double> GetPriorAlleleFreqs( int locus, int population )const;
  std::vector<int> GetAlleleCounts( int locus, int population )const;
  std::vector<double> getAlleleFreqsMAP( int locus, int population )const;
  std::vector<double> GetAlleleFreqs( int locus, int population )const;
  const double *GetAlleleFreqs(int locus)const;
  const array_of_allelefreqs& GetAlleleFreqs()const;
  const int *GetAlleleCounts(int locus)const;
  
  void UpdateAlleleCounts(int locus, const int h[2], const int ancestry[2], bool diploid, bool anneal );
#ifdef PARALLEL
  void SumAlleleCountsOverProcesses(unsigned K);
  void BroadcastAlleleFreqs();
#endif
  void ResetSumAlleleFreqs();
  void setAlleleFreqsMAP();

//   float getEtaRWSamplerAcceptanceRate(int k)const;
//   float getEtaRWSamplerStepsize(int k)const; 

  float getAlphaSamplerAcceptanceRate(int)const;
  float getAlphaSamplerStepsize(int)const;
  float getEtaSamplerAcceptanceRate(int)const;
  float getEtaSamplerStepsize(int)const;
  void OutputAlleleFreqSamplerAcceptanceRates(const char* filename);

  void resetStepSizeApproximator(int k);

protected:
  int Populations, NumberOfCompositeLoci;
  array_of_allelefreqs Freqs;// allele frequencies
  array_of_allelefreqs AlleleFreqsMAP; // posterior mode of allele freqs
  array_of_allelecounts AlleleCounts;
  array_of_allelecounts hetCounts;//counts of het individuals with distinct ancestry states at SNPs
#ifdef PARALLEL
  int* globalAlleleCounts;
  int* globalHetCounts;
  double * globalFreqs;
  int startlocus;
  int* LociCounts;
  int* displs;
  //MPI_Aint stride;//stride of Freqs array for datatype definition
  //MPI::Datatype AlleleFreqArrayType;//datatype for a 2d array
#endif
  int worker_rank;
  int NumWorkers;
  double **PriorParams;
  int FREQSAMPLER;// 1 = conjugate sampler, 2 = Hamiltonian sampler
  std::vector<AlleleFreqSampler*> FreqSampler;
  bool RandomAlleleFreqs;//indicator for whether allele freqs are fixed or random
  bool hapmixmodel;

  Genome *Loci;//pointer to Loci object
  std::ofstream allelefreqoutput;// object to output allele frequencies

  void LoadAlleleFreqs(AdmixOptions* const options, InputData* const data, LogWriter &Log);
  virtual void OpenOutputFile(const char* filename);
  virtual void SampleAlleleFreqs(int, const double coolness);

private:

  double *eta; //dispersion parameter
  double *SumEta;
  double *psi,*tau;// eta has Gamma prior with shape and scale parameters psi and tau
  double psi0;
  MuSampler *muSampler;
  //new sampler for eta
  DispersionSampler *EtaSampler;
 
  double **HistoricAlleleFreqs;

  double **HistoricAlleleCounts;
  bool calculateFST;
  double** Fst;
  double** SumFst;
  bool IsHistoricAlleleFreq;//indicator for dispersion model

  bool CorrelatedAlleleFreqs;

  ///adaptive RW sampler for eta
  StepSizeTuner *TuneEtaSampler;
  int NumberOfEtaUpdates;
  int *NumberAccepted; 
  double *etastep;
  double etastep0;
  int w;//The eta sampler is tuned every w updates.

  //sampler for univariate mu (beta proportion parameters)
  //    DARS SampleMu;
  // DARS SampleMu is now initialized each time the allele freqs at a locus in a population are sampled

  // sampler for Dirichlet proportion parameters 
  std::vector<StepSizeTuner> *MuProposal;
  
  std::ofstream outputstream;//outputs eta to paramfile
  std::ofstream fstoutputstream;

  void OpenFSTFile(const AdmixOptions* const options, LogWriter &Log); 

  void LoadInitialAlleleFreqs(const char* filename, LogWriter &Log);
  void LoadAlleleFreqs(const Matrix_s& NewFreqs, int i, unsigned row0, bool);
  void SetDefaultAlleleFreqs(int i);
  void SetDefaultPriorParams(int i, double defaultpriorparams);

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
