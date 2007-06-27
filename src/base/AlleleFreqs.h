// *-*-C++-*-*
/** 
 *   AlleleFreqs.h 
 *   header file for AlleleFreqs class
 *   Copyright (c) 2005 - 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef ALLELEFREQS_H
#define ALLELEFREQS_H 1

#define FREQ_CONJUGATE_SAMPLER 1
#define FREQ_HAMILTONIAN_SAMPLER 2

#include "InputData.h"
#include "Genome.h"
#include "Options.h"
#include "AlleleFreqSampler.h"
#include "common.h"
#include "FreqArrays.h"
#include "bcppcl/RObjectWriter.h"

class LogWriter;

/// Class to hold allele/haplotype frequencies and their priors.
class AlleleFreqs{

public:
  AlleleFreqs();
  virtual ~AlleleFreqs();
  //virtual void Initialise(Options* const options, InputData* const Data, Genome *pLoci, LogWriter &Log, bool MAP=false);
  void LoadInitialAlleleFreqs(const char* filename, LogWriter &Log);
  void AllocateAlleleCountArrays(unsigned K);
  virtual void PrintPrior(const Vector_s&, LogWriter& Log)const;
  virtual void Update(IndividualCollection*IC , bool afterBurnIn, double coolness);

  ///outputs ergodic averages of dispersion parameters (SumEta)  to ErgodicAverageFile
  virtual void OutputErgodicAvg( int iteration, std::ofstream *avgstream)const;

  void OutputAlleleFreqs();
  void OutputAlleleFreqs(const char* filename, LogWriter& Log);
  void CloseOutputFile(int iterations, const Vector_s& PopulationLabels);

  ///resets Allelecounts to zero at start of iteration
  void ResetAlleleCounts(unsigned K);
  bool IsRandom()const;

  const double *GetStatsForEta( int , int locus)const;
  double GetAlleleProbsMAP( int x, int ancestry , int locus)const;
  std::vector<double> GetPriorAlleleFreqs( int locus, int population )const;
  std::vector<int> GetAlleleCounts( int locus, int population )const;
  std::vector<double> getAlleleFreqsMAP( int locus, int population )const;
  std::vector<double> GetAlleleFreqs( int locus, int population )const;
  const double *GetAlleleFreqs(int locus)const;
  const FreqArray& GetAlleleFreqs()const;
  const int *GetAlleleCounts(int locus)const;
  
  void UpdateAlleleCounts(int locus, const int h[2], const int ancestry[2], bool diploid, bool anneal );
  void ResetSumAlleleFreqs();
  void setAlleleFreqsMAP();

//   float getEtaRWSamplerAcceptanceRate(int k)const;
//   float getEtaRWSamplerStepsize(int k)const; 

  void OutputAlleleFreqSamplerAcceptanceRates(const char* filename);

  virtual void resetStepSizeApproximator(int k);

  void WriteKLInfo(unsigned samples, ostream& os);
  void WriteLocusInfo(unsigned samples, const std::string& ResultsDir, const std::vector<std::string>& PopLabels);

protected:
  int Populations, NumberOfCompositeLoci;
  FreqArray Freqs;// allele frequencies
  FreqArray AlleleFreqsMAP; // posterior mode of allele freqs
  array_of_allelecounts AlleleCounts;
  array_of_allelecounts hetCounts;//counts of het individuals with distinct ancestry states at SNPs
  int worker_rank;
  int NumWorkers;
  double **PriorParams;
  int FREQSAMPLER;// 1 = conjugate sampler, 2 = Hamiltonian sampler
  std::vector<AlleleFreqSampler*> FreqSampler;
  bool RandomAlleleFreqs;//indicator for whether allele freqs are fixed or random
  bool hapmixmodel;

  Genome *Loci;//pointer to Loci object
  RObjectWriter allelefreqoutput;// object to output allele frequencies
  float* SumKLInfo;// to accumulate Kullback-Liebler info
  float** SumLocusInfo;//to accumulate locus information content (f)

  virtual void LoadAlleleFreqs(const Matrix_s& NewFreqs, int i, unsigned row0, bool);
  virtual void OpenOutputFile(const char* filename);
  virtual void SampleAlleleFreqs(int, const double coolness);
  void SetDefaultAlleleFreqs(int i);
  virtual void SetDefaultPriorParams(int i, double defaultpriorparams);

  void Initialise(bool OutputFreqs);
  void AccumulateKLInfo();
  void AccumulateLocusInfo();
private:
  //void LoadAlleleFreqs(Options* const options, InputData* const data, LogWriter &Log);

  static double muEnergyFunction(unsigned K, const double * const alpha, const double* const *args);
  static void muGradient(unsigned K, const double * const alpha, const double* const *args, double *g);

};

#endif
