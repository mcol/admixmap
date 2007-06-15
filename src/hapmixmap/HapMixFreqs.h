// *-*-C++-*-*
/** 
 *   HAPMIXMAP
 *   HapMixFreqs.h 
 *   header file for HapMixFreqs class, used to sample prior parameters of frequencies in a hapmixmodel
 *   Copyright (c) 2006, 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef HAPMIXFREQS_H
#define HAPMIXFREQS_H 1

#include "AlleleFreqs.h"
#include "bcppcl/AdaptiveRejection.h"
#include <vector>
#include <fstream>
class Genome;
class HapMixOptions;

typedef struct{
    unsigned K;
    double sumlogfreqs1;
    double sumlogfreqs2;
    double eta;

}hapmixmuargs;

/// Class to hold and sample prior parameters of frequencies in a hapmixmodel
class HapMixFreqs : public AlleleFreqs{

public:
  HapMixFreqs();
  ~HapMixFreqs();
  void Initialise(HapMixOptions* const options, InputData* const Data, Genome *pLoci, LogWriter &Log);
  void ReadInitialPriorParamsFromFile(const char* filename, LogWriter& Log);
  void setSampler(bool thermo, bool AllHaploid, bool /*DefaultPriors*/);
  void Update(IndividualCollection*IC , bool afterBurnIn, double coolness, bool AllHaploid);
  void PrintPrior(LogWriter& Log)const;
  void SamplePriorPrecision(unsigned locus, unsigned Populations, double sumlogfreqs1, double sumlogfreqs2);
  void SamplePriorProportions(unsigned locus, double sumlogfreqs1, double sumlogfreqs2);
  void OutputErgodicAvg( int samples, std::ofstream *avgstream)const;
  //void OutputPriorParams();
  void OutputPriorParams(bool tofile, bool toscreen);
  void OutputPosteriorMeans(const char* filename, LogWriter& Log)const;
  void OutputFinalValues(const char* filename, LogWriter& Log)const;
  //double getHapMixPriorRate()const{return HapMixPriorRate;};
  double getParams(unsigned locus)const;
  float getAcceptanceRate()const;
  float getStepSize()const;

  void AllocateDiploidGenotypeProbs();
  void SetDiploidGenotypeProbs();
  const FreqArray& getHaploidGenotypeProbs()const;
  const FreqArray& getDiploidGenotypeProbs()const;

private:
  /**
     In hapmixmodel, prior on allele freqs is Dirichlet with locus-specific mean mu and precision eta.
     Mu and the PriorParams are stored, but not eta.

     the etas each have the same Gamma prior with shape EtaPriorShape and rate EtaPriorRate, specifiable by the user.
     mu has a fixed uniform / beta prior.

     optional: EtaPriorRate has a conjugate gamma prior.
  */

  double* Eta;
  double* DirichletParams;//params of Dirichlet prior on frequencies
  double EtaShape;//params of Gamma prior on Dirichlet params
  double EtaRate;//
  double EtaRatePriorShape;// params of Gamma Prior on params of Gamma prior on params of Dirichlet prior on freqs
  double EtaRatePriorRate;

  bool etaHierModel;//indicates whether to fir a hierarchical model for eta
  bool accumulateEta;//indicates whether to accumulate Eta for output
  unsigned long NumEtaUpdates;
  double* SumEta;//cumulative sum of eta
  double SumLambda;// cumulative sum of EtaRate

  StepSizeTuner* EtaSampler;
  AdaptiveRejection MuSampler;
  hapmixmuargs HapMixMuArgs;

  //TODO: this should be static. No problem in practice as we only have one object of this class
  FreqArray DiploidGenotypeProbs;

  std::ofstream allelefreqprioroutput;//to output mean and variance of frequency prior precision

  void LoadAlleleFreqs(HapMixOptions* const options, InputData* const data_, LogWriter &Log);
  void InitialisePrior(unsigned Populations, unsigned L, const HapMixOptions* const options, LogWriter& Log);
  void OpenOutputFile(const char* filename);
  void SampleAlleleFreqs(int, const double coolness);
  void SampleEtaRate(bool afterburnin, double sum);
  static double fmu_hapmix(double, const void* const);
  static double dfmu_hapmix(double, const void* const);
  static double d2fmu_hapmix(double, const void* const);

};

#endif
