// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   HapMIxFreqs.h 
 *   header file for HapMixFreqs class, used to sample prior parameters of frequencies in a hapmixmodel
 *   Copyright (c) 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef HAPMIXFREQS_H
#define HAPMIXFREQS_H 1

#include "utils/LogWriter.h"
#include "samplers/StepSizeTuner.h"
#include "samplers/AdaptiveRejection.h"
#include <vector>
#include <fstream>
class Genome;

typedef struct{
    unsigned K;
    double sumlogfreqs1;
    double sumlogfreqs2;
    double eta;

}hapmixmuargs;

/// Class to hold and sample prior parameters of frequencies in a hapmixmodel
class HapMixFreqs{

public:
  HapMixFreqs();
  ~HapMixFreqs();
  void Initialise(unsigned Populations, unsigned L, const std::vector<double> &params);
  void PrintPrior(LogWriter& Log)const;
  void OpenOutputFile(const char* filename);
  void SamplePriorDispersion(unsigned locus, unsigned Populations, double sumlogfreqs1, double sumlogfreqs2);
  void SamplePriorProportions(unsigned locus, double sumlogfreqs1, double sumlogfreqs2);
  void OutputErgodicAvg( int samples, std::ofstream *avgstream)const;
  void OutputPriorParams();
  void OutputPriorParams(std::ostream& os, bool tofile);
  //double getHapMixPriorRate()const{return HapMixPriorRate;};
  double getParams(unsigned locus)const;
  float getAcceptanceRate()const;
  float getStepSize()const;

private:
  unsigned NumberOfCompositeLoci;
  double* HapMixPriorEta;
  double* HapMixPriorParams;//params of Dirichlet prior on frequencies
  double HapMixPriorShape;//params of Gamma prior on Dirichlet params
  double HapMixPriorRate;//
  double HapMixPriorRatePriorShape;// params of Gamma Prior on params of Gamma prior on params of Dirichlet prior on freqs
  double HapMixPriorRatePriorRate;
  double SumLambda;// cumulative sum of HapMixPriorRatePriorRate

  StepSizeTuner* HapMixPriorEtaSampler;
  AdaptiveRejection HapMixPriorMuSampler;
  hapmixmuargs HapMixMuArgs;

  std::ofstream allelefreqprioroutput;//to output mean and variance of frequency prior dispersion in hapmixmodel

  static double fmu_hapmix(double, const void* const);
  static double dfmu_hapmix(double, const void* const);
  static double d2fmu_hapmix(double, const void* const);

};

#endif
