// *-*-C++-*-*
/* 
 *   HAPMIXMAP
 *   HapmixOptions.h 
 *   header file for AdmixOptions class
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef HAPMIX_OPTIONS_H
#define HAPMIX_OPTIONS_H 1

#include "Options.h"
using namespace::std;

//TODO:check dispparamfile, allelefreqoutputfile, allelefreqprior, globalrhoprior



/// Class to hold program options
class HapMixOptions : public Options
{
public:
  HapMixOptions(int, char**);
  ~HapMixOptions();
  
  int checkOptions(LogWriter &Log, int NumberOfIndividuals);
  void PrintOptions();

  //main output files
  const char *getEtaOutputFilename() const;
  const char *getAlleleFreqPriorOutputFilename() const;
  bool OutputAlleleFreqPrior()const;
  const char* getHapMixLambdaOutputFilename()const;
  const char* getFinalLambdaFilename()const;
  const char* getFinalFreqPriorFilename()const;

  const std::vector<double> &getHapMixLambdaPrior()const;
  const std::vector<double> & getAlleleFreqPriorParams()const;
  const char* getInitialHapMixLambdaFilename()const;
  const char* getInitialAlleleFreqFilename()const;
  const char* getInitialFreqPriorFilename()const;
  const char* getCCGenotypesFilename()const;

  //indicators and model options
  bool getHapMixModelIndicator() const;
  int getPopulations() const;
  int getNumberOfBlockStates()const;
  void setPopulations(int num);
  bool getFixedAlleleFreqs() const;
  bool isFreqDispersionHierModel()const;
  const vector<float>& getLambdaSamplerParams()const;
  bool getTestOneIndivIndicator() const{return false;};//not supported in hapmixmodel
  bool isRandomMatingModel() const{return false;};//required to pass as Options object to some functions
  bool isGlobalRho() const{return false;};//                 "

  //Score test file names
  const char* getMHTestFilename()const;

  //score test indicators 
  bool getMHTest()const;  
  
  //
  const std::vector<unsigned>& getMaskedIndividuals()const;
  const std::vector<unsigned>& getMaskedLoci()const;
  bool OutputCGProbs()const;
  unsigned GetNumMaskedIndividuals()const;
  unsigned GetNumMaskedLoci()const;

private:
  int NumBlockStates;
  static const bool HapMixModelIndicator = true; //model haplotypes with mixture model
  bool FreqDispersionHierModel;

  std::vector<double> allelefreqprior;
  std::vector<double> hapmixlambdaprior;///< prior means of rho prior params in hapmixmodel

  string EtaOutputFilename;
  string AlleleFreqPriorOutputFilename;
  string MHTestFilename;
  string HapMixLambdaOutputFilename;
  string InitialHapMixLambdaFilename;
  string InitialAlleleFreqFilename;
  string InitialFreqPriorFile;
  string CCGenotypesFilename;//case-control genotypes file (hapmixmodel only)

  string FinalFreqPriorFilename;
  string FinalLambdaFilename;

  //indices for assessing prediction of missing genotypes in hapmixmodel
  std::vector<unsigned> MaskedIndividuals;
  std::vector<unsigned> MaskedLoci;
  
  void SetOptions(OptionMap& ProgOptions);

  void SetDefaultValues();  


  // UNIMPLEMENTED: to avoid use
  HapMixOptions();
  HapMixOptions(const HapMixOptions&);
  HapMixOptions& operator=(const HapMixOptions&);
};

#endif /* HAPMIX_OPTIONS_H */
