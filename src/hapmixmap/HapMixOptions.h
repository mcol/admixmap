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

//TODO:check allelefreqoutputfile, allelefreqprior, globalrhoprior



/// Class to hold program options
class HapMixOptions : public Options
{
public:
  HapMixOptions(int, char**, bool PrintOptionList);
  ~HapMixOptions();
  
  int checkOptions(LogWriter &Log, int NumberOfIndividuals);
  void PrintUserOptions();

  //main output files
  const char *getFreqPrecisionOutputFilename() const;
  const char *getAlleleFreqPriorOutputFilename() const;
  bool OutputAlleleFreqPrior()const;
  const char* getArrivalRateOutputFilename()const;
  const char* getFinalMixturePropsFilename()const;
  const char* getFinalLambdaFilename()const;
  const char* getFinalFreqPriorFilename()const;

  const std::vector<double>& getMixturePropsPrecisionPrior()const;
  const std::vector<double> &getLambdaPrior()const;
  const std::vector<double> & getAlleleFreqPriorParams()const;
  const char* getInitialArrivalRateFilename()const;
  const char* getInitialMixturePropsFilename()const;
  const char* getInitialAlleleFreqFilename()const;
  const char* getInitialFreqPriorFilename()const;
  const char* getCCGenotypesFilename()const;

  //indicators and model options
  bool getHapMixModelIndicator() const;
  int getPopulations() const;
  int getNumberOfBlockStates()const;
  void setPopulations(int num);
  bool getFixedAlleleFreqs() const;
  bool isFreqPrecisionHierModel()const;
  bool getFixedMixtureProps()const;
  bool getFixedMixturePropsPrecision()const;
  float getMixturePropsPrecision()const;
  const std::vector<float>& getLambdaSamplerParams()const;
  bool getTestOneIndivIndicator() const{return false;};//not supported in hapmixmodel
  bool isRandomMatingModel() const{return false;};//required to pass as Options object to some functions
  bool isGlobalRho() const{return false;};//                 "

  //Score test file names
  const char* getMHTestFilename()const;

  //score test indicators 
  bool getMHTest()const;  
  
  const std::vector<unsigned>& getMaskedIndividuals()const;
  const std::vector<unsigned>& getMaskedLoci()const;
  bool OutputCGProbs()const;
  unsigned GetNumMaskedIndividuals()const;
  unsigned GetNumMaskedLoci()const;

private:
  //model settings
  static const bool HapMixModelIndicator = true; 
  int NumBlockStates;
  bool FreqPrecisionHierModel;
  bool FixedMixtureProps;
  bool FixedMixturePropsPrecision;
  float MixturePropsPrecision;

  //prior parameters
  std::vector<double> allelefreqprecisionprior;
  std::vector<double> lambdaprior;///< parameters of gamma priors on arrival rate distribution
  std::vector<double> MixturePropsPrecisionPrior;// parameters of Gamma prior on mixture props dispersion

  //data files
  std::string CCGenotypesFilename;//case-control genotypes file

  //parameter output
  std::string FreqPrecisionOutputFilename;

  //posterior means
  std::string AlleleFreqPriorOutputFilename;
  std::string ArrivalRateOutputFilename;

  //tests
  std::string MHTestFilename;

  //intitial values
  std::string InitialArrivalRateFilename;
  std::string InitialMixturePropsFilename;
  std::string InitialAlleleFreqFilename;
  std::string InitialFreqPriorFile;

  //final values
  std::string FinalFreqPriorFilename;
  std::string FinalLambdaFilename;
  std::string FinalMixturePropsFilename;
  //allelefreqoutputfile used in lieu of FinalAlleleFreqFile

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
