// *-*-C++-*-*
/* 
 *   ADMIXMAP
 *   AdmixOptions.h 
 *   header file for AdmixOptions class
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef ADMIX_OPTIONS_H
#define ADMIX_OPTIONS_H 1

#include "common.h"
#include <getopt.h>    /* for getopt and getopt_long */
#include "LogWriter.h"
#include <string.h>
#include <string>
#include <map>
using namespace::std;

/// map to hold user options for output to file.
typedef map<const char*, char*> OptionMap;

/// Class to hold program options
class AdmixOptions
{
public:
  AdmixOptions();
  AdmixOptions(int, char**);
  ~AdmixOptions();
  
  void SetOptions(int, char**);
  int checkOptions(LogWriter &Log, int NumberOfIndividuals);
  void PrintOptions();

  long getBurnIn() const;
  long getTotalSamples() const;
  long getSampleEvery() const; 
  long getSeed() const;     

  //main output files
  const string getResultsDir() const;
  const char *getLogFilename() const;
  const char *getErgodicAverageFilename() const;
  const char *getParameterFilename() const;
  const char *getRegressionOutputFilename() const;
  const char* getResidualFilename()const;
  const char *getEtaOutputFilename() const;
  const char *getIndAdmixtureFilename() const;
  const char *getAlleleFreqOutputFilename() const;
  bool getOutputAlleleFreq() const;

  //input file names
  const char *getLocusFilename() const;
  const char *getGenotypesFilename() const;
  const char *getCovariatesFilename() const;  
  const char *getHistoricalAlleleFreqFilename() const;
  const char *getPriorAlleleFreqFilename() const;
  const char *getAlleleFreqFilename() const;
  const char *getReportedAncestryFilename() const;
  const char *getEtaPriorFilename() const;
  const char *getOutcomeVarFilename() const;  
  const char *getCoxOutcomeVarFilename() const;  
  int getTargetIndicator() const;
  double getRhoalpha() const;
  double getRhobeta() const;
  double getRhobetaShape()const;
  double getRhobetaRate()const;
  double getRhoPriorMean()const;
  const std::vector<double> & getAlleleFreqPriorParams()const;

  vector<double> getInitAlpha(int) const;
  std::vector<std::vector<double> > getInitAlpha()const;
  int sizeInitAlpha() const;
  double getEtaMean() const;
  double getEtaVar() const;
  double getRegressionPriorPrecision()const;
  const vector<float>& getrhoSamplerParams()const;
  const vector<float>& getrhoPriorParamSamplerParams() const;

  unsigned int getgenotypesSexColumn() const;
  void setgenotypesSexColumn(unsigned int i);
  
  //indicators and model options
  int getDisplayLevel() const;
  bool getFixedAlleleFreqs() const;
  bool getCorrelatedAlleleFreqs() const;
  bool isRandomMatingModel() const;
  int getNumberOfOutcomes() const;
  void setNumberOfOutcomes(int);
  void setRegType(RegressionType R);
  bool isGlobalRho() const;
  bool PopAdmixturePropsAreEqual()const;
  bool getIndAdmixHierIndicator() const;
  bool getHapMixModelIndicator() const;
  bool getChibIndicator() const;
  bool getThermoIndicator() const;
  bool getTestOneIndivIndicator() const;
  long getNumAnnealedRuns() const;
  int getPopulations() const; 
  void setPopulations(int num);
  bool isAdmixed(unsigned) const;
  bool isSymmetric()const;
  bool CheckData()const;

  //Score test file names
  const char *getAffectedsOnlyScoreFilename() const;
  const char *getHaplotypeAssociationScoreFilename() const;
  const char *getAllelicAssociationScoreFilename() const; 
  const char *getAncestryAssociationScoreFilename() const;
  const char *getAlleleFreqScoreFilename() const;
  const char *getAlleleFreqScoreFilename2() const;
  const char *getAssocScoreFilename() const;
  const char* getResidualAllelicAssocScoreFilename()const;
  const char *getHWTestFilename() const;
  const char* getLikRatioFilename() const;
  const char* getIndAdmixModeFilename()const;

  //score test indicators 
  bool getScoreTestIndicator() const; //indicator for any score test (except misspec allelefreqs) 
  bool getTestForAdmixtureAssociation() const;
  bool getTestForAffectedsOnly() const;
  void setTestForAffectedsOnly(bool);  
  bool getTestForAllelicAssociation() const;
  void setTestForAllelicAssociation(bool); 
  bool getTestForHaplotypeAssociation() const;
  void setTestForHaplotypeAssociation(bool);
  bool getTestForLinkageWithAncestry() const;
  void setTestForLinkageWithAncestry(bool);
  bool getTestForResidualAllelicAssoc()const;   
  bool getTestForMisspecifiedAlleleFreqs() const;
  bool getTestForMisspecifiedAlleleFreqs2() const;
  bool getHWTestIndicator() const;
  
  //other test file names
  const char *getStratTestFilename() const;
  const char *getDispersionTestFilename() const;
  const char *getFSTOutputFilename() const;
  
  //other test indicators
  bool getTestForDispersion() const;
  bool getStratificationTest() const;
  void setStratificationTest(bool);
  bool getOutputFST() const;
  bool getLocusForTestIndicator() const;
  int getLocusForTest() const;
  
private:
  long burnin;
  long TotalSamples;
  long SampleEvery;
  long Seed;
  int NumberOfOutcomes;
  RegressionType RegType;
  int TargetIndicator;
  int Populations;
  int displayLevel;
  bool OutputFST;
  unsigned int genotypesSexColumn;
  bool locusForTestIndicator;
  int LocusForTest;
  bool fixedallelefreqs;
  bool correlatedallelefreqs;
  bool RandomMatingModel;//random mating model
  bool GlobalRho;//indicator for global rho
  bool IndAdmixHierIndicator;//hierarchical model on ind admixture
  bool HapMixModelIndicator; //model haplotypes with mixture model
  bool chibIndicator;//calculate marginal likelihood using Chib method
  bool thermoIndicator;//calculate marginal likelihood using simulated annealing
  bool TestOneIndivIndicator;//calculate marginal likelihood for one individual only
  bool PopAdmixPropsAreEqual;
  long NumAnnealedRuns;
  bool ScoreTestIndicator; //indicator for any of the score tests in ScoreTests class
  bool TestForAdmixtureAssociation;
  bool StratificationTestIndicator;
  bool TestForAffectedsOnly;
  bool TestForAllelicAssociation;
  bool TestForHaplotypeAssociation;
  bool TestForResidualAllelicAssoc;
  bool TestForDispersion;
  bool TestForLinkageWithAncestry;
  bool TestForMisspecifiedAlleleFreqs;
  bool TestForMisspecifiedAlleleFreqs2;
  bool HWTest;
  bool OutputAlleleFreq;
  bool checkData;

  std::vector<bool> _admixed;
  bool _symmetric;         

  //priors
  //double Rhoalpha, Rhobeta;//gamma parameters for sumintensities
  //double RhobetaShape, RhobetaRate;//gamma parameters for prior on rhobeta
  std::vector<double> globalrhoPrior;
  std::vector<double> rhoPrior;
  bool hapmixrhoparamprior;// indicates whether globalrho
  std::vector<double> alpha0;
  std::vector<double> alpha1;
  std::vector< std::vector<double> > initalpha;
  std::vector<double> allelefreqprior;
  double etamean, etavar;//gamma parameters for dispersion parameter
  double regressionPriorPrecision;

  std::vector<float> rhoSamplerParams;
  std::vector<float> rhoPriorParamSamplerParams;
  string ResultsDir;
  string LogFilename;
  string AffectedsOnlyScoreFilename;
  string AlleleFreqOutputFilename;
  string AlleleFreqScoreFilename;
  string AlleleFreqScoreFilename2;
  string AssocScoreFilename;
  string alleleFreqFilename;
  string StratTestFilename;
  string ErgodicAverageFilename;
  string ParameterFilename;
  string RegressionOutputFilename;
  string EtaOutputFilename;
  string DispersionTestFilename;
  string ResidualFilename;
  string IndAdmixtureFilename;
  string FSTOutputFilename;
  string HaplotypeAssociationScoreFilename;
  string AllelicAssociationScoreFilename;
  string AncestryAssociationScoreFilename;
  string ResidualAllelicAssocScoreFilename;
  string HWTestFilename;
  string LikRatioFilename;
  string IndAdmixModeFilename;

  string LocusFilename;
  string GenotypesFilename;
  string HistoricalAlleleFreqFilename;
  string PriorAlleleFreqFilename;
  string CovariatesFilename;
  string OutcomeVarFilename;
  string CoxOutcomeVarFilename;
  string EtaPriorFilename;
  string ReportedAncestryFilename;

  OptionMap OptionValues;//to output user options
  
  void Initialise();  
  void SetOutputNames();
  void setInitAlpha(LogWriter &Log);
  bool CheckInitAlpha( const std::vector<double> &alphatemp)const;

  // UNIMPLEMENTED: to avoid use
  AdmixOptions(const AdmixOptions&);
  AdmixOptions& operator=(const AdmixOptions&);
};

#endif /* ADMIX_OPTIONS_H */
