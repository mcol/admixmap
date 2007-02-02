// *-*-C++-*-*
/* 
 *   ADMIXMAP
 *   AdmixOptions.h 
 *   header file for AdmixOptions class
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef ADMIX_OPTIONS_H
#define ADMIX_OPTIONS_H 1

#include "Options.h"
using namespace::std;

/// Class to hold program options
class AdmixOptions : public Options
{
public:
  AdmixOptions(int, char**);
  ~AdmixOptions();
  
  int checkOptions(LogWriter &Log, int NumberOfIndividuals);
  void PrintOptions();

  //main output files
  const char *getEtaOutputFilename() const;
  const char *getIndAdmixtureFilename() const;

  //input file names
  const char *getAlleleFreqFilename() const;
  const char *getHistoricalAlleleFreqFilename() const;
  const char *getReportedAncestryFilename() const;
  const char *getEtaPriorFilename() const;

  double getRhoalpha() const;
  double getRhobeta() const;
  double getRhobetaShape()const;
  double getRhobetaRate()const;
  double getRhoPriorMean()const;

  vector<double> getInitAlpha(int) const;
  std::vector<std::vector<double> > getInitAlpha()const;
  int sizeInitAlpha() const;
  double getEtaMean() const;
  double getEtaVar() const;
  const vector<float>& getrhoSamplerParams()const;
  const vector<float>& getPopAdmixSamplerParams()const;

  //indicators and model options
  bool getCorrelatedAlleleFreqs() const;
  bool isRandomMatingModel() const;
  bool isGlobalRho() const;
  bool PopAdmixturePropsAreEqual()const;
  bool getIndAdmixHierIndicator() const;
  bool getHapMixModelIndicator() const;
  bool getChibIndicator() const;
  bool getTestOneIndivIndicator() const;
  bool isAdmixed(unsigned) const;
  bool isSymmetric()const;
  int getPopulations() const;
  void setPopulations(int num);
  bool getFixedAlleleFreqs() const;

  //Score test file names
  const char *getAffectedsOnlyScoreFilename() const;
  const char *getHaplotypeAssociationScoreFilename() const;
  const char *getAncestryAssociationScoreFilename() const;
  const char *getAlleleFreqScoreFilename() const;
  const char *getAlleleFreqScoreFilename2() const;
  const char *getAssocScoreFilename() const;
  const char* getLikRatioFilename() const;
  const char* getIndAdmixModeFilename()const;

  //score test indicators 
  bool getTestForAdmixtureAssociation() const;
  bool getTestForAffectedsOnly() const;
  void setTestForAffectedsOnly(bool);  
  bool getTestForHaplotypeAssociation() const;
  void setTestForHaplotypeAssociation(bool);
  bool getTestForLinkageWithAncestry() const;
  void setTestForLinkageWithAncestry(bool);
  bool getTestForMisspecifiedAlleleFreqs() const;
  bool getTestForMisspecifiedAlleleFreqs2() const;
  
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
  int Populations;
  bool OutputFST;
  bool locusForTestIndicator;
  int LocusForTest;
  bool correlatedallelefreqs;
  bool RandomMatingModel;//random mating model
  bool GlobalRho;//indicator for global rho
  bool IndAdmixHierIndicator;//hierarchical model on ind admixture
  static const bool HapMixModelIndicator = false; //model haplotypes with mixture model
  bool chibIndicator;//calculate marginal likelihood using Chib method
  bool TestOneIndivIndicator;//calculate marginal likelihood for one individual only
  bool PopAdmixPropsAreEqual;
  bool TestForAdmixtureAssociation;
  bool StratificationTestIndicator;
  bool TestForAffectedsOnly;
  bool TestForHaplotypeAssociation;
  bool TestForDispersion;
  bool TestForLinkageWithAncestry;
  bool TestForMisspecifiedAlleleFreqs;
  bool TestForMisspecifiedAlleleFreqs2;

  std::vector<bool> _admixed;
  bool _symmetric;         

  //priors
  //double Rhoalpha, Rhobeta;//gamma parameters for sumintensities
  //double RhobetaShape, RhobetaRate;//gamma parameters for prior on rhobeta
  std::vector<double> globalrhoPrior;
  std::vector<double> rhoPrior;
  std::vector<double> hapmixlambdaprior;///< prior means of rho prior params in hapmixmodel
  std::vector<double> alpha0;
  std::vector<double> alpha1;
  std::vector< std::vector<double> > initalpha;
  double etamean, etavar;//gamma parameters for dispersion parameter

  std::vector<float> popAdmixSamplerParams;//parameters for sampler of population admixture
  string AffectedsOnlyScoreFilename;
  //string AlleleFreqPriorOutputFilename;
  string AlleleFreqScoreFilename;
  string AlleleFreqScoreFilename2;
  string AssocScoreFilename;
  string StratTestFilename;
  string EtaOutputFilename;
  string DispersionTestFilename;
  string IndAdmixtureFilename;
  string FSTOutputFilename;
  string HaplotypeAssociationScoreFilename;
  string AncestryAssociationScoreFilename;
  string LikRatioFilename;
  string IndAdmixModeFilename;
  string alleleFreqFilename;  
  string HistoricalAlleleFreqFilename;
  string EtaPriorFilename;
  string ReportedAncestryFilename;
  
  void SetOptions(OptionMap& ProgOptions);

  void SetDefaultValues();  
  void setInitAlpha(LogWriter &Log);
  bool CheckInitAlpha( const std::vector<double> &alphatemp)const;

  // UNIMPLEMENTED: to avoid use
  AdmixOptions();
  AdmixOptions(const AdmixOptions&);
  AdmixOptions& operator=(const AdmixOptions&);
};

#endif /* ADMIX_OPTIONS_H */
