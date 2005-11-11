// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   AdmixOptions.h 
 *   header file for AdmixOptions class
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
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

#ifndef ADMIX_OPTIONS_H
#define ADMIX_OPTIONS_H 1

#include "common.h"
#include <getopt.h>    /* for getopt and getopt_long */
#include "LogWriter.h"
#include <string.h>
#include <string>
#include <map>
using namespace::std;

typedef map<const char*, char*> OptionMap;

class AdmixOptions
{
public:
  AdmixOptions();
  AdmixOptions(int, char**);
  ~AdmixOptions();
  
  void SetOptions(int, char**);
  int checkOptions(LogWriter *Log, int NumberOfIndividuals);
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
  const char *getMLEFilename() const;
  const char *getHistoricalAlleleFreqFilename() const;
  const char *getPriorAlleleFreqFilename() const;
  const char *getAlleleFreqFilename() const;
  const char *getReportedAncestryFilename() const;
  const char *getEtaPriorFilename() const;
  const char *getOutcomeVarFilename() const;  
  int getTargetIndicator() const;
  double getRho() const;
  double getRhoalpha() const;
  double getRhobeta() const;
  double getRhobetaShape()const;
  double getRhobetaRate()const;
  bool RhoFlatPrior() const;
  bool logRhoFlatPrior() const;  

  vector<double> getInitAlpha(int) const;
  std::vector<std::vector<double> > getInitAlpha()const;
  int sizeInitAlpha() const;
  double getAlphamean()const;
  double getAlphavar() const;
  double getEtaMean() const;
  double getEtaVar() const;

  unsigned int getgenotypesSexColumn() const;
  void setgenotypesSexColumn(unsigned int i);
  
  //indicators and model options  
  int useCOUT() const;
  bool getFixedAlleleFreqs() const;
  bool getCorrelatedAlleleFreqs() const;
  bool isRandomMatingModel() const;
  int getNumberOfOutcomes() const;
  void setNumberOfOutcomes(int);
  void setRegType(RegressionType R);
  bool isGlobalRho() const;
  bool getIndAdmixHierIndicator() const;
  bool getMLIndicator() const;
  bool getAnnealIndicator() const;
  int getNumberOfAnnealedRuns() const;
  double getTruncPt() const;  
  int getPopulations() const; 
  void setPopulations(int num);
  bool IsPedFile() const;
  void IsPedFile(bool);
  bool getXOnlyAnalysis() const;
  bool isAdmixed(unsigned) const;
  bool isSymmetric()const;

  //Score test file names
  const char *getAffectedsOnlyScoreFilename() const;
  const char *getTestsForSNPsInHaplotypeOutputFilename() const;
  const char *getAllelicAssociationScoreFilename() const; 
  const char *getAncestryAssociationScoreFilename() const;
  const char *getAlleleFreqScoreFilename() const;
  const char *getAlleleFreqScoreFilename2() const;
  const char *getAssocScoreFilename() const;
  const char *getHWTestFilename() const;
  const char* getLikRatioFilename() const;

  //score test indicators 
  bool getScoreTestIndicator() const; //indicator for any score test (except misspec allelefreqs) 
  bool getTestForAdmixtureAssociation() const;
  bool getTestForAffectedsOnly() const;
  void setTestForAffectedsOnly(bool);  
  bool getTestForAllelicAssociation() const;
  void setTestForAllelicAssociation(bool); 
  bool getTestForSNPsInHaplotype() const;
  void setTestForSNPsInHaplotype(bool);
  bool getTestForLinkageWithAncestry() const;
  void setTestForLinkageWithAncestry(bool);   
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
  //int AnalysisTypeIndicator;
  int NumberOfOutcomes;
  RegressionType RegType;
  int TargetIndicator;
  double TruncPt;
  int Populations;

  int use_cout;
  bool OutputFST;
  bool XOnlyAnalysis;
  bool isPedFile;
  unsigned int genotypesSexColumn;
  bool locusForTestIndicator;
  int LocusForTest;
  bool fixedallelefreqs;
  bool correlatedallelefreqs;
  bool RandomMatingModel;//random mating model
  bool GlobalRho;//indicator for global rho
  bool IndAdmixHierIndicator;//hierarchical model on ind admixture
  bool MLIndicator;//calculate marginal likelihood 
  int AnnealedRuns;
  bool ScoreTestIndicator; //indicator for any of the score tests in ScoreTests class
  bool TestForAdmixtureAssociation;
  bool StratificationTestIndicator;
  bool TestForAffectedsOnly;
  bool TestForAllelicAssociation;
  bool TestForSNPsInHaplotype;
  bool TestForDispersion;
  bool TestForLinkageWithAncestry;
  bool TestForMisspecifiedAlleleFreqs;
  bool TestForMisspecifiedAlleleFreqs2;
  bool HWTest;
  bool OutputAlleleFreq;

  std::vector<bool> _admixed;
  bool _symmetric;         

  //priors
  double Rhoalpha, Rhobeta;//gamma parameters for sumintensities
  double RhobetaShape, RhobetaRate;//gamma parameters for prior on rhobeta
  std::vector<double> alpha0;
  std::vector<double> alpha1;
  std::vector< std::vector<double> > initalpha;
  double alphamean, alphavar;
  double etamean, etavar;//gamma parameters for dispersion parameter

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
  string TestsForSNPsInHaplotypeOutputFilename;
  string AllelicAssociationScoreFilename;
  string AncestryAssociationScoreFilename;
  string HWTestFilename;
  string LikRatioFilename;

  string LocusFilename;
  string GenotypesFilename;
  string HistoricalAlleleFreqFilename;
  string PriorAlleleFreqFilename;
  string CovariatesFilename;
  string OutcomeVarFilename;
  string MLEFilename;
  string EtaPriorFilename;
  string ReportedAncestryFilename;

  OptionMap OptionValues;//to output user options
  
  void Initialise();  
  void SetOutputNames();
  void setInitAlpha(LogWriter *Log);
  bool CheckInitAlpha( const std::vector<double> &alphatemp)const;

  // UNIMPLIMENTED: to avoid use
  AdmixOptions(const AdmixOptions&);
  AdmixOptions& operator=(const AdmixOptions&);
};

#endif /* ADMIX_OPTIONS_H */
