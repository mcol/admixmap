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
    ~AdmixOptions();

  void SetOptions(int, char**);
  int checkOptions(LogWriter *Log);
  void PrintOptions();

  long getBurnIn() const;
  long getTotalSamples() const;
  long getSampleEvery() const; 
  int getAnalysisTypeIndicator() const;
  long getSeed() const;     

  //main output files
  const char *getResultsDir() const;
  const char *getLogFilename() const;
  const char *getErgodicAverageFilename() const;
  const char *getParameterFilename() const;
  const char *getRegressionOutputFilename() const;
  const char *getEtaOutputFilename() const;
  const char *getIndAdmixtureFilename() const;
  const char *getAlleleFreqOutputFilename() const;
  bool getOutputAlleleFreq() const;

  //input file names
  const char *getGeneInfoFilename() const;
  const char *getGeneticDataFilename() const;
  const char *getInputFilename() const;  
  const char *getMLEFilename() const;
  const char *getHistoricalAlleleFreqFilename() const;
  const char *getPriorAlleleFreqFilename() const;
  const char *getAlleleFreqFilename() const;
  const char *getReportedAncestryFilename() const;
  const char *getEtaPriorFilename() const;
  const char *getTargetFilename() const;  
  int getTargetIndicator() const;
  double getRho() const;  

  Vector_d getInitAlpha(int) const;
  int sizeInitAlpha() const;

  unsigned int getgenotypesSexColumn() const;
  void setgenotypesSexColumn(unsigned int i);
  
  //indicators and model options  
  int useCOUT() const;
  bool getFixedAlleleFreqs() const;
  int getTextIndicator() const;
  bool isRandomMatingModel() const;
  bool getRhoIndicator() const;
  bool getIndAdmixHierIndicator() const;
  bool getMLIndicator() const;
  bool getAnnealIndicator() const;
  double getTruncPt() const;  
  int getPopulations() const; 
  void setPopulations(int num);
  unsigned int IsPedFile() const;
  void IsPedFile(unsigned int);
  bool getXOnlyAnalysis() const;  

  //Score test file names
  const char *getAffectedsOnlyScoreFilename() const;
  const char *getTestsForSNPsInHaplotypeOutputFilename() const;
  const char *getAllelicAssociationScoreFilename() const; 
  const char *getAncestryAssociationScoreFilename() const;
  const char *getAlleleFreqScoreFilename() const;
  const char *getAlleleFreqScoreFilename2() const;
  const char *getAssocScoreFilename() const;
  const char *getHWTestFilename() const;

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
  const char *getDICoutputFilename() const;
  const char *getDispersionTestFilename() const;
  const char *getFSTOutputFilename() const;
  
  //other test indicators
  bool getTestForDispersion() const;
  bool getStratificationTest() const;
  void setStratificationTest(bool);
  bool getOutputFST() const;
  bool getLocusForTestIndicator() const;
  int getLocusForTest() const;

  Vector_i getStratificationTestLoci() const;
  
  
private:
  long burnin;
  long TotalSamples;
  long SampleEvery;
  long Seed;
  int AnalysisTypeIndicator;
  int TargetIndicator;
  double TruncPt;
  int Populations;

  int use_cout;
  int TextIndicator;
  bool OutputFST;
  bool XOnlyAnalysis;
  unsigned int isPedFile;
  unsigned int genotypesSexColumn;
  bool locusForTestIndicator;
  int LocusForTest;
  bool fixedallelefreqs;
  bool RandomMatingModel;//random mating model
  bool RhoIndicator;//indicator for non-global rho
  bool IndAdmixHierIndicator;//hierarchical model on ind admixture
  bool MLIndicator;//calculate marginal likelihood - valid only for analysistypeindicator < 0
  bool AnnealIndicator;
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

  double Rho;
  std::vector<Vector_d> alpha;

  string ResultsDir;
  string LogFilename;
  string AffectedsOnlyScoreFilename;
  string AlleleFreqOutputFilename;
  string AlleleFreqScoreFilename;
  string AlleleFreqScoreFilename2;
  string AssocScoreFilename;
  string alleleFreqFilename;
  string DICoutputFilename;
  string ErgodicAverageFilename;
  string ParameterFilename;
  string RegressionOutputFilename;
  string EtaOutputFilename;
  string DispersionTestFilename;
  string IndAdmixtureFilename;
  string FSTOutputFilename;
  string TestsForSNPsInHaplotypeOutputFilename;
  string AllelicAssociationScoreFilename;
  string AncestryAssociationScoreFilename;
  string HWTestFilename;

  string GeneInfoFilename;
  string GeneticDataFilename;
  string HistoricalAlleleFreqFilename;
  string PriorAlleleFreqFilename;
  string InputFilename;
  string TargetFilename;
  string MLEFilename;
  string EtaPriorFilename;
  string ReportedAncestryFilename;

  OptionMap OptionValues;//to output user options
  
  // UNIMPLIMENTED: to avoid use
  AdmixOptions(const AdmixOptions&);
  AdmixOptions& operator=(const AdmixOptions&);
  
  void SetOutputNames();
  };

#endif /* ADMIX_OPTIONS_H */
