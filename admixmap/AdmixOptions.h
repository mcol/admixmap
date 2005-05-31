// *-*-C++-*-*
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

  unsigned int genotypesSexColumn() const;    
  void genotypesSexColumn(unsigned int i);
  
  //indicators and model options  
  int useCOUT() const;
  bool getFixedAlleleFreqs() const;
  int getTextIndicator() const;
  bool isRandomMatingModel() const;
  bool getRhoIndicator() const;
  bool getIndAdmixHierIndicator() const;
  bool getMLIndicator() const;
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
  
  
private: // members
  
  struct Options;     // class declaration only
  Options *imp;       // implementation
  OptionMap OptionValues;
  
  // UNIMPLIMENTED: to avoid use
  AdmixOptions(const AdmixOptions&);
  AdmixOptions& operator=(const AdmixOptions&);
  
  void SetOutputNames();
  };

#endif /* ADMIX_OPTIONS_H */
