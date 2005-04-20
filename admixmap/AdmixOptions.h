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

  const char *getResultsDir() const;
  
  const char *getAffectedsOnlyScoreFilename() const;
  const char *getAlleleFreqFilename() const;
  
  long getBurnIn() const;
  
  const char *getDICoutputFilename() const;
  const char *getDispersionTestFilename() const;
  
  bool getTestForDispersion() const;
  bool getStratificationTest() const;
  void setStratificationTest(bool);
  bool getFixedAlleleFreqs() const;
  
  Vector_i getStratificationTestLoci() const;
  
  const char *getErgodicAverageFilename() const;
  
  const char *getGeneInfoFilename() const;
  const char *getGeneticDataFilename() const;
  const char *getIndAdmixtureFilename() const;
  
  const char *getOutputFilename1() const;
  const char *getOutputFilename2() const;
  const char *getOutputFilename3() const;
  
  const char *getParameterFilename() const;
  const char *getRegressionOutputFilename() const;
  const char *getEtaOutputFilename() const;
  const char *getReportedAncestryFilename() const;
  
  long getTotalSamples() const;
  int getTargetIndicator() const;
  
  int useCOUT() const;
  
  const char *getAlleleFreqOutputFilename() const;
  bool getOutputAlleleFreq() const;

  const char *getFSTOutputFilename() const;
  
  bool getOutputFST() const;
  
  const char *getAlleleFreqScoreFilename() const;
  const char *getAlleleFreqScoreFilename2() const;
  
  int getAnalysisTypeIndicator() const;
  
  const char *getAssocScoreFilename() const;
  
  long getSampleEvery() const;    
  const char *getHistoricalAlleleFreqFilename() const;
  const char *getInitialParameterFilename() const;
  
  const char *getInputFilename() const;  
  const char *getMLEFilename() const;
  
  bool getLocusForTestIndicator() const;
  int getLocusForTest() const;
  
  const char *getLocusScoreFilename() const; 
  const char *getAncestryAssociationScoreFilename() const;
  const char *getLogFilename() const;  
  
  bool getModelIndicator() const;
  bool getRhoIndicator() const;
  bool getIndAdmixHierIndicator() const;
  bool getMLIndicator() const;
  double getTruncPt() const;  
  
  const char *getPosteriorSummariesFilename() const; 
  
  int getPopulations() const; 
  void setPopulations(int num);
  
  const char *getPrecisionFilename() const;    
  const char *getPredictionFilename() const;
  const char *getPriorAlleleFreqFilename() const;
  
  bool getScoreTestIndicator() const;  
  long getSeed() const;  
  double getRho() const;  
  
  const char *getTargetFilename() const;  
  
  bool getTestForAffectedsOnly() const;
  void setTestForAffectedsOnly(bool);  
  bool getTestForAllelicAssociation() const;
  void setTestForAllelicAssociation(bool); 
  bool getTestForSNPsInHaplotype() const;
  bool getTestForLinkageWithAncestry() const;
  void setTestForLinkageWithAncestry(bool);   
  bool getTestForMisspecifiedAlleleFreqs() const;
  bool getTestForMisspecifiedAlleleFreqs2() const;
  
  const char *getTestInputFilename() const;  
  const char *getTestsForSNPsInHaplotypeOutputFilename() const;
  const char *getEtaPriorFilename() const;
  const char *getTestTargetFilename() const;  
  
  int getTextIndicator() const;
  
  Vector_d getInitAlpha(int) const;
  int sizeInitAlpha() const;
  
  unsigned int IsPedFile() const;
  void IsPedFile(unsigned int);
  
  unsigned int genotypesSexColumn() const;    
  void genotypesSexColumn(unsigned int i);
  
  bool getXOnlyAnalysis() const;
  
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
