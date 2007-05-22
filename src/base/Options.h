// *-*-C++-*-*
/* 
 *   Options.h 
 *   header file for Options base class
 *   Copyright (c) 2006, 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef OPTIONS_H
#define OPTIONS_H 1

#include "common.h"
#include "bcppcl/LogWriter.h"
#include "bcppcl/OptionReader.h"

/// Class to hold program options
class Options: public bcppcl::OptionReader
{
public:
  virtual ~Options();
  virtual bool SetOptions();
  bool ReadUserOptions(const char* fileargIndicator = 0);
  virtual int checkOptions(LogWriter &Log, int NumberOfIndividuals);
  ///print user options to file
  virtual void PrintUserOptions(const char* filename);

  long getBurnIn() const;
  long getTotalSamples() const;
  long getSampleEvery() const;
  long getNumAnnealedRuns() const; 
  long getSeed() const; 
  bool getThermoIndicator() const;    
  int getDisplayLevel() const;
  bool CheckData()const;
  virtual int getPopulations() const = 0;
  virtual void setPopulations(int) = 0;

  //main output files
  const string getResultsDir() const;
  const char *getLogFilename() const;
  const char *getErgodicAverageFilename() const;
  const char *getParameterFilename() const;
  const char *getRegressionOutputFilename() const;
  const char* getEYFilename()const;
  const char *getAlleleFreqOutputFilename() const;
  bool getOutputAlleleFreq() const;

  //input file names
  const char *getLocusFilename() const;
  const char *getGenotypesFilename() const;
  const char *getOutcomeVarFilename() const;  
  const char *getCoxOutcomeVarFilename() const;  
  int getTargetIndicator() const;
  const char *getCovariatesFilename() const;  
  const char *getPriorAlleleFreqFilename() const;

  double getRegressionPriorPrecision()const;
  unsigned int getgenotypesSexColumn() const;
  void setgenotypesSexColumn(unsigned int i);

  //indicators and model options
  virtual bool getHapMixModelIndicator() const = 0;
  virtual bool isRandomMatingModel() const = 0;
  virtual bool isGlobalRho() const = 0;
  virtual bool getFixedAlleleFreqs() const;
  int getNumberOfOutcomes() const;
  void setNumberOfOutcomes(int);
  void setRegType(RegressionType R);
  virtual bool getTestOneIndivIndicator() const = 0;
  bool getDeleteOldResultsIndicator()const;
  bool doPrintBuildInfo()const;

  //Score test file names
  const char *getAllelicAssociationScoreFilename() const; 
  const char* getResidualAllelicAssocScoreFilename()const;
  const char *getHWTestFilename() const;

  //score test indicators 
  bool getTestForAllelicAssociation() const;
  void setTestForAllelicAssociation(bool); 
  bool getTestForResidualAllelicAssoc()const; 
  bool getHWTestIndicator() const;

  //filenames in derived classes
  //these are required by InputData. Can be removed if InpuData and ScoreTets are forked successfully
  virtual const char* getCCGenotypesFilename()const{return "";};
  virtual const char* getAlleleFreqFilename ()const{return "";};
  virtual const char *getHistoricalAlleleFreqFilename() const{return "";};
  virtual const char *getReportedAncestryFilename() const{return "";};
  virtual const char *getEtaPriorFilename() const{return "";};
  virtual const char *getAffectedsOnlyScoreFilename() const{return "";};
  virtual const char *getHaplotypeAssociationScoreFilename() const{return "";};
  virtual const char *getAncestryAssociationScoreFilename() const{return "";};
  virtual const char *getAlleleFreqScoreFilename() const{return "";};
  virtual const char *getAlleleFreqScoreFilename2() const{return "";};
  virtual const char *getAssocScoreFilename() const{return "";};
  virtual const char* getLikRatioFilename() const{return "";};
  virtual const char* getIndAdmixModeFilename()const{return "";};
  virtual bool getTestForAdmixtureAssociation() const{return false;};
  virtual bool getTestForAffectedsOnly() const{return false;};
  virtual void setTestForAffectedsOnly(bool){};  
  virtual bool getTestForHaplotypeAssociation() const{return false;};
  virtual void setTestForHaplotypeAssociation(bool){};
  virtual bool getTestForLinkageWithAncestry() const{return false;};
  virtual void setTestForLinkageWithAncestry(bool){};
  virtual bool getTestForMisspecifiedAlleleFreqs() const{return false;};
  virtual bool getTestForMisspecifiedAlleleFreqs2() const{return false;};


protected:
  long burnin;
  long TotalSamples;
  long SampleEvery;
  long Seed;
  bool thermoIndicator;//calculate marginal likelihood using simulated annealing
  long NumAnnealedRuns;
  bool checkData;
  int displayLevel;
  int NumberOfOutcomes;
  RegressionType RegType;
  int TargetIndicator;
  bool fixedallelefreqs;
  bool OutputAlleleFreq;
  bool TestForAllelicAssociation;
  bool TestForResidualAllelicAssoc;
  bool HWTest;
  double regressionPriorPrecision;
  bool DeleteOldResultsIndicator;//indicates whether to delete contents of resultsdir
  bool PrintBuildInfo;

  string ResultsDir;
  string LogFilename;
  string LocusFilename;
  string GenotypesFilename;
  unsigned int genotypesSexColumn;
  string AlleleFreqOutputFilename;
  string ErgodicAverageFilename;
  string ParameterFilename;
  string RegressionOutputFilename;
  string AllelicAssociationScoreFilename;
  string ResidualAllelicAssocScoreFilename;
  string HWTestFilename;
  string PriorAlleleFreqFilename;
  string CovariatesFilename;
  string OutcomeVarFilename;
  string CoxOutcomeVarFilename;
  string EYFilename;

  std::vector<float> rhoSamplerParams;//parameters for sampler of population sumintensities or arrival rate

  virtual void SetDefaultValues();  
  virtual void DefineOptions();
  Options();
  bool ReadUserOptions(int, char**, const char* fileargIndicator = 0);
private:

  // UNIMPLEMENTED: to avoid use
  Options(const Options&);
  Options& operator=(const Options&);
};

#endif /* OPTIONS_H */
