// *-*-C++-*-*
/* 
 *   Options.h 
 *   header file for Options base class
 */

/*
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
#include "bclib/LogWriter.h"
#include "bclib/OptionReader.h"


/** \addtogroup base
 * @{ */



/// Class to hold program options
class Options: public bclib::OptionReader
{
private:
  bool usePedForInd	; ///< Use Pedigree rather than Individual objects for individuals (single-member pedigrees).
  bool warningsAreErrors; ///< Abort execution without running model if input data contains warnings.
  int  maxCPUsToUse	; ///< Maximum CPUs (cores) to use in parallel.  0 (default) for all available.
  unsigned maxPedigreeSize; ///< Maximum number of organisms in a pedigree (larger pedigrees will be reduced by removing non-founders).

  /// When parsing pedigree files, turn parent-IDs that refer to non-existent
  /// individuals into unknown-parent.  Otherwise there are errors.
  bool ignoreInvalidParents;

  bool excludeMendelError;
  bool excludeUnaffectedSibs;

public:
  virtual ~Options();
  virtual bool SetOptions();
  bool ReadUserOptions(const char* fileargIndicator = 0);
  virtual int checkOptions(bclib::LogWriter &Log, int NumberOfIndividuals);
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
  virtual bool outputParams()const = 0;
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
  const char *getCovariatesFilename() const;  
  const std::vector<unsigned>& getOutcomeVarColumns()const;
  const std::vector<unsigned>& getCovariateColumns()const;
  const char *getPriorAlleleFreqFilename() const;

  double getRegressionPriorPrecision()const;

  //indicators and model options
  virtual bool getHapMixModelIndicator() const = 0;
  virtual bool isRandomMatingModel() const = 0;
  virtual bool isGlobalRho() const = 0;
  virtual bool getFixedAlleleFreqs() const;
  int getNumberOfOutcomes() const;
  void setNumberOfOutcomes(unsigned);
  void setRegType(RegressionType R);
  virtual bool getTestOneIndivIndicator() const = 0;
  bool getDeleteOldResultsIndicator()const;

  //score test indicators 
  bool getTestForAllelicAssociation() const;
  void setTestForAllelicAssociation(bool); 
  bool getTestForResidualAllelicAssoc()const; 
  bool getHWTestIndicator() const;

  const genepi::cvector<float> & getrhoSamplerParams()const;

protected:
  long burnin;
  long TotalSamples;
  long SampleEvery;
  long Seed;
  bool thermoIndicator;//calculate marginal likelihood using simulated annealing
  long NumAnnealedRuns;
  bool checkData;
  int displayLevel;
  RegressionType RegType;
  std::vector<unsigned> OutcomeVarColumns;//columns in outcomevarfile to be used
  std::vector<unsigned> CovariateColumns;//columns in covariatesfile to be used
  bool fixedallelefreqs;
  bool OutputAlleleFreq;
  bool TestForAllelicAssociation;
  bool TestForResidualAllelicAssoc;
  bool HWTest;
  double regressionPriorPrecision;
  bool DeleteOldResultsIndicator;//indicates whether to delete contents of resultsdir

  string ResultsDir;
  string LogFilename;
  string LocusFilename;
  string GenotypesFilename;
  string AlleleFreqOutputFilename;
  string ErgodicAverageFilename;
  string ParameterFilename;
  string RegressionOutputFilename;
  string PriorAlleleFreqFilename;
  string CovariatesFilename;
  string OutcomeVarFilename;
  string CoxOutcomeVarFilename;
  string EYFilename;
  unsigned NumberOfOutcomes;

  genepi::cvector<float> rhoSamplerParams;//parameters for sampler of population sumintensities or arrival rate

  virtual void SetDefaultValues();  
  virtual void DefineOptions();
  Options();
  bool ReadUserOptions(int, char**, const char* fileargIndicator = 0);
private:

  // UNIMPLEMENTED: to avoid use
  Options(const Options&);
  Options& operator=(const Options&);


  public:

    bool getUsePedForInd     () const { return usePedForInd	; } ///< Use Pedigree rather than Individual objects for
								    ///< unrelated individuals (i.e. single-member pedigrees).
    bool getWarningsAreErrors() const { return warningsAreErrors; } ///< Abort execution without running model if
								    ///< the input data generates and warnings.
    int	 getMaxCPUsToUse     () const { return maxCPUsToUse	; } ///< Maximum CPUs (cores) to use in parallel.
								    ///< 0 (default) for all available.
    unsigned getMaxPedigreeSize() const { return maxPedigreeSize  ; } ///< Maximum number of organisms in a pedigree (larger
								    ///< pedigrees will be reduced by removing non-founders).

    /// When parsing pedigree files, turn parent-IDs that refer to non-existent
    /// individuals into unknown-parent.  Otherwise there are errors.
    bool getIgnoreInvalidParents() const { return ignoreInvalidParents; }

    /// Should pedigrees with Mendelian inconsistencies be excluded from the input dataset?
    bool getExcludeMendelError() const { return excludeMendelError; }

    /// Should exclude unaffected non-founders?
    bool getExcludeUnaffectedSibs() const { return excludeUnaffectedSibs; }

};


/** @} */


#endif /* OPTIONS_H */
