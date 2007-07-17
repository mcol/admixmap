// *-*-C++-*-*
/** 
 *   InputData.h 
 *   header file for InputData class
 *   Copyright (c) 2005 - 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef INPUT_DATA_H
#define INPUT_DATA_H 1
#include "common.h"
#include "bclib/DataMatrix.h"
#include "GeneticDistanceUnit.h"
#include "GenotypeLoader.h"
#include <vector>

/**
 *  Forward declarations.
 */    
class Options;
class LogWriter;
class Genome;
class Chromosome;

///Class to read and check all input data files
class InputData{
public:    

  /**
   *  Constructor.
   */    
  InputData();
    
  /**
   *  Destructor.
   */    
  virtual ~InputData();
    
  ///erases (nearly) all memory used by object
  virtual void Delete();

  /*
   *  Getters to retrieve data (in string form).
   */
  ///return contents of locusfile
  const Matrix_s& getLocusData() const;
  ///return contents of covariatesfile
  const Matrix_s& getCovariatesData()    const;
  ///return contents of outcomevarfile
  const Matrix_s& getOutcomeVarData()   const;
  ///return contents of priorallelefreqfile
  const Matrix_s& getPriorAlleleFreqData() const;

  /*
   *  Getters to retrieve data (converted to DataMatrix).
   */    
  const DataMatrix& getLocusMatrix() const;
  const DataMatrix& getOutcomeVarMatrix() const;
  const DataMatrix& getCoxOutcomeVarMatrix() const;
  const DataMatrix& getCovariatesMatrix() const;

  ///
  void getOutcomeTypes(DataType*)const;
  ///returns the DataType of an outcome variable
  DataType getOutcomeType(unsigned)const;
  ///returns population/hidden state labels
  const Vector_s& GetHiddenStateLabels() const;
  ///returns outcome variable labels
  Vector_s getOutcomeLabels()const;
  ///returns locus labels (from locusfile)
  const Vector_s& getLocusLabels()const;
  ///returns covariete labels
  const Vector_s getCovariateLabels()const;

  /**
     Determines if an individual is female.
     Determined from sex column of genotypesfile, if there is one.
     \param i Individual number (count from 1)
     \return true if female, false if male or unknown
  */
  bool isFemale(int i)const;

  ///returns number of individuals
  int getNumberOfIndividuals()const;

  ///returns total number of loci in locusfile
  int getNumberOfSimpleLoci()const;
  ///returns number of composite loci 
  unsigned getNumberOfCompositeLoci()const{return NumCompositeLoci;};
  ///returns unit of distance used in locusfile (defaults to Morgans)
  GeneticDistanceUnit getUnitOfDistance()const;
  ///returns a string representing the unit of distance in locusfile
  const std::string& getUnitOfDistanceAsString()const;
  float getLocusDistanceThreshold()const;

protected:
  /** Extracts population labels from header line of allelefreq file.
      \param data header line of file.
      \param Populations number of populations/hidden states
      \param labels vector to fill with labels
  */
  void getPopLabels(const Vector_s& data, size_t Populations, Vector_s& labels);

  /**
     Checks contents of locusfile and extracts locus labels.
     If check = true, checks that loci labels in locusfile are unique and that they match the names in the genotypes file.
     \param sexColumn the index of the sex column in genotypesfile if there is one; 0 is not.
     \param threshold maximum allowable distance between loci on a chromosome. If a distance exceeds this, the chromosome is broken up and a warning is printed.
     \param check determines whether to check locus labels
  */
  bool checkLocusFile(LogWriter& Log);
  void SetLocusLabels();

  ///determine number of composite loci from locusfile
  unsigned determineNumberOfCompositeLoci()const;
  ///determine unit of distance in locusfile
  GeneticDistanceUnit DetermineUnitOfDistance();
  //virtual void CheckAlleleFreqs(Options *options, LogWriter &Log);
  ///check contents of outcomevarfile
  void CheckOutcomeVarFile(unsigned N, Options * const options, LogWriter &Log);
  ///check contsnts of coxoutcomevarfile
  void CheckCoxOutcomeVarFile(LogWriter &log)const;
  ///check contents of covariatesfile
  void CheckCovariatesFile(unsigned NumIndividuals, Options* const options, LogWriter &Log);
  /**
   *  Read input data and store in internal structures.
   */    
  void ReadData(Options *options, LogWriter &log);    

  Matrix_s locusData_;
  Matrix_s covariatesData_;
  Matrix_s outcomeVarData_;
  Matrix_s coxOutcomeVarData_;
  Matrix_s priorAlleleFreqData_;

  DataMatrix locusMatrix_;
  DataMatrix covariatesMatrix_;
  DataMatrix outcomeVarMatrix_;
  DataMatrix coxOutcomeVarMatrix_;
  DataMatrix priorAlleleFreqMatrix_;

  GenotypeLoader* genotypeLoader;
  int NumSimpleLoci;
  unsigned NumCompositeLoci;
  GeneticDistanceUnit distanceUnit;
  Vector_s HiddenStateLabels;
  std::vector<DataType> OutcomeType;
  Vector_s CovariateLabels;
  Vector_s OutcomeLabels;

private:    

  Vector_s LocusLabels;

  //void CheckData(Options *options, LogWriter &Log);

  /*
   *  UNIMPLEMENTED: to avoid undesired copying.
   */
  ///copy constructor, not implemented    
  InputData(const InputData&);
  //assignment operator, not implemented
  void operator=(const InputData&);

};

void getLabels(const Vector_s& data, std::string *labels);
#endif /* !defined INPUT_DATA_H */
