// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   InputData.h 
 *   header file for InputData class
 *   Copyright (c) 2005, 2006 David O'Donnell and Paul McKeigue
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
#include "DataMatrix.h"

/**
 *  Forward declarations.
 */    
class AdmixOptions;
class LogWriter;
class Chromosome;
class Genome;

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
  ~InputData();
    
  void Delete();//erases (nearly) all memory used by object
  /**
   *  Read all input data and store in internal structures.
   */    
  void readData(AdmixOptions *options, LogWriter &log, int rank);    


  /*
   *  Getters to retrieve data (in string form).
   */
  const Matrix_s& getLocusData() const;
  const Matrix_s& getGeneticData()  const;
  const Matrix_s& getInputData()    const;
  const Matrix_s& getOutcomeVarData()   const;
  const Matrix_s& getAlleleFreqData() const;
  const Matrix_s& getHistoricalAlleleFreqData() const;
  const Matrix_s& getPriorAlleleFreqData() const;
  const Matrix_s& getEtaPriorData() const;
  const Matrix_s& getReportedAncestryData() const;

  /*
   *  Getters to retrieve data (converted to DataMatrix).
   */    
  const DataMatrix& getLocusMatrix() const;
  const DataMatrix& getOutcomeVarMatrix() const;
  const DataMatrix& getCoxOutcomeVarMatrix() const;
  const DataMatrix& getCovariatesMatrix() const;
  //const DataMatrix& getAlleleFreqMatrix() const;
  //const DataMatrix& getHistoricalAlleleFreqMatrix() const;
  //const DataMatrix& getPriorAlleleFreqMatrix() const;
  const DataMatrix& getEtaPriorMatrix() const;
  const DataMatrix& getReportedAncestryMatrix() const;

  static void convertMatrix(const Matrix_s& data, DataMatrix& m, size_t row0, size_t col0, size_t ncols);

  void getOutcomeTypes(DataType*)const;
  DataType getOutcomeType(unsigned)const;
  const Vector_s& GetPopLabels() const;
  Vector_s getOutcomeLabels()const;
  const Vector_s& getLocusLabels()const;
  const Vector_s getCovariateLabels()const;
  //TODO: make all above const

  //void convertGenotypesToIntArray(AdmixOptions *options);
  //void convertToVectorsOverCLoci(Genome & Loci, Chromosome **chrm);

  bool isFemale(int i)const;
  int getNumberOfIndividuals()const;
  int getNumberOfSimpleLoci()const;
  unsigned getNumberOfCompositeLoci()const{return NumCompositeLoci;};
  bool distancesAreInCentiMorgans()const;
  void GetGenotype(int i, int SexColumn, const Genome &Loci, std::vector<genotype>* genotypes, bool **Missing)const;

private:    
  Matrix_s locusData_;
  Matrix_s geneticData_;
  Matrix_s inputData_;
  Matrix_s outcomeVarData_;
  Matrix_s coxOutcomeVarData_;
  Matrix_s alleleFreqData_;
  Matrix_s historicalAlleleFreqData_;
  Matrix_s priorAlleleFreqData_;
  Matrix_s etaPriorData_;
  Matrix_s reportedAncestryData_;

  DataMatrix locusMatrix_;
  DataMatrix covariatesMatrix_;
  DataMatrix outcomeVarMatrix_;
  DataMatrix coxOutcomeVarMatrix_;
  DataMatrix alleleFreqMatrix_;
  DataMatrix historicalAlleleFreqMatrix_;
  DataMatrix priorAlleleFreqMatrix_;
  DataMatrix etaPriorMatrix_;
  DataMatrix reportedAncestryMatrix_;

  Vector_s PopulationLabels;
  std::vector<std::string> LocusLabels;
  Vector_s OutcomeLabels;
  std::vector<DataType> OutcomeType;
  Vector_s CovariateLabels;
  int NumIndividuals;
  int NumSimpleLoci;
  unsigned NumCompositeLoci;
  bool IsPedFile;

  void getPopLabels(const Vector_s& data, size_t Populations, Vector_s& labels);
  void readFile(const char *fname, Matrix_s& data, LogWriter &Log);
  void readGenotypesFile(const char *fname, Matrix_s& data);
  void CheckGeneticData(AdmixOptions *options)const;
  void checkLocusFile(int sexColumn, double threshold, bool);
  unsigned determineNumberOfCompositeLoci()const;
  void CheckAlleleFreqs(AdmixOptions *options, LogWriter &Log);
  void CheckOutcomeVarFile(AdmixOptions * const options, LogWriter &Log);
  void CheckCoxOutcomeVarFile(LogWriter &log)const;
  void CheckCovariatesFile(LogWriter &log, bool usePopLabels);
  void CheckRepAncestryFile(int populations, LogWriter &Log)const;
  void throwGenotypeError(int ind, int locus, std::string label, int g0, int g1, int numalleles)const;
  bool determineIfPedFile()const;
  std::vector<unsigned short> GetGenotype(unsigned locus, int individual, int SexColumn)const;
  void CheckData(AdmixOptions *options, LogWriter &Log, int rank);

  /*
   *  UNIMPLEMENTED: to avoid undesired copying.
   */    
  InputData(const InputData&);
  void operator=(const InputData&);

};

#endif /* !defined INPUT_DATA_H */
