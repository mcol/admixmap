// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   InputData.h 
 *   header file for InputData class
 *   Copyright (c) 2005, 2006 LSHTM
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
class InputData
{
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
  void readData(AdmixOptions *options, LogWriter &log);    


  /**
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
  const Matrix_s& getMLEData() const;
  const Matrix_s& getReportedAncestryData() const;

  /**
   *  Getters to retrieve data (converted to DataMatrix).
   */    
  const DataMatrix& getLocusMatrix() const;
  const DataMatrix& getOutcomeVarMatrix() const;
  const DataMatrix& getCovariatesMatrix() const;
  const DataMatrix& getAlleleFreqMatrix() const;
  const DataMatrix& getHistoricalAlleleFreqMatrix() const;
  const DataMatrix& getPriorAlleleFreqMatrix() const;
  const DataMatrix& getEtaPriorMatrix() const;
  const DataMatrix& getMLEMatrix() const;
  const DataMatrix& getReportedAncestryMatrix() const;

  void getOutcomeTypes(DataType*)const;
  std::string *GetPopLabels() const;
  Vector_s getOutcomeLabels()const;
  std::vector<std::string> getLocusLabels()const;

  bool determineIfPedFile(AdmixOptions *options)const;
  void convertGenotypesToIntArray(AdmixOptions *options);
  void convertToVectorsOverCLoci(Genome & Loci, Chromosome **chrm);

  Sex GetSexValue(int i)const;
  int getNumberOfIndividuals()const;
  int getNumberOfSimpleLoci()const;
  unsigned getNumberOfCompositeLoci()const{return NumCompositeLoci;};  

  void GetGenotype(int i, int SexColumn, const Genome &Loci, std::vector<genotype>* genotypes)const;

private:    
  Matrix_s locusData_;
  Matrix_s geneticData_;
  Matrix_s inputData_;
  Matrix_s outcomeVarData_;
  Matrix_s alleleFreqData_;
  Matrix_s historicalAlleleFreqData_;
  Matrix_s priorAlleleFreqData_;
  Matrix_s etaPriorData_;
  Matrix_s MLEData_;
  Matrix_s reportedAncestryData_;

  DataMatrix locusMatrix_;
  DataMatrix covariatesMatrix_;
  DataMatrix outcomeVarMatrix_;
  DataMatrix alleleFreqMatrix_;
  DataMatrix historicalAlleleFreqMatrix_;
  DataMatrix priorAlleleFreqMatrix_;
  DataMatrix etaPriorMatrix_;
  DataMatrix MLEMatrix_;
  DataMatrix reportedAncestryMatrix_;

  std::string* PopulationLabels;
  std::vector<std::string> LocusLabels;
  Vector_s OutcomeLabels;
  DataType* OutcomeType;
  int NumIndividuals;
  int NumSimpleLoci;
  unsigned NumCompositeLoci;
  bool IsPedFile;

  void readFile(const char *fname, Matrix_s& data, LogWriter &Log);
  void CheckGeneticData(AdmixOptions *options)const;
  void checkLocusFile(int sexColumn);
  unsigned determineNumberOfCompositeLoci()const;
  void CheckAlleleFreqs(AdmixOptions *options, LogWriter &Log);
  RegressionType CheckOutcomeVarFile(int, int, LogWriter &Log);
  void CheckCovariatesFile(LogWriter &log)const;
  void CheckRepAncestryFile(int populations, LogWriter &Log)const;
  void throwGenotypeError(int ind, int locus, std::string label, int g0, int g1, int numalleles)const;

  /**
   *  UNIMPLEMENTED: to avoid undesired copying.
   */    
  InputData(const InputData&);
  void operator=(const InputData&);

};

#endif /* !defined INPUT_DATA_H */
