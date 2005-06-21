// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   InputData.h 
 *   header file for InputData class
 *   Copyright (c) 2005 LSHTM
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
    

  /**
   *  Read all input data and store in internal structures.
   */    
  void readData(AdmixOptions *options, LogWriter *log);    
  int getNumberOfIndividuals();
  int getNumberOfSimpleLoci(); 

  /**
   *  Getters to retrieve data (in string form).
   */
  const Matrix_s& getGeneInfoData() const;
  const Matrix_s& getGeneticData()  const;
  const Matrix_s& getInputData()    const;
  const Matrix_s& getTargetData()   const;
  const Matrix_s& getAlleleFreqData() const;
  const Matrix_s& getHistoricalAlleleFreqData() const;
  const Matrix_s& getPriorAlleleFreqData() const;
  const Matrix_s& getEtaPriorData() const;
  const Matrix_s& getMLEData() const;
  const Matrix_s& getReportedAncestryData() const;

  /**
   *  Getters to retrieve data (converted to Matrix_d).
   */    
  const Matrix_d& getGeneInfoMatrix() const;
  const Matrix_d& getTargetMatrix() const;
  const Matrix_d& getInputMatrix() const;
  const Matrix_d& getAlleleFreqMatrix() const;
  const Matrix_d& getHistoricalAlleleFreqMatrix() const;
  const Matrix_d& getPriorAlleleFreqMatrix() const;
  const Matrix_d& getEtaPriorMatrix() const;
  const Matrix_d& getMLEMatrix() const;
  const Matrix_d& getReportedAncestryMatrix() const;

  bool determineIfPedFile(AdmixOptions *options);
  void convertGenotypesToIntArray(AdmixOptions *options);
  void convertToVectorsOverCLoci(Genome & Loci, Chromosome **chrm);

  int GetSexValue(int i);
 
  void GetGenotype(int i, int SexColumn, Genome &Loci,unsigned short ****genotype);
  void CheckAlleleFreqs(AdmixOptions *options, int NumberOfCompositeLoci, int NumberOfStates);
private:    
  Matrix_s geneInfoData_;
  Matrix_s geneticData_;
  Matrix_s inputData_;
  Matrix_s targetData_;
  Matrix_s alleleFreqData_;
  Matrix_s historicalAlleleFreqData_;
  Matrix_s priorAlleleFreqData_;
  Matrix_s etaPriorData_;
  Matrix_s MLEData_;
  Matrix_s reportedAncestryData_;

  Matrix_d geneInfoMatrix_;
  Matrix_d inputMatrix_;
  Matrix_d targetMatrix_;
  Matrix_d alleleFreqMatrix_;
  Matrix_d historicalAlleleFreqMatrix_;
  Matrix_d priorAlleleFreqMatrix_;
  Matrix_d etaPriorMatrix_;
  Matrix_d MLEMatrix_;
  Matrix_d reportedAncestryMatrix_;

  int NumIndividuals;
  int NumSimpleLoci;
  int NumCompositeLoci;
  bool IsPedFile;
  LogWriter *Log;

  Matrix_g genotypes_g; // 3-way array of genotypes: indivs, simple loci, alleles
  
  std::vector<std::vector< std::vector<unsigned int> > > genotypes_c; 
  // 3-way array of genotypes: indivs, composite loci, 2*simple loci    


private:
  /**
   *  UNIMPLEMENTED: to avoid undesired copying.
   */    
  InputData(const InputData&);
  void operator=(const InputData&);

  void readFile(const char *fname, Matrix_s& data);
  void CheckGeneticData(int genotypesSexColumn);
  void checkLociNames(AdmixOptions *options);
  void CheckOutcomeVarFile(bool);
  void CheckCovariatesFile();
  void CheckRepAncestryFile(int populations);
  void throwGenotypeError(int ind, int locus, std::string label, int g0, int g1, int numalleles);
};

#endif /* !defined INPUT_DATA_H */
