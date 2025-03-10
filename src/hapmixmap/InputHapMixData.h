// *-*-C++-*-*
/**
 *   HAPMIXMAP 
 *   InputHapMixData.h 
 *   Class to read HAPMIXMAP data
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef INPUTHAPMIXDATA_H
#define INPUTHAPMIXDATA_H

#include "InputData.h"
#include "HapMixGenotypeLoader.h"

class HapMixOptions;

///Class to read HAPMIXMAP data
class InputHapMixData : public InputData{
public:
  /**
     Constructor
     \param options pointer to options object
     \param Log LogWriter, for mesages
  */
  InputHapMixData(HapMixOptions *options, bclib::LogWriter &log); 

  ///check contents of data files, after they are read in
  bool CheckData(HapMixOptions *options, bclib::LogWriter &Log);

  /**
     Retrieves an individual's genotypes (obsolete).
     retrieves ADMIXMAP-style genotypes.
     \see InputAdmixData::GetGenotype
  */ 
  void GetGenotype(int i, const Genome &Loci, std::vector<genotype>* genotypes, bool **Missing)const;

  /**
     Retrieves an individual's genotypes.
     \param i Individual number (count from 1)
     \param Loci Genome object, to count through loci
     \param genotypes pointer to individual's genotypes
     \param Mising individual's array of missing genotype indicators
     \return true if individual is haploid, false if diploid
     \see HapMixGenotypeLoader::GetHapMixGenotype
  */
  bool GetHapMixGenotype(int i, const Genome &Loci, std::vector<unsigned short>* genotypes, bool** Missing);

  /**
     Determines if an Individual is under test or not.
     \param i Individual number (count from 1)
     \return true if under test, false if not
  */
  bool IsTestIndividual(int i)const;

  /**
     Determines if a locus is typed (case/control)
     \param locus index of locus
     \return true if locus is typed, false if not
  */
  bool isTypedLocus(unsigned locus)const;

  ///returns number of typed loci
  unsigned getNumTypedLoci()const;

  ///returns the number of test individuals
  int getNumberOfTestIndividuals()const;

private:

  ///object to read and assign hapmix genotypes
  HapMixGenotypeLoader* hGenotypeLoader;

  ///check outcome variables and covariates and extract required variables 
  void CheckRegressionData(HapMixOptions *options, bclib::LogWriter &Log);

  /** read block state labels.
      Read from priorallelefreqfile, if specified. 
      Otherwise assign default names - BlockState1, BlockState2 etc
  */
  void ReadBlockStateLabels(HapMixOptions *options);

  /**
     Checks contents of priorallelefreqfile
  */
  bool CheckAlleleFreqs(HapMixOptions *options, bclib::LogWriter &Log);
  /*
   *  UNIMPLEMENTED: to avoid undesired copying.
   */
  ///default constructor, not implemented
  InputHapMixData();
  ///copy constructor, not implemented    
  InputHapMixData(const InputHapMixData&);
  ///assignment operator, not implemented
  void operator=(const InputHapMixData&);

};


#endif
