// *-*-C++-*-*
/* 
 *   HAPMIXMAP
 *   HapMixGenotypeLoader.h 
 *   class to load and assign genotypes in hapmixmodel
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "GenotypeLoader.h"

///class to read genotypes from file and load into Individuals' arrays
class HapMixGenotypeLoader : public GenotypeLoader{

public:
  ///constructor
  HapMixGenotypeLoader();
  ///read case-control genotypesfile
  void ReadCaseControlGenotypes(const char* filename, bclib::LogWriter& Log);

  /**
     Retrieve HAPMIXMAP-style genotype.
     Extracts genotypes read as strings from file, converts to 0, 1, 2(haploid) or 0, 1, 2, 3, 4 (diploid). 
     0 denotes missing and assigns to genotype vector.
     \param i individual number (count from 1)
     \param Loci Genome object, to count through loci
     \param genotypes pointer to vector to store genotypes
     \param Missing missing-genotype indicators to be filled
     \return true if individual is haploid, false if diploid
  */
  bool GetHapMixGenotype(int i, const Genome &Loci,
			 std::vector<unsigned short>* genotypes, bool** Missing);

  ///returns number of individuals (including cases and controls)
  unsigned getNumberOfIndividuals()const;
  ///returns number of cases and controls
  unsigned getNumberOfCaseControlIndividuals()const;
  ///returns the number of typed loci in a hapmix case-control analysis
  unsigned getNumTypedLoci()const;
  ///returns number of case-control loci (columns in ccgenotypesfile)
  unsigned NumCaseControlLoci()const;
  ///determines if an individual is a case-control or not
  bool IsCaseControl(unsigned i)const;
  ///determines if a locus is typed (in case-control file) or not
  bool isTypedLocus(unsigned locus)const;
  ///free memory
  void clear();

  /**
     get ADMIXMAP-style genotype (obsolete).
     \see GenotypeLoader::GetGenotype
  */
  void GetGenotype(int i, const Genome &Loci, 
		   std::vector<genotype>* genotypes, bool **Missing)const;
private:
  ///string matrix to store case-control genotypes for hapmixmodel
  Matrix_s CCgeneticData_;
  ///vector of indicators of typed loci
  std::vector<bool> isCaseControlSNP;
  ///number of cases/controls
  int NumCCIndividuals;

  /**
     gets a hapmix case-control genotype from the ccgenotypes file.
     \param locus locus index
     \param cclocus case-control locus index (incremented if the locus is typed)
     \param individual individual number (count from 1)
     \return vector of length 1(haploid) or 2(diploid) of alles
  */
  std::vector<unsigned short> GetCaseControlGenotype(unsigned locus, unsigned* cclocus, 
						     int individual)const;
  ///determine which loci in genotypesfile are also in case-control genotypes file
  void FindCaseControlLoci();
  /// retrieve ADMIXMAP-style case-control genotype (obsolete)
  void GetCaseControlGenotype(int i, const Genome &Loci, std::vector<genotype>* genotypes, bool** Missing)const;
};
