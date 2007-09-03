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
  void ReadTestGenotypes(const char* filename, bclib::LogWriter& Log);

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

  ///returns number of individuals (including typed)
  unsigned getNumberOfIndividuals()const;
  ///returns number of test individuals
  unsigned getNumberOfTestIndividuals()const;
  ///returns the number of typed loci in a hapmix test analysis (columns in testgenotypesfile)
  unsigned getNumTypedLoci()const;

  ///determines if an individual is under test or not
  bool IsTestIndividual(unsigned i)const;
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
  Matrix_s testGeneticData_;
  ///vector of indicators of typed loci
  std::vector<bool> isTypedSNP;
  ///number of cases/controls
  int NumTestIndividuals;

  /**
     gets a hapmix test genotype from the testgenotypes file.
     \param locus locus index
     \param test test locus index (incremented if the locus is typed)
     \param individual individual number (count from 1)
     \return vector of length 1(haploid) or 2(diploid) of alles
  */
  std::vector<unsigned short> GetTestGenotype(unsigned locus, unsigned* testlocus, 
					      int individual)const;
  ///determine which loci in genotypesfile are also in test genotypes file
  void FindTypedLoci();
  /// retrieve ADMIXMAP-style case-control genotype (obsolete)
  void GetTestGenotype(int i, const Genome &Loci, std::vector<genotype>* genotypes, bool** Missing)const;
};
