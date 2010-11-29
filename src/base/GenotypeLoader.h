//=============================================================================
//
// Copyright (C) 2007  David O'Donnell and Paul McKeigue
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file GenotypeLoader.h
/// Definition of the GenotypeLoader class.
//=============================================================================

#ifndef GENOTYPELOADER_H
#define GENOTYPELOADER_H


#warning GenotypeLoader is deprecated -- replaced by GenotypeParser


#include <vector>
#include <string>
#include "common.h"



/** \addtogroup base
 * @{ */



class Genome;
namespace bclib{
  class DataMatrix;
  class LogWriter;
}

typedef std::vector<std::vector<unsigned short> > genotype;

/// @deprecated Replaced by genepi::GenotypeParser
///
/// class to read genotypes from file and load into Individuals' arrays
class GenotypeLoader{
  
public:
  ///constructor
  GenotypeLoader();
  ///destructor
  virtual ~GenotypeLoader(){};
  ///read genotypesfile into a string array
  void Read(const char* filename, unsigned NumLociInLocusfile, bclib::LogWriter& Log);
  /**
     retrieve an Individual's genotype.
     reads genotypes as strings from string matrix (read from file) , converts them to vectors of unsigned short ints and allocates them to the Individual's genotype vector.
     \param i Individual number (count from 1)
     \param SexColumn index of sex column or 0 if there is none.
     \param Loci Genome object, to count through through loci
     \param genotypes pointer to vector to place the genotypes
     \param Missing indicators for missing genotypes, to be set while assigning genotypes
  */
  virtual void GetGenotype(int i, const Genome &Loci, 
		   std::vector<genotype>* genotypes, bool **Missing)const;

  ///check genotypes for unobserved alleles
  /// [DDF: apparently only called from hapmixmap]
  bool CheckForUnobservedAlleles(const bclib::DataMatrix& LocusData, bclib::LogWriter& Log);

  ///returns number of individuals
  virtual unsigned getNumberOfIndividuals()const;
  ///reurns number of loci (columns) in genotypesfile
  unsigned NumLoci()const;
  ///determines if an individual is female
  bool isFemale(unsigned i)const;
  ///determines if genotypesfile is in pedfile format
  bool isPedFile()const;
  ///returns the header of genotypesfile
  const std::vector<std::string>& getHeader()const;
  ///frees memory
  virtual void clear();
  ///returns index of sex column or 0 if there is none
  int getSexColumn()const;
protected:
  Matrix_s geneticData_;
  
  bool IsPedFile;
  int NumIndividuals, numDiploid, SexColumn;

  std::vector<unsigned short> GetGenotype(unsigned locus, int individual)const;
  void throwGenotypeError(int ind, int locus, std::string label, int g0, int g1, int numalleles)const;
/**
   check an Individual's genotypes are valid.
   writes error messages to cerr as logfile is not available yet.
   \param numhaploid number of haploid autosomal genotypes
   \param numdiploid number of diploid autosomal genotypes
   \param numhaploidX number of haploid X genotypes
   \param numdiploidX number of diploid X genotypes
   \param i individual number
   \param ID individual ID
*/
  void CheckGenotypes(unsigned long numObserved, unsigned long numhaploid, unsigned long numdiploid, 
		      unsigned long numhaploidX, unsigned long numdiploidX, unsigned i, const std::string& ID)const;

  std::vector<unsigned short> GetGenotype(const std::string & genostring)const;
  bool determineIfPedFile()const;
  ///looks for sex column 
  void DetermineSexColumn(unsigned NumLociInLocusfile, bclib::LogWriter& Log);

};



/** @} */



#endif
