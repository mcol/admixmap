// *-*-C++-*-*
/** 
 *   GenotypeLoader.h 
 *   class to load and assign genotypes
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef GENOTYPELOADER_H
#define GENOTYPELOADER_H
#include <fstream>
#include <vector>
#include <string>
#include "common.h"

class LogWriter;
class Genome;

typedef std::vector<std::vector<unsigned short> > genotype;

///class to read genotypes from file and load into Individuals' arrays
class GenotypeLoader{
  
public:
  ///constructor
  GenotypeLoader();
  ///destructor
  virtual ~GenotypeLoader(){};
  ///read genotypesfile into a sring array
  void Read(const char* filename, LogWriter& Log);
  /**
     retrieve an Individual's genotype.
     reads genotypes as strings from string matrix (read from file) , converts them to vectors of unsigned short ints and allocates them to the Individual's genotype vector.
     \param i Individual number (count from 1)
     \param SexColumn index of sex column or 0 if there is none.
     \param Loci Genome object, to count through through loci
     \param genotypes pointer to vector to place the genotypes
     \param Missing indicators for missing genotypes, to be set while assigning genotypes
  */
  virtual void GetGenotype(int i, int SexColumn, const Genome &Loci, 
		   std::vector<genotype>* genotypes, bool **Missing)const;

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

protected:
  Matrix_s geneticData_;
  
  bool IsPedFile;
  int NumIndividuals, numDiploid;

  std::vector<unsigned short> GetGenotype(unsigned locus, int individual, int SexColumn)const;
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

  std::vector<unsigned short> GetGenotype(const std::string genostring)const;
  bool determineIfPedFile()const;

};
#endif
