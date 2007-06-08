// *-*-C++-*-*
/**
   \file UserGenotypes.h
   class to read user genotypes and write in HAPMIXMAP-format
   This file is part of FPHD

   This program is free software distributed WITHOUT ANY WARRANTY. 
   You can redistribute it and/or modify it under the terms of the GNU General Public License, 
   version 2 or later, as published by the Free Software Foundation. 
   See the file COPYING for details.

   Copyright (c) David O'Donnell 2007
*/

#ifndef USERGENOTYPES_H
#define USERGENOTYPES_H
#include <fstream>
#include <string>
#include <vector>

class HapMapLegend;

///class to read user genotypes and write in HAPMIXMAP-format
class UserGenotypes{
public:
  ///Constructor
  UserGenotypes(HapMapLegend& Legend,const char* infilename);
  ///Destructor
  ~UserGenotypes();
  unsigned getNumberOfTypedLoci()const;
  const std::vector<std::string>& getTypedLoci()const;
  ///returns first typed locus by map position
  const std::string& getFirstTypedLocus()const;
  ///returns last typed locus by map position
  const std::string& getLastTypedLocus()const;
  /**
     Encodes genotypes given in infilename, sorts and writes to outfilename
     \param Legend HapMapLegend object
     \param outfileprefix output filename, without extension
     \param Missing character denoting missing genotypes
     \return number of individuals in file
  */
  unsigned FormatUserGenotypes(HapMapLegend& Legend, const char* outfileprefix, char MISSING);

private:
  ///input genotypesfile
  std::ifstream RawGenofile;
  std::string inFileName;
  ///ID label from infile header
  std::string ID;
  ///vector of names of typed loci, including those not in HapMap
  std::vector<std::string> TypedLoci;
  ///number of typed loci, excluding those not in HapMap
  unsigned NumTypedHapMapLoci;
  ///ranks
  std::vector<unsigned> ranks;

  UserGenotypes();
  UserGenotypes(const UserGenotypes&);
  ///rank loci by map position
  void RankLoci(const std::vector<unsigned long>& TypedPos);
  ///write header to output file(s)
  void WriteHeader(  std::vector<std::ofstream*>& OutGeno, HapMapLegend& Legend );
  unsigned EncodeGenotypes(std::vector<std::ofstream*>& OutGeno, HapMapLegend& Legend, char missing);


};
#endif
