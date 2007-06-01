// *-*-C++-*-*
/**
   \file HapMapLegend.h
   Class to read and store a HapMap legendfile.
   This file is part of FPHD

   This program is free software distributed WITHOUT ANY WARRANTY. 
   You can redistribute it and/or modify it under the terms of the GNU General Public License, 
   version 2 or later, as published by the Free Software Foundation. 
   See the file COPYING for details.

   Copyright (c) David O'Donnell 2007
*/
#ifndef HAPMAPLEGEND_H
#define HAPMAPLEGEND_H

#include <vector>
#include <string>
#include <map>

//class std::ofstream;

///struct to hold locus information
typedef struct{
  std::string rsnumber;
  unsigned long position;
  std::pair<char, char> alleles;
  std::vector<unsigned> subchr;

  void print(std::ostream& os)const;
}LocusInfo;

///class to represent a HapMap legend
class HapMapLegend{
  
public:
  /**
     Constructor.
     \param filename name of file containing HapMap legend.
  */
  HapMapLegend(const char* filename);
  ///Destructor
  ~HapMapLegend(){
    clear();
  }
  ///free memory
  void clear();

  ///determines if a locus is listed in the legend
  bool isInHapMap(const std::string s)const;

  ///access locus info by index
  const LocusInfo& operator[](unsigned long i)const;

  ///access locus info by rsnumber
  const LocusInfo& operator[](const std::string& s);

  ///get the rs number of ith locus
  const std::string& getRSNumber(unsigned long i)const;

  ///print contents to stream
  void print(std::ostream& os);

  ///returns the number of loci in the legend
  unsigned size()const;
  
  ///set chromosome limits at indices lo and l1
  void setLimits(unsigned long l0, unsigned long l1);

  ///set chromosome limits at indices of loci with names l0 and l1
  void setLimits(const std::string& l0, const std::string& l1);
  
  ///get index of first locus
  unsigned long getFirstIndex()const;
  ///get index of last locus
  unsigned long getLastIndex()const;
  ///get name of first locus
  const std::string& getFirst()const;
  ///get name of last locus
  const std::string& getLast()const;
  
  ///expand limits by offset either side
  void OffsetLimits(float offset);

  /**
     Determine points to divide up the chromosome.
     \param maxLoci maximum number of loci per sub-chromosome
     \param minOverlap minimum length in bp of overlap between subchromosomes
  */
  void DetermineCutPoints(unsigned maxLoci, unsigned long minOverlap);

  ///return number of sub-chromosomes
  unsigned getNumSubChromosomes()const;

  ///return sizes of subchromosomes
  const std::vector<unsigned>& getSubChromosomeSizes()const;

  /**
     adjust flanking regions on each sub-chromosome according to any typed loci.
     NOT TESTED.
     \param TypedLoci vector of names of typed loci
     \param FlankSize desired flanking region size in bp
  */
  void AdjustFlankingRegions(const std::vector<std::string>& TypedLoci, float FlankSize);

private:
  ///map rsnumber to index
  std::map<std::string, unsigned long > RSmap;
  ///vector of locus information
  std::vector<LocusInfo> LocusVector;
  unsigned long first, last;
  unsigned numSubChrs;
  std::vector<unsigned> SubChrSizes;

  HapMapLegend();
  HapMapLegend(const HapMapLegend&);
  const HapMapLegend& operator=(const HapMapLegend&);

  void AdjustSubChrUpperLimits(const std::vector<std::string>& TypedLoci, float FlankSize);
  void AdjustSubChrLowerLimits(const std::vector<std::string>& TypedLoci, float FlankSize);
};

#endif
