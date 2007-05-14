// *-*-C++-*-*
/**
   HapMapLegend.h
   Class to read and store a HapMap legendfile
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

  void print(std::ostream& os)const;
}LocusInfo;

///class to represent a HapMap legend
class HapMapLegend{
  
public:
  HapMapLegend(std::ifstream& legendfile);
  ~HapMapLegend(){
    clear();
  }
  void clear(){
    RSmap.clear();
    LocusVector.clear();
  };
  bool isInHapMap(const std::string s)const{
    return (RSmap.find(s) != RSmap.end());
  }
  const LocusInfo& operator[](unsigned long i){
    return LocusVector[i];
  }
  const LocusInfo& operator[](const std::string& s){
    return LocusVector[RSmap[s]];
  }
  const std::string& getRSNumber(unsigned long i){
    return  LocusVector[i].rsnumber;
  }
  void print(std::ostream& os);
  unsigned size()const{return LocusVector.size();};
  
  void setLimits(unsigned long l0, unsigned long l1);
  void setLimits(const std::string& l0, const std::string& l1);
  
  unsigned long getFirstIndex()const{return first;}
  unsigned long getLastIndex()const{return last;}
  const std::string& getFirst()const{return LocusVector[first].rsnumber;}
  const std::string& getLast()const{return LocusVector[last].rsnumber;}
  
  ///expand limits by offset either side
    void OffsetLimits(float offset);
private:
  std::map<std::string, unsigned long > RSmap;
  std::vector<LocusInfo> LocusVector;
  unsigned long first, last;
  
  HapMapLegend();
  HapMapLegend(const HapMapLegend&);
  const HapMapLegend& operator=(const HapMapLegend&);
};

#endif
