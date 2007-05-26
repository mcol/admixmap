/**
   HapMapLegend.cpp
   This file is part of FPHD
   Class to read and store a HapMap legend file.

   This program is free software distributed WITHOUT ANY WARRANTY. 
   You can redistribute it and/or modify it under the terms of the GNU General Public License, 
   version 2 or later, as published by the Free Software Foundation. 
   See the file COPYING for details.

   Copyright (c) David O'Donnell 2007
*/
#include "HapMapLegend.h"
#include <fstream>

using namespace::std;

HapMapLegend::HapMapLegend(ifstream& LegendFile) {
  string header;
  LocusInfo info;
  getline(LegendFile, header);                    //skip header
  LegendFile >> info.rsnumber >> info.position >> info.alleles.first >> info.alleles.second;

  unsigned index = 0;
  first = 0;
  last  = 0;
  while(!LegendFile.eof()) {
    LocusVector.push_back(info);
    RSmap[info.rsnumber] = index++;
    ++last;
    LegendFile >> info.rsnumber >> info.position >> info.alleles.first >> info.alleles.second;
  }
  if(last == 0 ){
    throw "Empty legend file";
  }
}

void HapMapLegend::print(ostream& os) {
  for(vector<LocusInfo>::const_iterator i = LocusVector.begin(); i != LocusVector.end(); ++i) {
    i->print(os);
  }

}

void HapMapLegend::setLimits(unsigned long l0, unsigned long l1) {
  if(l0 <= l1 && l1 < LocusVector.size()) {
    first = l0;
    last = l1;
  }
  else throw string("Error in HapMapLegend::setLimits: invalid limits");
}

void HapMapLegend::setLimits(const std::string& l0, const std::string& l1) {
  setLimits(RSmap[l0], RSmap[l1]);
}

void HapMapLegend::OffsetLimits(float offset) {
  //offset lower limit
  //vector<LocusInfo>::const_iterator i = LocusInfo.begin()+first;
  unsigned long distance = 0;
  while( (first >0) && (distance < 1e5) ) {
    --first;
    distance += LocusVector[first+1].position - LocusVector[first].position;
  }
  //offset lower limit
  distance = 0;
  while( (last < LocusVector.size()-1) && (distance < 1e5) ) {
    ++last;
    distance += LocusVector[last].position - LocusVector[last-1].position;
  }

}

void LocusInfo::print(ostream& os)const
{
  os << rsnumber << " " << position << " " << alleles.first << " " << alleles.second << endl;
}
