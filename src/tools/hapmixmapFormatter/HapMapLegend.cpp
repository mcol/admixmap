/**
   \file HapMapLegend.cpp
   Class to read and store a HapMap legend file.
   This file is part of FPHD

   This program is free software distributed WITHOUT ANY WARRANTY. 
   You can redistribute it and/or modify it under the terms of the GNU General Public License, 
   version 2 or later, as published by the Free Software Foundation. 
   See the file COPYING for details.

   Copyright (c) David O'Donnell 2007
*/
#include "HapMapLegend.h"

using namespace::std;

// void LocusInfo::print(ostream& os)const
// {
//   os << rsnumber << " " << position << " " 
//      << alleles.first << " " << alleles.second;
// }

HapMapLegend::HapMapLegend(const char* filename) {
  ifstream LegendFile(filename);
  if (!LegendFile.is_open()) {
    string errstring("Could not open legend file, ");
    errstring.append(filename); 
    throw errstring;
  }

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
  LegendFile.close();
  numSubChrs = 1;
}

void HapMapLegend::print(ostream& os) {
  for(vector<LocusInfo>::const_iterator i = LocusVector.begin(); i != LocusVector.end(); ++i) {
    i->print(os);
    os << '\n';
  }

}

void HapMapLegend::setLimits(unsigned long l0, unsigned long l1) {
  if(l0 <= l1 && l1 < LocusVector.size()) {
    first = l0;
    last = l1+1;
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
  while( (first >0) && (distance < offset) ) {
    --first;
    distance += LocusVector[first+1].position - LocusVector[first].position;
  }
  //offset lower limit
  distance = 0;
  while( (last < LocusVector.size()-1) && (distance < offset) ) {
    ++last;
    distance += LocusVector[last].position - LocusVector[last-1].position;
  }

}

void HapMapLegend::DetermineCutPoints(unsigned maxLoci, unsigned long minOverlap){
  /*
    divides chromosome as follows
       <-------   flanking region   ------------------>
   ... |  st   | ov |   st  | ov | .... //  ...   | st |...
       <-  c 1    ->        <-   c 3  ->
               <-     c 2      ->      
    st = straight (belongs to a single subchr)
    ov = overlap  (belongs to 2 subchrs)
  */
  //Segments.clear();
  if(maxLoci-1 > last-first)
    throw string("ERROR: locus limit exceeds length of chromosome");

  if(maxLoci == last-first){
    for(unsigned i = first; i < last; ++i)
      LocusVector[i].subchr.push_back(0);
    return;
  }

  pair<unsigned, unsigned> st, ov;
  st.first = st.second = first;
  ov.first = ov.second = first + maxLoci-1;
  SubChrSizes.push_back(0);

  numSubChrs = 0;
  while(ov.second < last){
    for(unsigned i = st.first; i <= ov.second; ++i)
      LocusVector[i].subchr.push_back(numSubChrs);
    *(SubChrSizes.end()-1) += ov.second + 1 - st.first;

    //set pos to position of last locus in overlap
    unsigned long pos = LocusVector[ov.second].position;
    unsigned end_st = ov.second;
    //determine max position of start of overlap to satisfy minOverlap
    const unsigned long target = pos - minOverlap;

    if(target < LocusVector[st.first].position){
      throw string("ERROR: overlap is too large for specified chromosome size");
    }

    do{
      //cout << "pos = " << pos << " b= " << end_st << " " << target << endl;
      LocusVector[end_st].subchr.push_back(numSubChrs+1);
      end_st--;
      pos = LocusVector[end_st].position;
    }while(pos > target && end_st > st.first);

    //set middle co-ordinates and add to vector
    st.second = end_st;
    ov.first = end_st+1;
    //Segments.push_back(st);
    //Segments.push_back(ov);
    //cout << " >> " <<  st.first << " " << st.second << " " << ov.first << " " << ov.second << endl;
    SubChrSizes.push_back(ov.second+1 - ov.first);

    //set start of next non-overlap 
    st.first = st.second = ov.second+1;
    //set end of next overlap
    ov.second = ov.first + maxLoci - 1;
    ov.first = ov.second;
    ++numSubChrs;
  }
  //add last segment
  for(unsigned i = st.first; i < last; ++i)
    LocusVector[i].subchr.push_back(numSubChrs);
  *(SubChrSizes.end()-1) += last - st.first;
  

  ++numSubChrs;
//   for(vector<LocusInfo>::const_iterator i = LocusVector.begin(); i!= LocusVector.end(); ++i){
//     cout << i->rsnumber << " ";
//     for(vector<unsigned>::const_iterator j = i->subchr.begin(); j != i->subchr.end(); ++j)
//       cout << *j << " ";
//     cout << endl;
//   }
//   cout << endl << endl;
}

//the following 3 functions not yet tested
void HapMapLegend::AdjustFlankingRegions(const std::vector<std::string>& TypedLoci, float FlankSize){
  AdjustSubChrUpperLimits(TypedLoci, FlankSize);
  AdjustSubChrLowerLimits(TypedLoci, FlankSize);

}
void HapMapLegend::AdjustSubChrUpperLimits(const std::vector<std::string>& TypedLoci, float FlankSize){
  for(vector<string>::const_iterator t = TypedLoci.begin(); t != TypedLoci.end(); ++t){
    if(LocusVector[RSmap[*t]].subchr.size() > 1){//if in an overlap
      //find last typed locus in the overlapping region
      do{
	++t;
      }while(LocusVector[RSmap[*t]].subchr.size() > 1); 
      --t;

      const LocusInfo& TypedLocus = LocusVector[RSmap[*t]];
      //find end of overlap
      vector<LocusInfo>::iterator j = LocusVector.begin() + RSmap[*t];
      while(j->subchr.size()>1){
	++j;
      }
      //if distance between current typed locus and end of flanking region is too short, extend the overlap
      while((j-1)->position - TypedLocus.position < FlankSize){
	j->subchr.push_back(j->subchr[0]-1);
	++j;
      }
    }
  }
}
void HapMapLegend::AdjustSubChrLowerLimits(const std::vector<std::string>& TypedLoci, float FlankSize){
  for(vector<string>::const_iterator t = TypedLoci.begin(); t != TypedLoci.end(); ++t){
    if(LocusVector[RSmap[*t]].subchr.size() > 1){//if in an overlap
      //find first locus in overlap
      do{
	--t;
      }while(LocusVector[RSmap[*t]].subchr.size() > 1); 
      ++t;

      const LocusInfo& TypedLocus = LocusVector[RSmap[*t]];
      //find start of overlap
      vector<LocusInfo>::iterator j = LocusVector.begin() + RSmap[*t];
      while(j->subchr.size()>1){
	--j;
      }
      //if distance between current typed locus and start of flanking region is too short, extend the overlap
      while((j+1)->position - TypedLocus.position < FlankSize){
	j->subchr.push_back(j->subchr[0]+1);
	--j;
      }

    }
  }
}

void HapMapLegend::clear(){
  RSmap.clear();
  LocusVector.clear();
}

bool HapMapLegend::isInHapMap(const std::string s)const{
  return (RSmap.find(s) != RSmap.end());
}
const LocusInfo& HapMapLegend::operator[](unsigned long i)const{
  return LocusVector[i];
}
const LocusInfo& HapMapLegend::operator[](const std::string& s){
  return LocusVector[RSmap[s]];
}
const std::string& HapMapLegend::getRSNumber(unsigned long i)const{
  return LocusVector[i].rsnumber;
}
unsigned HapMapLegend::getIndex(const std::string& s){
  return RSmap[s];
}
unsigned HapMapLegend::size()const{
  return LocusVector.size();
}
unsigned long HapMapLegend::getFirstIndex()const{
  return first;
}
unsigned long HapMapLegend::getLastIndex()const{
  return last;
}
const std::string& HapMapLegend::getFirst()const{
  return LocusVector[first].rsnumber;
}
const std::string& HapMapLegend::getLast()const{
  return LocusVector[last].rsnumber;
}
  
unsigned HapMapLegend::getNumSubChromosomes()const{
  return numSubChrs;
}
const std::vector<unsigned>& HapMapLegend::getSubChromosomeSizes()const{
  return SubChrSizes;
}
