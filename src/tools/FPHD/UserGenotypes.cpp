/**
   \file UserGenotypes.cpp
   class to read user genotypes and write in HAPMIXMAP-format
   This file is part of FPHD

   This program is free software distributed WITHOUT ANY WARRANTY. 
   You can redistribute it and/or modify it under the terms of the GNU General Public License, 
   version 2 or later, as published by the Free Software Foundation. 
   See the file COPYING for details.

   Copyright (c) David O'Donnell 2007
*/

#include "UserGenotypes.h"
#include "GenotypeEncoding.h"
#include "HapMapLegend.h"
#include "bclib/ranker.h"
#include "gsl/gsl_permutation.h"
#include "bclib/StringSplitter.h"
#include <iostream>
#include <vector>
#include <sstream>
#include <string.h>

using namespace std;
#define REMOTE_POSITION 999999999

UserGenotypes::UserGenotypes(HapMapLegend& Legend, const char* infilename){
  if(infilename == 0 || strlen(infilename)==0)
    return;

  RawGenofile.open(infilename);
  if(!RawGenofile.is_open()){
    cerr << "** ERROR: could not open " << infilename << endl;
    exit(1);
  }
  inFileName = infilename;
  vector<string> NonHapMapLoci;//typed loci not in HapMap (usually because they are monomorphic in the HapMap sample)
  string line;

  //Read Header of raw genotypes file
  RawGenofile >> ID;

  getline(RawGenofile, line);
  if(line.find_first_not_of(" \t") == string::npos){
    cerr << "** ERROR: " << infilename << " header is empty!" << endl;
    exit(1);
  }
  vector<string> tokens;//vector of typed loci
  bclib::StringSplitter::Tokenize(line, TypedLoci, " \t");

   //create vector of positions of typed loci
  vector<unsigned long> TypedPos;
  for(vector<string>::const_iterator i = TypedLoci.begin(); i!= TypedLoci.end(); ++i) {
    //check the locus in is the HapMap legend file
    if(!Legend.isInHapMap(*i)){
      NonHapMapLoci.push_back(*i);
      //assign a position off end of genome to ensure a rank higher than any locus in HapMap
      TypedPos.push_back(REMOTE_POSITION);
    }
    else{
      //assign position from legend file
      TypedPos.push_back(Legend[*i].position);
    }
  }
  NumTypedHapMapLoci = TypedLoci.size() - NonHapMapLoci.size();
  if( NumTypedHapMapLoci == 0){
    throw string("ERROR: no HapMap loci found in" + inFileName );
  }
  //report loci not in HapMap
  else if(NonHapMapLoci.size()){
    cerr << "** Warning: " << NonHapMapLoci.size() 
	 << " loci are not in the HapMap for this chromosome and will be omitted";
    //     if(NonHapMapLoci.size()<=10){
    //       cerr << ":" << endl;
    //       copy(NonHapMapLoci.begin(), NonHapMapLoci.end(), ostream_iterator<string>(cerr, " "));
    //     }
    cerr << endl;
  }

  RankLoci(TypedPos);
  TypedPos.clear();
}
UserGenotypes::~UserGenotypes(){
  RawGenofile.close();
}


void UserGenotypes::RankLoci(const vector<unsigned long>& TypedPos){
  const unsigned Num = TypedPos.size();
  //sort typed loci by position
  ranks.assign(Num, 0);
  ::rank(TypedPos, ranks, "default");

  gsl_permutation * p = gsl_permutation_alloc (Num);
  for (uint i = 0; i < Num; ++i) {
    p->data[i] = (int)ranks[i]-1;
  }
  //get inverse of permuation
  //if(!gsl_permutation_valid (p))exit(2);
  gsl_permutation * q = gsl_permutation_alloc (Num);
  gsl_permutation_inverse (q, p);
  //gsl_permutation_fprintf (stdout, q, " %u");

  //clear ranks and start again, this time exclusing loci not in HapMap
  ranks.clear();
  for (uint i = 0; i < NumTypedHapMapLoci; ++i) {
    //ranks[i] = index of ith locus by map position (if TypedLoci were sorted)
    ranks.push_back(q->data[i]);
  }

  //tidy up
  gsl_permutation_free(p);
  gsl_permutation_free(q);
}

unsigned UserGenotypes::getNumberOfTypedLoci()const{
  return NumTypedHapMapLoci;
}
const string& UserGenotypes::getFirstTypedLocus()const{
  return TypedLoci[ranks[0]];
}
const string& UserGenotypes::getLastTypedLocus()const{
  return TypedLoci[ranks[ranks.size()-1]];
}
const vector<string>& UserGenotypes::getTypedLoci()const{
  return TypedLoci;
}
unsigned UserGenotypes::FormatUserGenotypes(HapMapLegend& Legend, const char* outfileprefix, char Missing){
  //open output file(s)
  vector<ofstream*> outgeno;
  for(unsigned j = 0; j < Legend.getNumSubChromosomes(); ++j){
    stringstream filename;
    if(Legend.getNumSubChromosomes()==1)
      filename << outfileprefix << ".txt";
    else
      filename << outfileprefix << "_" << j+1 << ".txt";
    
    outgeno.push_back(new ofstream(filename.str().c_str()));
    if (!outgeno[j]->is_open()) {
      cerr << "** ERROR: could not open " << filename << endl;
      exit(1);
    }
  }

  WriteHeader(outgeno, Legend);
  unsigned NumTypedIndivs = EncodeGenotypes(outgeno, Legend, Missing);

  //clean up
  for(vector<ofstream*>::iterator i = outgeno.begin(); i != outgeno.end(); ++i){
    (*i)->close();
    delete *i;
  }
  outgeno.clear();

  return NumTypedIndivs;
}

void UserGenotypes::WriteHeader(  vector<ofstream*>& OutGeno, HapMapLegend& Legend ){
  //start with ID
  for(unsigned j = 0; j < Legend.getNumSubChromosomes(); ++j)
    *(OutGeno[j]) << ID;

  //write locus names
  for(vector<unsigned>::const_iterator i = ranks.begin(); i !=ranks.end(); ++i) {
    if(Legend.isInHapMap(TypedLoci[*i])){
      const vector<unsigned>& subchrs = Legend[TypedLoci[*i]].subchr;
      for(vector<unsigned>::const_iterator j = subchrs.begin(); 
	  j != subchrs.end(); ++j){
	*(OutGeno[*j]) << " " <<  TypedLoci[*i];
      }
    }
  }
  //finish with newline
  for(unsigned j = 0; j < Legend.getNumSubChromosomes(); ++j)
    *(OutGeno[j]) << endl;
}

unsigned UserGenotypes::EncodeGenotypes(vector<ofstream*>& OutGeno, HapMapLegend& Legend, char MISSING) {
  string IndivID, line;
  unsigned NumInd = 0;
  vector<pair<string, unsigned> > AllMissingIndivs;
  string ValidChars = "ACTG:,; \t";
  ValidChars.append(1, MISSING);

  //read ID of 1st individual
  RawGenofile >> IndivID;
  while(getline(RawGenofile, line)) {
    ++NumInd;

    //check for invalid characters 
    string::size_type badchar = line.find_first_not_of(ValidChars);
    if( badchar != string::npos){
      cerr << "** ERROR: invalid character on line " << NumInd +1 << " of " << inFileName << ": " << line[badchar] << endl;
      exit(1);
    }

    //check there is at least one observed genotype
    if(line.find_first_of("ACTG") == string::npos){
      AllMissingIndivs.push_back(pair<string, unsigned>(IndivID, NumInd));
      RawGenofile >> IndivID;
      --NumInd;   
      continue;
    }

    //Tokenize line into strings
    vector<string> genotypesline;
    bclib::StringSplitter::Tokenize(line, genotypesline, " \t");
    //check individual has the correct number of genotypes
    if(genotypesline.size() != TypedLoci.size()){
      cerr << "** ERROR: inconsistent row lengths in " << inFileName << ". " 
	   << TypedLoci.size() << " loci in header but " << genotypesline.size() << " in row " << NumInd;
      exit(1); 
    }

    //write ID
    for(unsigned j = 0; j < Legend.getNumSubChromosomes(); ++j)
      *(OutGeno[j]) << IndivID;

    //convert genotypes and write to file
    for(vector<unsigned>::const_iterator i = ranks.begin(); i !=ranks.end(); ++i) {
      if(Legend.isInHapMap(TypedLoci[*i])){
	const vector<unsigned>& subchrs = Legend[TypedLoci[*i]].subchr;
	for(vector<unsigned>::const_iterator j = subchrs.begin(); 
	    j != subchrs.end(); ++j){
	  *(OutGeno[*j]) << " " <<  getGenotype(genotypesline[*i], Legend[TypedLoci[*i]].alleles, MISSING);
	}
      }
    }
   
    //finish with newline
    for(unsigned j = 0; j < Legend.getNumSubChromosomes(); ++j)
      *(OutGeno[j]) << endl;
    //read ID
    RawGenofile >> IndivID;
  }
  RawGenofile.close();

  //warn user of individuals with all missing genotypes
  if(AllMissingIndivs.size()){
    cerr << "** Warning: The following individuals have no observed genotypes and have been omitted." << endl;
    for(vector<pair<string, unsigned> >::const_iterator i = AllMissingIndivs.begin(); i != AllMissingIndivs.end(); ++i){
      cout << i->first << "(" << i->second << ")" << endl;
    }
    cout << endl;
  }
  return NumInd;
}



