/**
   GenotypeEncoding.cpp
   Functions to turn genotypes coded as pairs of alleles to strings
   eg A/C ->"1,2"
   This file is part of FPHD
   
   This program is free software distributed WITHOUT ANY WARRANTY. 
   You can redistribute it and/or modify it under the terms of the GNU General Public License, 
   version 2 or later, as published by the Free Software Foundation. 
   See the file COPYING for details.
   
   Copyright (c) David O'Donnell 2007
*/

#include "GenotypeEncoding.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "bcppcl/ranker.h"
#include "gsl/gsl_permutation.h"
#include "bcppcl/StringSplitter.h"

using namespace::std;

#define MISSING 'N'
#define NA_HAPLOID 'X'

//small function to determine if a given string is in a list of strings
bool isListedSNP(const string HapMapSNP, const vector<string>SNPs) {
  bool b = false;
  if(find(SNPs.begin(), SNPs.end(), HapMapSNP) < SNPs.end())b = true;
  return b;
}

char getComplement(char a){
  char b;
  switch(a){
  case 'A':{
    b = 'T';
    break;
  }
  case 'T':{
    b = 'A';
    break;
  }
  case 'C':{
    b = 'G';
    break;
  }
  case 'G':{
    b = 'C';
    break;
  }
  default:{
    //cerr << "Error: bases must be A, C, T or G" << endl;
    //exit(1);
    throw string("Error: bases must be A, C, T or G");
  }

  }
  return b;
}

pair<char, char> getComplement(pair<char, char> a) {
  pair<char, char> b;
  b.first = getComplement(a.first);
  b.second = getComplement(a.second);
  return b;
}

///converts a base pair string g, eg "CT" or "A/G", to an admixmap/hapmixmap format genotype eg "1,2"
///a is the reference string of alleles eg "C/T"
string getGenotype(const string& g, pair<char, char> a) {

  return getGenotype(GenotypeString2Pair(g), a);

}

///converts a base pair to an admixmap/hapmixmap format genotype eg "1,2"
///a is the reference string of alleles eg "C/T"
string getGenotype(const pair<char, char>& g, pair<char, char> a) {
  int allele1 = (int)(g.first  == a.first)+ 2*(int)(g.first  == a.second);
  int allele2 = (int)(g.second == a.first)+ 2*(int)(g.second == a.second);

  //check for reverse strand coding
  //NB: assumes each base is not paired with its complement ie a is never "AT" or "CG"
  if((allele1 == 0)  && (g.first != MISSING)) {
    //throw(1);
    pair< char, char> b = getComplement(a);
    allele1 = (int)(g.first  == b.first)+ 2*(int)(g.first  == b.second);
    allele2 = (int)(g.second == b.first)+ 2*(int)(g.second == b.second);
  }

  stringstream s;
  s << "\"" << allele1;
  //write second allele if diploid
  //TODO: check genotypes are consistent, ie all diploid or all haploid
  if( g.second != NA_HAPLOID) 
    s << "," << allele2;
  s << "\"";
  return s.str();
}

///turns a genotype encoded as a string of 2 bases into a pair of bases
///assumes genotype has form XY, X:Y, X;Y, or X,Y (no quotes)
///or X if haploid. X and Y may be A,C,T,G or N.
pair<char, char> GenotypeString2Pair(const string& g) {
  pair<char, char> gPair;
  
  if(!g.length())throw string("Error in GenotypeString2Pair: empty string!\n");
  
  gPair.first = g[0];
  switch(g.length()) {
    case 1:{                                       //haploid genotype
      gPair.second = NA_HAPLOID;                   //use special character as second base
      break;
    }
    case 2:{
      if(g.find_first_of(":;,")==1)                        //separator but no second base
	gPair.second = NA_HAPLOID;
      else                                       //no separator
	gPair.second = g[1];
      break;
    }
    case 3:{                                       //separator used
      gPair.second = g[2];
      break;
    }
    default:{                                      //too long
      throw string("Error in GenotypeString2Pair: string is too long!");
    }
  }
  
  return gPair;
}

unsigned EncodeGenotypes(HapMapLegend& Legend, const char* infilename, const char* outfilename) {
  ifstream RawGeno(infilename);
  if(!RawGeno.is_open()){
    cerr << "** ERROR: could not open " << infilename << endl;
    exit(1);
  }
  ofstream OutGeno(outfilename);
  if(!OutGeno.is_open()){
    cerr << "** ERROR: could not open " << outfilename << endl;
    exit(1);
  }

  vector<string> TypedLoci;
  string line, ID;
  //Read Header of raw genotypes file
  RawGeno >> ID;                                  //skip first col in header
  OutGeno << ID;
  getline(RawGeno, line);
  if(line.find_first_not_of(" \t") == string::npos){
    cerr << "** ERROR: " << infilename << " header is empty!" << endl;
    exit(1);
  }
  StringSplitter::Tokenize(line, TypedLoci, " \t");

  const unsigned NumTypedLoci = TypedLoci.size();

  //create vector of positions of typed loci
  vector<unsigned long> TypedPos;
  for(vector<string>::const_iterator i = TypedLoci.begin(); i!= TypedLoci.end(); ++i) {
    //check the locus in is the HapMap legend file
    if(!Legend.isInHapMap(*i)){
      cerr << "** ERROR: " << *i << " is not the name of a HapMap locus on this chromosome" << endl;
      exit(1);
    }
      TypedPos.push_back(Legend[*i].position);
  }

  //sort typed loci by position
  const unsigned num = TypedPos.size();
  vector<unsigned> ranks(num);
  rank(TypedPos, ranks, "default");

  gsl_permutation * p = gsl_permutation_alloc (num);
  for (uint i = 0; i < num; ++i) {
    p->data[i] = (int)ranks[i]-1;
  }
  //get inverse of permuation
  //if(!gsl_permutation_valid (p))exit(2);
  gsl_permutation * q = gsl_permutation_alloc (num);
  gsl_permutation_inverse (q, p);
  //gsl_permutation_fprintf (stdout, q, " %u");
  for (uint i = 0; i < num; ++i) {
    ranks[i] = q->data[i];
  }

  //tidy up
  TypedPos.clear();
  gsl_permutation_free(p);
  gsl_permutation_free(q);

  //set first and last loci in Legend
  Legend.setLimits(TypedLoci[ranks[0]], TypedLoci[ranks[ranks.size()-1]]);

  //finish writing header
  for(vector<unsigned>::const_iterator i = ranks.begin(); i !=ranks.end(); ++i) {
    //if(Legend.isInHapMap(TypedLoci[*i]))
      OutGeno << " " <<  TypedLoci[*i];
  }
  OutGeno << endl;

  unsigned NumInd = 0;
  vector<pair<string, unsigned> > AllMissingIndivs;
  //read ID of 1st individual
  RawGeno >> ID;
  while(getline(RawGeno, line)) {
    ++NumInd;

    //check for invalid characters 
    string::size_type badchar = line.find_first_not_of("ACTGN:,; \t");
    if( badchar != string::npos){
      cerr << "** ERROR: invalid character on line " << NumInd +1 << " of " << infilename << ": " << line[badchar] << endl;
      exit(1);
    }

    //check there is at least one observed genotype
    if(line.find_first_of("ACTG") == string::npos){
      AllMissingIndivs.push_back(pair<string, unsigned>(ID, NumInd));
      RawGeno >> ID;
      --NumInd;   
      continue;
    }

    //Tokenize line into strings
    vector<string> genotypesline;
    StringSplitter::Tokenize(line, genotypesline, " \t");
    //check individual has the correct number of genotypes
    if(genotypesline.size() != NumTypedLoci){
      cerr << "** ERROR: inconsistent row lengths in " << infilename << ". " 
	   << NumTypedLoci << " loci in header but " << genotypesline.size() << " in row " << NumInd;
      exit(1); 
    }

    OutGeno << ID;
    //convert genotypes and write to file
    for(vector<unsigned>::const_iterator i = ranks.begin(); i !=ranks.end(); ++i) {
      OutGeno << " " <<  getGenotype(genotypesline[*i], Legend[TypedLoci[*i]].alleles);
    }

    OutGeno << endl;
    RawGeno >> ID;                                //read ID
  }
  RawGeno.close();
  OutGeno.close();

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


