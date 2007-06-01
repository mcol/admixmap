/**
   \file GenotypeEncoding.cpp
   Functions to turn genotypes coded as pairs of alleles to strings.
   eg A/C ->"1,2"
   This file is part of FPHD
   
   This program is free software distributed WITHOUT ANY WARRANTY. 
   You can redistribute it and/or modify it under the terms of the GNU General Public License, 
   version 2 or later, as published by the Free Software Foundation. 
   See the file COPYING for details.
   
   Copyright (c) David O'Donnell 2007
*/

#include "GenotypeEncoding.h"
#include <sstream>

using namespace::std;

//#define MISSING 'N'
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
string getGenotype(const string& g, pair<char, char> a, char MISSING) {

  return getGenotype(GenotypeString2Pair(g), a, MISSING);

}

///converts a base pair to an admixmap/hapmixmap format genotype eg "1,2"
///a is the reference string of alleles eg "C/T"
string getGenotype(const pair<char, char>& g, pair<char, char> a, char MISSING) {
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


