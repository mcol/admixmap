/***
GenotypeEncoder.cpp
Functions to turn genotypes coded as pairs of alleles to strings
eg A/C ->"1,2"

Copyright (c) David O'Donnell 2007
*/

#include "Formatting.h"
#include <fstream>
#include <sstream>
#include "ranker.h"
#include "gsl/gsl_permutation.h"
#include "StringSplitter.h"

using namespace::std;

#define MISSING 'N'

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
  s << "\"" << allele1 << "," << allele2 << "\"";
  return s.str();
}

///turns a genotype encoded as a string of 2 bases into a pair of bases
pair<char, char> GenotypeString2Pair(const string& g) {
  pair<char, char> gPair;
  
  //TODO: dequote, requires StringConvertor
  if(!g.length())throw string("Error in GenotypeString2Pair: empty string!\n");
  
  gPair.first = g[0];
  switch(g.length()) {
    case 1:{                                       //haploid genotype, not fully supported yet
      gPair.second = MISSING;                     //use missing-value character as second base
      break;
    }
    case 2:{                                       //no separator
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
  ofstream OutGeno(outfilename);

  vector<string> rslist;
  string line, ID;
  //Read Header of raw genotypes file
  RawGeno >> ID;                                  //skip first col in header
  OutGeno << ID;
  getline(RawGeno, line);
  StringSplitter::Tokenize(line, rslist, " \t");

  const unsigned NumTypedLoci = rslist.size();

  //create vector of positions of typed loci
  vector<unsigned long> TypedPos;
  for(vector<string>::const_iterator i = rslist.begin(); i!= rslist.end(); ++i) {
    TypedPos.push_back(Legend[*i].position);
  }

  //sort typed loci by position
  const unsigned num = TypedPos.size();
  vector<unsigned> ranks;
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

  //copy(ranks.begin(), ranks.end(), ostream_iterator<unsigned>(cout, " "));

  //for(vector<unsigned>::const_iterator i = ranks.begin(); i !=ranks.end(); ++i){
  //cout << rslist[*i] << " " << TypedPos[*i] << " ";
  //Legend[rslist[*i]].print(cout);
  //}

  //set first and last loci in Legend
  Legend.setLimits(rslist[ranks[0]], rslist[ranks[ranks.size()-1]]);

  //write genotypes in order
  vector<string> genotypesline(NumTypedLoci);
  //finish writing header
  for(vector<unsigned>::const_iterator i = ranks.begin(); i !=ranks.end(); ++i) {
    OutGeno << " " <<  rslist[*i];
  }
  OutGeno << endl;

  RawGeno >> ID;                                  //read ID
  for(vector<string>::iterator i = genotypesline.begin(); i != genotypesline.end(); ++i)
    RawGeno >> *i;
  //TODO: use copy to read genotypes
  //TODO: check correct number of genotypes in each row
  unsigned NumInd = 0;
  while(!RawGeno.eof()) {
    ++NumInd;
    OutGeno << ID;                                // << " ";
    //copy(genotypesline.begin(), genotypesline.end(), ostream_iterator<string>(OutGeno, " "));
    for(vector<unsigned>::const_iterator i = ranks.begin(); i !=ranks.end(); ++i) {
      OutGeno << " " <<  getGenotype(genotypesline[*i], Legend[rslist[*i]].alleles);
    }

    OutGeno << endl;
    RawGeno >> ID;                                //read ID
    for(vector<string>::iterator i = genotypesline.begin(); i != genotypesline.end(); ++i)
      RawGeno >> *i;
  }
  RawGeno.close();
  OutGeno.close();
  return NumInd;
}

void LocusInfo::print(ostream& os)const
{
  os << rsnumber << " " << position << " " << alleles.first << " " << alleles.second << endl;
}

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

void HapMapLegend::OffsetLimits() {
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
