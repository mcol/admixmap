/***
Formatting .h

Prototypes of various HAPMIXMAP data formatting functions

Copyright (c) 2007 David O'Donnell

*/
#ifndef FORMATTING_H
#define FORMATTING_H
#include <vector>
#include <string>
#include <map>

//class std::ofstream;
typedef struct
{
  std::string rsnumber;
  unsigned long position;
  std::pair<char, char> alleles;

  void print(std::ostream& os)const;
}LocusInfo;

class HapMapLegend
{

  public:
    HapMapLegend(std::ifstream& legendfile);
    ~HapMapLegend() {
      clear();
    }
    void clear() {
      RSmap.clear();
      LocusVector.clear();
    };
    const LocusInfo& operator[](unsigned long i) {
      return LocusVector[i];
    }
    const LocusInfo& operator[](const std::string& s) {
      return LocusVector[RSmap[s]];
    }
    const std::string& getRSNumber(unsigned long i) {
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

    ///expand limits by 100kb either side
    void OffsetLimits();
  private:
    std::map<std::string, unsigned long > RSmap;
    std::vector<LocusInfo> LocusVector;
    unsigned long first, last;

    HapMapLegend();
    HapMapLegend(const HapMapLegend&);
    const HapMapLegend& operator=(const HapMapLegend&);
};

//genotype conversion

///determines if a given string is one of a vector of strings
bool isListedSNP(const std::string HapMapSNP, const std::vector<std::string>SNPs);

///returns the complement of a given base (A/T, C/G)
char getComplement(char a);

///returns the complement of a pair of bases
std::pair<char, char> getComplement(std::pair<char, char> a);

///turns a genotype encoded as a string of 2 bases into a pair of bases
std::pair<char, char> GenotypeString2Pair(const std::string& g);

///Encodes a string containing a pair of bases as a genotype string
std::string getGenotype(const std::string& g, std::pair<char, char> a);

///Encodes a pair of bases as a genotype string
std::string getGenotype(const std::pair<char, char>& g, std::pair<char, char> a);

//Encodes genotypes given in infilename, sorts and writes to outfilename
unsigned EncodeGenotypes(HapMapLegend& Legend,
const char* infilename, const char* outfilename);
#endif
