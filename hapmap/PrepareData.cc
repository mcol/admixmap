/**
 C++ Program to convert raw ColonCancer data to HAPMIXMAP format.
 Program is very ad hoc, according to the format of the raw data (see README file)
 but could be modified to be distributed.

writes 2 genotypesfiles: one with genotypes from the 300 array, one with only those from the 240 array.
Author: David O'Donnell, 2006  

TODO: write outcome file
*/
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
//#include <map>

using namespace::std;

//small function to determine if a given string is in a list of strings
bool isListedSNP(const string HapMapSNP, const vector<string>SNPs){
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
    cerr << "Error: bases must be A, C, T or G" << endl;
    exit(1);
  }

  }
  return b;
}

pair<char, char> getComplement(pair<char, char> a){
  pair<char, char> b;
  b.first = getComplement(a.first);
  b.second = getComplement(a.second);
  return b;
}

string getGenotype(string g, pair<char, char> a){
  //converts a base pair string g, eg "CT", to an admixmap format genotype eg "1,2"
  //a is the reference string of alleles eg "C/T"
  int allele1 = (int)(g[0]==a.first)+2*(int)(g[0]==a.second);
  int allele2 = (int)(g[1]==a.first)+2*(int)(g[1]==a.second);

  //check for reverse strand coding
  //NB: assumes each base is not paired with its complement ie a is never "AT" or "CG" 
  if((allele1 == 0)  && (g[0] != 'N')){
    //throw(1);
    pair< char, char> b = getComplement(a);
    int allele1 = (int)(g[0]==b.first)+2*(int)(g[0]==b.second);
    int allele2 = (int)(g[1]==b.first)+2*(int)(g[1]==b.second);
  }

  stringstream s;
  s << "\"" << allele1 << "," << allele2 << "\"";
  return s.str();
}

int main(){
  ifstream rawgenofile("Paul_550_22.txt");
  ifstream locusdatafile("Rs_550_22.txt");
  ifstream legendfile("../chr22data/chr22_legend.txt");

  ofstream genotypes240file("genotypes240.txt"); 
  ofstream genotypes300file("genotypes300.txt"); 

  vector<string> loci;
  vector<bool> includeLocus;
  vector<bool> Study240;
  vector< pair<char, char> > legend;

  string s, scrap;
  pair<char , char>alleles;

  vector<string> HapMapLoci;
  getline(legendfile, scrap);//skip header
  legendfile >> s;
  while(!legendfile.eof()){
    HapMapLoci.push_back(s);
    getline(legendfile, scrap);
    legendfile >> s;
  }
  legendfile.close();


  //read in list of typed loci, excluding monomorphic loci
  cout << "Reading typed loci ..." << endl;
  unsigned i1, i2, position, index;
  locusdatafile >> i1 >> i2 >> position >> s >> index;
  while(!locusdatafile.eof()){
    if(isListedSNP(s, HapMapLoci)){//exclude loci not in legend file
      loci.push_back(s);
      includeLocus.push_back(true);
      Study240.push_back((bool)(index==2));
    }
    else includeLocus.push_back(false);
    locusdatafile >> i1 >> i2 >> position >> s >> index;
  }
  locusdatafile.close();
  HapMapLoci.clear();
  const unsigned NLoci = includeLocus.size();
  cout << "Found " << NLoci << " loci, Using " << loci.size() << " loci" << endl;

  //read in legend
  legendfile.clear();
  legendfile.open("chr22_legend.txt");
  cout << "Reading in allele codes ..." << endl;
  getline(legendfile, scrap);//skip header
  legendfile >> s >> scrap >> alleles.first >> alleles.second;
  while(!legendfile.eof()){
    if(isListedSNP(s, loci)){//skip untyped loci
      legend.push_back(alleles);
    }
    legendfile >> s >> scrap >> alleles.first >> alleles.second;
  }
  legendfile.close();


  //write table of loci and bases to file
  ofstream locustable("LocusTable.txt");
  locustable << "Locus\tallele1\tallele2" << endl;
  for(  unsigned counter = 0;counter < loci.size(); ++counter){
    locustable << loci[counter] << "\t" << legend[counter].first << " " << legend[counter].second << endl;
    //     cout << loci[counter] << " "  << legend[counter].first << " " << legend[counter].second << endl;
  }
  locustable.close();

//convert raw genotypes

  cout << "Converting genotypes ..." << endl;
  genotypes240file << "ID";
  genotypes300file << "ID";

  vector<bool>::const_iterator b = Study240.begin();
  for(vector<string>::const_iterator v = loci.begin(); v!= loci.end(); ++v, ++b){
    if(*b) genotypes240file << "\t" << *v;
    genotypes300file << "\t" << *v; 
  }
  genotypes240file << endl;
  genotypes300file << endl;

  rawgenofile >> i1 >> i2;
  unsigned indivindex = 0;
  while(!rawgenofile.eof()){
    //write indiv id
    ++indivindex;
    //cout << " Indiv " << indivindex << endl;
    genotypes240file << "Ind" << indivindex;
    genotypes300file << "Ind" << indivindex;

    //write genotypes line by line
    vector<bool>::const_iterator b = Study240.begin();
    vector<pair<char, char> >::iterator leg = legend.begin();
    unsigned legindex = 0;
    for(unsigned l = 0; l < NLoci; ++l){
      rawgenofile >> s;
      if(includeLocus[l]){
	//cout << "\rLocus     " << l+1 << " " << s << " " << leg->first << " " << leg->second;
	string g;
	try{
	  g = getGenotype(s, *leg);
	}
	catch(int){//found locus with mismatched alleles ie on reverse strand
	  //NB: assumes complementary bases are never paired
	  if(indivindex==1){//so that replacement is only done once
	    *leg = getComplement(*leg);//replace with reverse strand bases
	    g = getGenotype(s, *leg);//and try again
	  }
	  else {
	    cerr << "Error: anomalous genotype coding" << endl;
	    exit(1);
	  }
	}

// 	if(g[1]=='0' && s[0]!='N'){
// 	  cout << "ind " << indivindex << ", locus " << l << ", " << loci[legindex] << " " << s << " " << g << " " << legindex << " " << leg->first << " " << leg->second << endl; 
// 	  system("pause");
// 	}
	++leg;
	++b;
	++legindex;
	if(*b)genotypes240file << "\t" << g;
	genotypes300file << "\t" << g;
      }
    }
    //cout << endl;
    genotypes240file << endl;
    genotypes300file << endl;
    rawgenofile >> i1 >> i2;//skip first 2 columns
  }
  rawgenofile.close();
  genotypes240file.close();
  genotypes300file.close();
  cout << "Wrote genotypes for " << indivindex << " individuals" << endl;

  return 0;
}
