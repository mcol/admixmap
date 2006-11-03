
/*
Purpose: to set genotypes of some of the HapMap individuals Illumina/Affy loci to missing.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <getopt.h>
using namespace::std;

//small function to determine if a given string is in a list of strings
bool isListedSNP(const string HapMapSNP, const vector<string>SNPs){
  bool b = false;
  if(find(SNPs.begin(), SNPs.end(), HapMapSNP) < SNPs.end())b = true;
  return b;
}

void Tokenize(const string& str,
                      vector<string>& tokens,
                      const string& delimiters = " ")
{
  // Skip delimiters at beginning.
  string::size_type lastPos = 0;
  // Find first "non-delimiter".
  string::size_type pos     = 0;
  string::size_type paren = 0;
  string::size_type closeparen = 0;
  do    {
    // Skip delimiters.  Note the "not_of"
    lastPos = str.find_first_not_of(delimiters, pos);
    
    // Find next "non-delimiter"
    pos = str.find_first_of(delimiters, lastPos);
    
    paren = str.find_first_of("\"", lastPos);
    if(string::npos != paren ){
      closeparen = str.find_first_of("\"", paren+1);
      if(string::npos == closeparen) throw string("missing closing quotes");
      if(paren < pos && pos < closeparen){//skip delimiters within quotes
	//++lastPos;
	pos =closeparen+1;
      }
    }
    if(string::npos != pos || string::npos != lastPos)
      // Found a token, add it to the vector.
      tokens.push_back(str.substr(lastPos, pos - lastPos));
  }
  
  while (string::npos != pos || string::npos != lastPos);
}

void PrintHelpText(){
    cout << "makegenotypes: sets genotypes of some individuals in a HapMap genotypes file" 
	 << " to missing at loci not included in a given list\n";
    cout << "Usage: makegenotypes -g<genotypesfile> -l<SNPlist> [-n<#individuals>] [-o<outputfile>]\n";
}

int main(int argc, char** argv){
    if(argc < 3)PrintHelpText();
    ifstream SNPFile;//("Ill300Chr22rsnumbers.txt");
    ifstream ingenotypes;//("Eur/chr22data/genotypes.txt");
    ofstream outgenofile;//("Eur/chr22data/Illumina300_genotypes.txt");
    //ofstream outgenofile("chr22data/Affy_genotypes.txt");
    string s;
    int NUMCASES=10;

    int ich;
    while((ich = getopt(argc, argv, "hn:o:g:l:")) != EOF){
	switch(ich){
	    case 'h':{
		PrintHelpText();
		break;
	    }
	    case 'n':{
		NUMCASES = atoi(optarg);
		break;
	    }
	    case 'o':{//outputfile
		outgenofile.open(optarg);
		cout << "Writing output to " << optarg << endl;
		break;
	    }
	    case 'g':{//input genotypes file
		ingenotypes.open(optarg);
		break;
	    }
	    case 'l':{//list of SNPs
		SNPFile.open(optarg);
		break;
	    }
	    default:{
		cout << "Invalid args\n";
		break;
	    }
	}
    }
    if(!ingenotypes.is_open()){
	cout << "Error: must specify HapMap genotypesfile\n";
	exit(1);
    }
    if(!outgenofile.is_open()){
	outgenofile.open("genotypes.txt");
	cout << "Writing output to genotypes.txt" << endl;
    }
    if(!SNPFile.is_open()){
	cout << "Error: must specify a list of SNPs\n";
	exit(1);
    }


    if(optind < argc){
	cout << "non-option ARGV elements: ";
	while(optind < argc)cout << argv[optind++] << "\t";
	cout << endl;
	exit(1);
    }

  // 1: read supplied rs numbers
  string SNP;
  vector<string> SNPList;
  SNPFile >> SNP;

  while(!SNPFile.eof()){
    SNPList.push_back(SNP);
    SNPFile >> SNP;
  }

  //2: read header of HapMap genotypes file and store as vector 
  ingenotypes >> s;//discard label in first column
  outgenofile << s << "\t"; 

  vector<string> HapMapSNPs;
  getline(ingenotypes, s);//read in rest of header
  outgenofile << s << endl;//write to file
  Tokenize(s, HapMapSNPs, "\t");//break into seperate strings

  vector<bool> ListedSNPs;
  for(vector<string>::const_iterator i = HapMapSNPs.begin(); i != HapMapSNPs.end(); ++i){
      bool b = isListedSNP(*i, SNPList);
      ListedSNPs.push_back(b);
      //if(b)cout << " --> " << HapMapSNPs[i] << endl; 
  }

  const size_t NumSNPs = HapMapSNPs.size();
  cout << NumSNPs << " HapMap SNPs, " << SNPList.size() << " listed SNPs" << endl;
  HapMapSNPs.clear();
  SNPList.clear();

  //system("PAUSE");
//determine number of individuals (lines in file)
  ofstream::pos_type line1 = ingenotypes.tellg();
  getline(ingenotypes, s);
  int NUMHAPMAPINDS = 0;
  while(!ingenotypes.eof()){
      ++NUMHAPMAPINDS;
      getline(ingenotypes, s);
  }
  cout << NUMHAPMAPINDS << " individuals in genotypes file. Modifying genotypes of last " << NUMCASES << endl;
  if(NUMHAPMAPINDS <= NUMCASES) {
      cout << "Error: NUMCASES > #individuals in genotyes file\n";
      exit(1);
  }
  ingenotypes.clear();
  ingenotypes.seekg(line1);

  //3: read all but last 10 individuals and write straight to file
  for(unsigned line  = 0; line < NUMHAPMAPINDS-NUMCASES; ++line){
      //cout << "Individual " << line+1 << endl;
    getline(ingenotypes, s);
    outgenofile  << s << endl;
  }

  //4: read remaining individuals, replacing unlisted SNPS by missing genotypes
  const string missing = "\"0,0\"";

  for(unsigned line = 0; line < NUMCASES; ++line){
      //cout << "Individual " << NUMHAPMAPINDS-NUMCASES+line+1 << endl;
    ingenotypes >> s; //read individual ID
    outgenofile << s << "\t";//write to new file

    for(size_t  i = 0; i < NumSNPs; ++i){
      ingenotypes >> s; //read in genotype
      if(ListedSNPs[i]){//if SNP in list
	  outgenofile << s;//output genotype to file
      }
      else{//output missing genotype
	  outgenofile << missing;
      }
      outgenofile << "\t";
    }
    outgenofile << endl;
  }

  cout << "output complete\n\n";
  return 0;
}
