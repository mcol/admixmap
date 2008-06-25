/**
 * prog to convert hapmap data to ADMIXMAP format
 * supply chromosome number as arg. 0 (default) means all.
 * also requires infofile and data directory (Eur, Afr or Asian).
 * Optionally supply genotypes file and locus file names
 * (default to genotypes.txt) loci.txt
 *
 * NB: input data file must be named "chr#.txt",
 * where # is a number from 1 to 22, and all be located
 * in the data directory
 *
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 *
 * Copyright (c) David O'Donnell 2006
 *
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <getopt.h>
#include <cstdlib>
#include <algorithm>

unsigned NUMIND = 90;
unsigned NUMVALIDINDS =0;
unsigned long NUMLOCI = 0;
unsigned long NUMBADLOCI = 0;
unsigned long NUMBADLOCIONCHR = 0;
unsigned NUMALLELES = 2;
bool be_quiet = false;

using namespace::std;
string getGenotype(string g, string a){
  //converts a base pair string g, eg "CT", to an admixmap format genotype eg "1,2"
  //a is the reference string of alleles eg "C/T"
  stringstream s;
  s << "\"" << (int)(g[0]==a[0])+2*(int)(g[0]==a[2]) << "," << (int)(g[1]==a[0])+2*(int)(g[1]==a[2]) << "\"";
  return s.str();
}

void PrintHelpText(){
    cout << "----------------------------------------------------------------------------\n" 
	 << "This program converts HapMap data to ADMIXMAP format"
	 << "Usage: convertdata -c<chrm> -i<infofile> [-p<pop>] [-g<genotypesfile>]"
	 << " [-l<locusfile>] [-n<numloci>" << endl
	 << "where pop is the directory containing data files ('Eur', 'Afr' or 'Asian')," 
	 << " defaults to current directory; " << endl
	 << "chrm is a chromosome number. -c0 converts all chromosomes." << endl
	 << "numloci is the number of loci per chromosome to use, default is all\n"
	 << "genotypesfile and locusfile default to 'genotypes.txt' and 'loci.txt'" << endl
	 << "Note that HapMap input files must be named 'chrN.txt', where N is an integer from 1 to 22\n\n";
}

int main(int argc, char **argv){
    if(argc<3){
	PrintHelpText();
	exit(0);
    }
    int ich;
    unsigned CHRNUM = 0;//chromosome number
    string popprefix = ".";
    ifstream infofile;
    ifstream genotypesin;
    ofstream genotypesfile;
    ofstream locusfile;
    unsigned long userloci = 1000000000;//some number > number of HapMap loci

    while((ich=getopt(argc, argv, "hi:c:g:l:p:n:q"))!=EOF){
	switch(ich){
	case 'h' :{
	    PrintHelpText();
	}
	    case 'c':{
		CHRNUM = atoi(optarg);
		break;
	    }
	    case 'p':{
		popprefix = optarg;
		break;
	    }
	    case 'i':{
		infofile.open(optarg);
		if(!infofile.is_open()){
		    cout << "Error: cannot open infofile\n";
		    exit(1);
		}
		break;
	    }
	    case 'g':{
		genotypesfile.open(optarg);
		if(!genotypesfile.is_open()){
		    cout << "Error: cannot open genotypesfile\n";
		    exit(1);
		}
		cout << "Writing genotypes to " << optarg << endl;
		break;
	    }
	    case 'l':{
		locusfile.open(optarg);
		if(!locusfile.is_open()){
		    cout << "Error: cannot open locusfile\n";
		    exit(1);
		}
		cout << "Writing locusfile to " << optarg << endl;
		break;
	    }
	    case 'n':{//number of loci
		const unsigned temp = atoi(optarg);
		if(temp>0)userloci = temp;
		break;
	    }
	    case 'q':{//quiet mode
		be_quiet = true;
		break;
	    }
	    default:{
		cout << "Invalid args specified: " << ich << endl;
		exit(1);
	    }

	}
    }
    if(!infofile.is_open()){
	cout << "Error: must specify infofile" << endl;
	exit(1);
    }
    if(!genotypesfile.is_open())genotypesfile.open("genotypes.txt");
    if(!locusfile.is_open())locusfile.open("loci.txt");
    
    
//     CHRNUM = atoi(argv[1]);
//     ifstream infofile(argv[2]);//HapMap file with info on individuals' sex and validity (only unrelated inds are valid)
//   ifstream genotypesin;
//   ofstream genotypesfile;
//   if(argc>3)genotypesfile.open(argv[3], ios::out);
//   else genotypesfile.open("genotypes.txt", ios::out);
//   ofstream locusfile;
//   if(argc>4)locusfile.open(argv[4]);
//   else locusfile.open("loci.txt");

  ofstream badlocusfile("excludedSNPs.txt");
  //ofstream outcomefile("outcome_halfmissing.txt");
  string alleles;
  string genotype;
  string obs;//observed genotype
  string SNPID;
  vector<bool> BadLoci;
  double position = 0.0, prev = 0.0;
  string scrap;
  unsigned int locus = 0;
  unsigned int indiv = 0;
  // unsigned chromosome = 0;
  unsigned lastchr = 22;
  if(CHRNUM==0)CHRNUM = 1;
  else lastchr = CHRNUM;

// ************************** PHASE 1:  read info file for indiv id, sex, pedigree *********************
/*
  Columns of the info file are:
  1: family PED ID (ignored) (omitted in unrelateds)
  2: individual PED ID (ignored)
  3: father PED ID (0 if founder)
  4: mother PED ID (0 if founder)
  5: sex (1=male, 2=female)
  6: individual LSID (ignored)
  7: sample LSID (assigned by DCC based on Coriell catalog no.; used to assign individual ID in genotypes file)
  
  We want unrelated individuals i.e. with 0 in cols 3 and 4
*/
  if(!infofile.is_open()){
      cout << "Unable to open infofile\n";
      exit(1);
  }

  vector<string> INDIVID;
  //vector<int> indisvalid;
  //vector<int> sex;//1 = male, 2 = female
  map<string, pair<int, int> > IndivInfo; 

  NUMIND = 0;
  bool related;
  if(popprefix == "Eur" || popprefix =="Afr")related = true;
  else if(popprefix == "Asian")related = false;//Asians are all unrelated
  else{
//determine if individuals in file are related by counting tab characters in first row
      getline(infofile, scrap);
      string::size_type pos = 0;
      unsigned ncols = 0;
      while(pos != string::npos){
	  pos = scrap.find_first_not_of("\t ", pos);
	  if(pos != string::npos){
	      ++ncols;
	      pos = scrap.find_first_of("\t ", pos);//reposition at start of next whitespace
	  }
      }
      if(ncols==6)related = false;
      else if(ncols==7)related = true;
      else {
	  cout << "wrong number of cols in infofile "<< endl;
	  exit(1);
      }
      //cout << ncols << "individuals are " << (related? "related" : "unrelated")<< endl;
      infofile.seekg(0);//reset to start of file
  }

  while(!infofile.eof()){
    int _sex, valid;
    string _id;
    int father_ped, mother_ped;
    if(related)infofile >> scrap;//skip family info ID
    infofile >> scrap >> father_ped >> mother_ped >> _sex >> scrap >> _id;

    if(_id.find_first_not_of(" \t\n\r")!=string::npos){//check for empty lines
	if(father_ped ==0 && mother_ped ==0)valid = 1;else valid = 0;
	string::size_type pos0 = _id.find("N");//ID begins with NA and ends with :1
	string::size_type pos1 = _id.find_last_of(":");

	INDIVID.push_back(_id.substr(pos0, pos1-pos0));
	//sex.push_back(_sex);
	//indisvalid.push_back(valid);
	IndivInfo[_id.substr(pos0, pos1-pos0)] = pair<int, int>(_sex, valid);
	//cout << _id.substr(pos0, pos1-pos0) << " " << valid << endl;
	++indiv;
	if(valid)++NUMVALIDINDS;
    }
  }
  sort(INDIVID.begin(), INDIVID.end());
//NB: assuming individual IDs in genotypes files are are sorted numerically
  NUMIND = indiv;
  infofile.close();
  cout << NUMIND << " individuals" <<endl;
  cout << NUMVALIDINDS << " unrelated individuals" << endl;

  //exit(0);

// ************************** PHASE 2: write locus file and header of genotypesfile, count number of loci  *********************
  genotypesfile << "\"Individ\"\t";
  locusfile << "\"SNPid\"\t\"NumAlleles\"\t\"DistanceinMb\"\n";//locusfile header
  locusfile << setiosflags(ios::fixed) << setprecision(8);

  long unsigned abslocus = 0;//absolute locus number
  for(unsigned chr = CHRNUM; chr <= lastchr; ++chr){
      cout << "\nChromosome " << chr << "  " <<endl;
      NUMBADLOCIONCHR = 0;
      locus = 0;
      stringstream ss;

/// unzip data
    //ss << "gzip -d chr" << chr << ".gz";
    //cout << "\nUnzipping..."<<flush;
    //system(ss.str().c_str());
    //ss.clear();
///open HapMap genotypes file
/// NB: files should be named chr1, chr2 etc
    ss << popprefix<< "/chr" << chr << ".txt";//"_sorted";
    genotypesin.clear();//to clear fail status at eof
    genotypesin.open(ss.str().c_str());
    if(!genotypesin.is_open()){
	cout << "Could not open genotypes file " << ss.str() <<endl;
	exit(1); 
    }

    //cout << "\nReading raw genotypes...\n"<<flush;
    getline(genotypesin, scrap);//skip header
    while ( (!genotypesin.eof())  && (locus < userloci) ){
	////use this to split up large locusfiles
	//if(!(locus%10000)){
      //locusfile.close();
      //stringstream s; s << "SNP" << locus +10000  << ".txt";
      //locusfile.open(s.str().c_str());
      //}
      cout << "\rLocus    " << locus+1 << flush;
      //we want cols: 0(snpid), 1(alleles), 3(position in basepairs), 11 onwards(indiv genotypes)
      genotypesin >> SNPID >> alleles
		  >> scrap >> position;
      unsigned indiv_index = 0;
      if(SNPID.find_first_not_of(" \t\n\r")!=string::npos){//check for empty lines
//check for and skip monomorphic loci (where all genotypes are the same and homozygous)
	  for(unsigned ii = 0; ii < 7; ++ii)//skip 7 columns
	      genotypesin >> scrap;
	  bool ismono = true;
	  string firstgeno;
	  do{//find first non-missing genotype
	      genotypesin >> firstgeno;
	      ++indiv_index;
	  }
//keep going until we find a founder with nonmissing genotype or run out of individuals
//(NN denotes missing genotype)
	  while( (!IndivInfo[INDIVID[indiv_index]].second || firstgeno[0]=='N') && (indiv_index < NUMIND));
	  if(indiv_index< NUMIND && (firstgeno[0] == firstgeno[1])){
	      string nextgeno;
	      do{
		  genotypesin >> nextgeno;
		  if(IndivInfo[INDIVID[indiv_index]].second)//only consider founders
		      ismono = ismono & ((nextgeno[0]=='N') || (nextgeno == firstgeno));
		  ++indiv_index;
	      }
//keep going until we find an individual with nonmissing genotype different from the first or run out of individuals
	      while((indiv_index < NUMIND) && ismono);
	  }

	//write locusfile
	  if(!ismono && (position-prev)>0){//strictly greater than to avoid having comp loci
	      locusfile << SNPID << "\t" <<  2 << "\t";
	      if(locus==0) locusfile << "#" ;//missing value for first distance on chr
	      else 
//converts position in basepairs to distance in centiMorgans 
		  //locusfile << 0.0000013 * (position-prev);
//convert position to distance in megabases
		  locusfile << 0.000001 * (position-prev);
	      locusfile << /*chrnumber <<*/ endl;
	      //write header of genotypesfile
	      genotypesfile << SNPID << "\t";
	      ++locus;
	      BadLoci.push_back(false);
	      prev = position;
	  }
	  else{//locus out of sequence
	    ++NUMBADLOCIONCHR;
	    ++NUMBADLOCI;
	    BadLoci.push_back(true);//excludes loci out of sequence
	    badlocusfile << SNPID << "\t" << chr << endl;
	}
	  ++abslocus;
      }
      SNPID.clear();       
      //skip rest of line
      if(indiv_index< NUMIND)getline(genotypesin, scrap);
    }
    genotypesin.close();
    cout << endl << locus << " Loci, " << NUMBADLOCIONCHR  << " loci excluded";
    NUMLOCI += locus;//count number of loci
  }
   locusfile.close();
   cout << "\nFinished writing locusfile. " << NUMLOCI <<  " loci" << endl;
   //if(BadLoci.size()){
       cout << NUMBADLOCI/*BadLoci.size()*/ << " loci excluded: ";
     //for(vector<string>::const_iterator i = BadLoci.begin(); i < BadLoci.end(); ++i)cout << *i << " "; 
     cout << endl;
     //}
     badlocusfile.close();

   genotypesfile << endl;

// ************************** PHASE 3: Write genotypesfile  *********************
   
   for(indiv = 0; indiv < NUMIND; ++indiv)
     //if individual is a founder
     if(IndivInfo[INDIVID[indiv]].second) {
       if(!be_quiet)
	   cout << "\n" << indiv+1 << " " << INDIVID[indiv] << " "  << flush;
      //write indiv id and sex
      genotypesfile << INDIVID[indiv] << "\t";// << sex[indiv] <<"\t";
      position = 0.0;     
//       if(indiv%2) outcomefile << "#\n";
//       else outcomefile << 9 << endl;
      abslocus = 0;
      for(unsigned chr = CHRNUM; chr <= lastchr; ++chr){
	  if(!be_quiet)
	      cout << "\nChromosome " << chr << "  "<< flush ;
	  stringstream ss;
	  ss << popprefix << "/chr" << chr << ".txt";//"_sorted";
	  genotypesin.open(ss.str().c_str());
 
	  //open genotypes input file
	  genotypesin.clear(); // to clear fail status caused by eof
	  //genotypesin.seekg(0);
	  getline(genotypesin, scrap);//skip header
      
	  //read genotypes, one locus (row) at a time
	  //for(locus = 0; locus < NUMLOCI+ NUMBADLOCI/*BadLoci.size()*/; ++locus){
	  prev = position;
	  genotypesin >> SNPID >> alleles >> scrap >> position;
	  locus = 0;
	  while ( (!genotypesin.eof())  && locus < userloci){
//	      for(unsigned _locus = 0; _locus < NUMLOCI; ++_locus){
//	      if(position-prev >=0.0){//skip loci out of sequence
	      if(!BadLoci[abslocus]){
		  //cout << "  locus" << locus+1 << " " << SNPID << " " <<position << " " << prev << endl;
		  for(unsigned int col = 0; col < 7+indiv; ++col) {
                    genotypesin >> scrap;//skip to col for this individual, genotypes start at col 11
                  }
		  genotypesin >> obs;
		  
		  //write to genotypesfile in admixmap format
		  //set 50% to missing at half the loci
// 		  if(indiv%2 && locus%2){
// 		      genotypesfile << "\"0,0\"" << "\t";
// 		  }
// 		  else{
 		      genotypesfile << getGenotype(obs, alleles) << "\t";
// 		  }
		      ++locus;
	      }
	      getline(genotypesin, scrap);//skip remaining individuals in line
	      prev = position;//start reading next line
	      genotypesin >> SNPID >> alleles >> scrap >> position;
	      ++abslocus;
	  }
	  genotypesin.close();
      }
      
      genotypesfile << endl;
    }
  genotypesfile.close();
//  outcomefile.close();
  cout << endl << "Finished writing genotypesfile" << endl;
}
