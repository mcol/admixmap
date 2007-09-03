/**
   \file InitialValues.cpp
   functions to split HAPMIXMAP initial value files
   This file is part of FPHD

   This program is free software distributed WITHOUT ANY WARRANTY. 
   You can redistribute it and/or modify it under the terms of the GNU General Public License, 
   version 2 or later, as published by the Free Software Foundation. 
   See the file COPYING for details.

   Copyright (c) David O'Donnell 2007
*/

#include "InitialValues.h"
#include "HapMapLegend.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iterator>

using namespace std;

//assuming K = 8
#define NUMBLOCKSTATES 8

void OpenOutputFiles(vector<ofstream*>& outfiles, const char* infilename, unsigned NumSubChrs){
  //strip extension from filename
  string infilestring = infilename;
  infilestring.erase(infilestring.find_first_of("."));

  for(unsigned j = 0; j < NumSubChrs; ++j){

    stringstream outfilename;
    //add extension
    if(NumSubChrs==1)
      outfilename << infilestring << ".txt";
    else
      outfilename << infilestring << "_" << j+1 << ".txt";
    
    outfiles.push_back(new ofstream(outfilename.str().c_str()));
    if (!outfiles[j]->is_open()) {
      stringstream err;
      err << "Error: cannot open " << outfilename << "\n";
      throw(err.str());
    }
  }
}

void WriteFiles(vector<ofstream*>& outfiles, ifstream& infile, HapMapLegend& Legend, unsigned stride)throw (string){

  double temp[stride];
  for(unsigned locus = 0; locus < Legend.size(); ++locus){
    const vector<unsigned>& subchrs = Legend[locus].subchr;
    for(unsigned s = 0; s < stride; ++s){
      if(!(infile >> temp[s])){
	throw string("too few entries in file");
      }
    }//end stride loop
    //copy read values to file
    for(vector<unsigned>::const_iterator sub = subchrs.begin(); sub != subchrs.end(); ++sub)
      copy(temp, temp + stride, ostream_iterator<double>(*(outfiles[*sub]), " "));

  }//end locus loop

  //see if anything more than whitespace left in file
  string test;
  infile >> test;
  if(test.find_first_not_of(" \t\n\r") != string::npos){
    throw string("ERROR: Too many entries in file\n");
    exit(1);
  }
}
void SplitInitialArrivalRateFile(const char* filename, HapMapLegend& Legend){
  if(filename == 0 || strlen(filename) == 0)
    return;
  ifstream infile(filename);
  if(!infile.is_open()){
    cerr << "ERROR: cannot open " << filename << endl;
    return;
  }
  try{
    vector<ofstream*> outfiles;
    OpenOutputFiles(outfiles, filename, Legend.getNumSubChromosomes());
    
    double h, beta;
    if(!(infile >> h >> beta)){
      stringstream err;
      err << "ERROR: too few values in " << filename << endl;
      throw err.str();
    }
    
    for(vector<ofstream*>::iterator i = outfiles.begin(); i != outfiles.end(); ++i)
      *(*i) << h << " " << beta << " ";
    try{
      WriteFiles(outfiles, infile, Legend, 1);
    }catch(string s){
      throw ("Error occurred while writing initial arrival rates: " + s); 
    }
    
    //clean up
    infile.close();
    for(vector<ofstream*>::iterator i = outfiles.begin(); i != outfiles.end(); ++i){
      (*i)->close();
      delete *i;
    }
    outfiles.clear();
  }
  catch(string s){
    cout << s << endl;
    return;
  }
}

void SplitInitialMixturePropsFile(const char* filename, HapMapLegend& Legend){
  if(filename ==0 || strlen(filename) == 0)
    return;

  ifstream infile(filename);
  if(!infile.is_open()){
    cerr << "ERROR: cannot open " << filename << endl;
    exit(1);
  }
  try{
    vector<ofstream*> outfiles;
    OpenOutputFiles(outfiles, filename, Legend.getNumSubChromosomes());
    
    try{
      WriteFiles(outfiles, infile, Legend, 1);
    }catch(string s){
      throw ("Error occurred while writing initial mixture props: " + s); 
    }
    
    //clean up
    infile.close();
    for(vector<ofstream*>::iterator i = outfiles.begin(); i != outfiles.end(); ++i){
      (*i)->close();
      delete *i;
    }
    outfiles.clear();
  }
  catch(string s){
    cout << s << endl;
    return;
  }
}

void SplitInitialAlleleFreqFile(const char* filename, HapMapLegend& Legend){
  if(filename ==0 || strlen(filename) == 0)
    return;

  ifstream infile(filename);
  if(!infile.is_open()){
    cerr << "ERROR: cannot open " << filename << endl;
    exit(1);
  }
  try{
    vector<ofstream*> outfiles;
    OpenOutputFiles(outfiles, filename, Legend.getNumSubChromosomes());
    try{
      WriteFiles(outfiles, infile, Legend, NUMBLOCKSTATES);
    }catch(string s){
      throw ("Error occurred while writing initial allele freqs: " + s); 
    }
    
    //clean up
    infile.close();
    for(vector<ofstream*>::iterator i = outfiles.begin(); i != outfiles.end(); ++i){
      (*i)->close();
      delete *i;
    }
    outfiles.clear();
  } catch(string s){
    cout << s << endl;
    return;
  }
}

void SplitInitialFreqPriorFile(const char* filename, HapMapLegend& Legend){
  if(filename ==0 || strlen(filename) == 0)
    return;

  ifstream infile(filename);
  if(!infile.is_open()){
    cerr << "ERROR: cannot open " << filename << endl;
    exit(1);
  }
  try{
    vector<ofstream*> outfiles;
    OpenOutputFiles(outfiles, filename, Legend.getNumSubChromosomes());
    try{
      WriteFiles(outfiles, infile, Legend, 2);
    }catch(string s){
      throw ("Error occurred while writing initial allele freq priors: " + s); 
    }
    
    //clean up
    infile.close();
    for(vector<ofstream*>::iterator i = outfiles.begin(); i != outfiles.end(); ++i){
      (*i)->close();
      delete *i;
    }
    outfiles.clear();
  } catch(string s){
    cout << s << endl;
    return;
  }
}
