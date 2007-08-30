/**
   hapmix2impute.cpp
   Program to convert data from HAPMIXMAP format to IMPUTE format

   This program is free software distributed WITHOUT ANY WARRANTY. 
   You can redistribute it and/or modify it under the terms of the GNU General Public License, 
   version 2 or later, as published by the Free Software Foundation. 
   See the file COPYING for details.

   Copyright (c) David O'Donnell 2007
*/
#include <iostream>
#include <fstream>
#include "HapMapLegend.h"
#include "bclib/OptionReader.h"
#include "bclib/StringSplitter.h"
#include <string>

using namespace::std;
using namespace::bclib;

void PrintCopyrightInfo(){
  cout << "hapmix2impute version 1.0" << endl
       << "converts HAPMIXMAP data to IMPUTE format" << endl 
       << "This program is free software distributed WITHOUT ANY WARRANTY." 
       << endl 
       << "You can redistribute it and/or modify it under the terms of the " 
       << endl
       << "GNU General Public License, version 2 or later, as published " 
       << endl
       << "by the Free Software Foundation. See the file COPYING for details." 
       << endl << endl
       << "Copyright (c) David O'Donnell 2007"<< endl
       << endl;
}

void PrintUsage(){
  cout << "Usage (* denotes required argument):" << endl
       << "  hapmix2impute [arg=value ...]" << endl
       << "shortname longname       description" << endl
       << "-------------------------------------------------------" << endl
       << "  h       help           print this help and exit" << endl
       << "  v       verbose        print verbose output" << endl
       << "  V       version        print version number and exit" << endl
    //<< "  i       inputdir       directory with input data" << endl
       << "  o       outputdir      directory to write output data" << endl
       << "* l       legendin       input hapmap legend file" << endl
       << "  m       mapin          input map file" << endl
       << "* g       traingenotypes HAPMIXMAP genotypesfile with haploid training data" << endl
       << "* t       testgenotypes  HAPMIXMAP testgenotypesfile with diploid genotypes" << endl
       << "  H       haplotypefile  IMPUTE haplotype file to write" << endl
       << "  G       genotypefile   IMPUTE genotype file to write" << endl
       << "  L       legendout      output legend file" << endl
       << "  M       mapout         output map file" << endl
       << endl << endl;
}

int main(int argc, char** argv){

  OptionReader options;
  //some flags
  options.addFlag('h', "help");
  options.addFlag('v', "verbose");
  options.addFlag('V', "version");

  //input file names
  //string inputDir=".";
  //options.addOption('i', "inputdir", stringOption, &inputDir);
  string inLegendFilename;
  options.addOption('l', "legendin", stringOption, &inLegendFilename, true);
  string inMapFilename;
  options.addOption('m', "mapin", stringOption, &inMapFilename);
  string hapmixGenotypesFilename;
  options.addOption('g', "traingenotypes", stringOption, &hapmixGenotypesFilename, true);
  string testGenotypesFilename;
  options.addOption('t', "testgenotypes", stringOption, &testGenotypesFilename, true);

  //output file names (optional)
  string outputDir=".";
  options.addOption('o', "outputdir", stringOption, &outputDir);
  string haplotypeFilename="haplo.txt";
  options.addOption('H', "haplotypefile", stringOption, &haplotypeFilename);
  string genotypeFilename="geno.txt";
  options.addOption('G', "genotypefile", stringOption, &haplotypeFilename);
  string outMapFilename="map.txt";
  options.addOption('M', "mapout", stringOption, &outMapFilename);
  string outLegendFilename="legend.txt";
  options.addOption('L', "legendout", stringOption, &outLegendFilename);

  if(!options.ReadCommandLineArgs(argc, argv))
    exit(1);
  if(options.getFlag( "help" ) || options.getFlag( 'h' )){
    PrintUsage();
    exit(1);
  }
  if(options.getFlag("version")){
    PrintCopyrightInfo();
    exit(0);
  }

  if(!options.SetOptions() | !options.CheckRequiredOptions())
    exit(1);

  string scrap, header;

  // haplotype file
  ifstream inHaplotypes((/*inputDir + "/" +*/hapmixGenotypesFilename).c_str());
  if(!inHaplotypes.is_open()){
    cerr << "ERROR: cannot open " 
	 << /*inputDir << "/" <<*/ hapmixGenotypesFilename << endl;
    exit(1);
  }
  ofstream outHaplotypes((outputDir + "/" + haplotypeFilename).c_str());
  if(!outHaplotypes.is_open()){
    cerr << "ERROR: cannot open " 
	 << outputDir << "/" << haplotypeFilename << endl;
    exit(1);
  }

  if(options.getFlag("verbose"))
    cout << "Writing haplotype file ... ";
  vector<string> locusnames;
  inHaplotypes >> scrap;//scrap header in first col
  getline(inHaplotypes, header);
  StringSplitter::Tokenize(header, locusnames, " \t");

  //names of first and last loci. 
  //We need these later for the map and legend files
  const unsigned NumLoci = locusnames.size();
  const string firstlocus = locusnames[0];
  const string lastlocus = locusnames[NumLoci-1];
  locusnames.clear();

  //we now read the haplotypes and write them transposed
  vector<unsigned short > haploArray;//array to store haplotypes
  unsigned short g;
  unsigned NumHaplotypes = 0;
  while(!inHaplotypes.eof()){
    //scrap ID
    inHaplotypes >> scrap;
    if(inHaplotypes.eof())
      break;
    for(unsigned i = 0; i < NumLoci && !inHaplotypes.eof(); ++i){
      inHaplotypes >> g;
      haploArray.push_back(g);
    }
    ++NumHaplotypes;
  };
  inHaplotypes.close();

  //format of haplotypefile is: one line per SNP, one col per gamete
  for(unsigned j = 0; j < NumLoci; ++j){
    for(unsigned i = 0; i < NumHaplotypes; ++i)
      outHaplotypes << haploArray[i*NumLoci +j]-1 << " ";
    outHaplotypes << '\n';
  }

  outHaplotypes.close();
  haploArray.clear();
  if(options.getFlag("verbose")){
    cout << "Done. " 
	 << NumLoci << " loci, " << NumHaplotypes << " haplotypes" << endl;
  }

  //write legend file
  HapMapLegend Legend((/*inputDir + "/" + */inLegendFilename).c_str());

  ofstream outLegend((outputDir + "/" + outLegendFilename).c_str());
  if(!outLegend.is_open()){
    cerr << "ERROR: cannot open " 
	 << outputDir << "/" << outLegendFilename << endl;
    exit(1);
  }

  if(options.getFlag("verbose"))
    cout << "Writing legend file ... ";

  const unsigned firstindex = Legend.getIndex(firstlocus);
  const unsigned lastindex = Legend.getIndex(lastlocus);
  //write header
  outLegend << "rs position X0 X1\n";
  for(unsigned j = firstindex; j <= lastindex; ++j){
    Legend[j].print(outLegend);
  }

  outLegend.close();
  if(options.getFlag("verbose"))
    cout << "Done.\n" ;

  //write genotype file
  ifstream inGeno((/*inputDir + "/" + */testGenotypesFilename).c_str());
  if(!inGeno.is_open()){
    cerr << "ERROR: cannot open " 
	 << /*inputDir << "/" <<*/ testGenotypesFilename << endl;
    exit(1);
  }
  ofstream outGeno((outputDir + "/" + genotypeFilename).c_str());
  if(!outGeno.is_open()){
    cerr << "ERROR: cannot open " 
	 << outputDir << "/" << genotypeFilename << endl;
    exit(1);
  }
  if(options.getFlag("verbose"))
    cout << "Writing genotype file ... ";

  //read header
  inGeno >> scrap;
  getline(inGeno, header);
  //count number of typed loci
  StringSplitter::Tokenize(header,locusnames, " \t");
  const unsigned NumTypedLoci = locusnames.size();
 
  //read genotypes into a temp array as they must be transposed and converted on output
  vector<string> genoArray;
  unsigned NumTypedIndividuals = 0;
  while(!inGeno.eof()){
    inGeno >> scrap;//scrap ID
    if(inGeno.eof())
      break;
    ++NumTypedIndividuals;
    string temp;
    for(unsigned j = 0; j < NumTypedLoci && !inGeno.eof(); ++j){
      inGeno >> temp;
      genoArray.push_back(temp);
    }
  }
  inGeno.close();

  //define mapping from genotypes as pairs of integers to triplets
  map<string, string> genomap;
  //unquoted versions
  genomap["0,0"] = "0 0 0";
  genomap["1,1"] = "1 0 0";
  genomap["1,2"] = "0 1 0";
  genomap["2,2"] = "0 0 1";
  //quoted versions (saves having to dequote)
  genomap["\"0,0\""] = "0 0 0";
  genomap["\"1,1\""] = "1 0 0";
  genomap["\"1,2\""] = "0 1 0";
  genomap["\"2,2\""] = "0 0 1";


  //write to output file. Format is :
  // locus_ID rs# position allele0 allele1 p0 p1 p2 p0 p1 p2 ...
  for(unsigned locus = 0; locus < NumTypedLoci; ++locus){
    outGeno << "Locus" << locus+1 << " ";
    Legend[locusnames[locus]].print(outGeno);
    for(unsigned i = 0; i < NumTypedIndividuals; ++i){
      outGeno << " " << genomap[genoArray[i*NumTypedLoci + locus]];
    }
    outGeno << '\n';
  }
  outGeno.close();

  if(options.getFlag("verbose"))
    cout << "Done. " 
	 << NumTypedLoci << " typed loci, " 
	 << NumTypedIndividuals << " typed individuals\n" ;

  genomap.clear();
  genoArray.clear();
  locusnames.clear();
  Legend.clear();

  //write map file
  if(inMapFilename.size()){
    ifstream inMap((/*inputDir + "/" + */inMapFilename).c_str());
    if(!inMap.is_open()){
      cerr << "ERROR: cannot open " 
	   << /*inputDir << "/" <<*/ inMapFilename << endl;
      exit(1);
    }
    ofstream outMap((outputDir + "/" + outMapFilename).c_str());
    if(!outMap.is_open()){
      cerr << "ERROR: cannot open " 
	   << outputDir << "/" << outMapFilename << endl;
      exit(1);
    }
    if(options.getFlag("verbose"))
      cout << "Writing map file ... ";

    //read and write header  
    getline(inMap, header);
    outMap << header;

    unsigned linenumber = 0;
    string line;
    while(linenumber <= lastindex && !inMap.eof()){
      getline(inMap,line);
      if(linenumber >= firstindex)
	outMap << line << '\n'; 
      ++linenumber;
    }

    outMap.close();
    inMap.close();
    if(options.getFlag("verbose"))
      cout << "Done.\n " ;
  }

  //finished
  if(options.getFlag("verbose"))
     cout << "Finished" << endl;
  return 0;
}

