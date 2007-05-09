/*
  FPHDOptions.cpp
  Options for FPHD
  This file is part of FPHD

  This program is free software distributed WITHOUT ANY WARRANTY. 
  You can redistribute it and/or modify it under the terms of the GNU General Public License, 
  version 2 or later, as published by the Free Software Foundation. 
  See the file COPYING for details.

  Copyright (c) David O'Donnell 2007
*/
#include "FPHDOptions.h"
#include <iostream>
#include <getopt.h>


using std::endl;
using std::cout;
#define MAXCHROMOSOMES 22//maximum number of chromosomes
#define PROGNAME "FPHD"

FPHDOptions::FPHDOptions(int argc, char*const* argv, std::ofstream& genotypesfile, std::ofstream& locusfile){
  //set defaults
  beVerbose = false;
  prefix = ".";
  incasecontrolfilename = 0;
  outcasecontrolfilename = 0;
  userloci = 1000000000;  //some number > number of HapMap loci
  LimitLoci = false;
  FirstChr = 0;
  flankLength = 1e5;//10Kb

  ParseOptions(argc, argv, genotypesfile, locusfile);
  if (FirstChr == 0){//all chromosomes (-c0 option)
    FirstChr = 1;
    lastchr = MAXCHROMOSOMES;
  }
  else//only one chromosome
    lastchr = FirstChr;

}

void FPHDOptions::ParseOptions(int argc, char*const* argv, std::ofstream& genotypesfile, std::ofstream& locusfile){
  int ich;
  while ((ich = getopt(argc, argv, "hc:g:l:p:i:o:n:f:qv")) != EOF) {
    switch (ich) {
    case 'h':{
      PrintHelpText();
    }
    case 'c':{
      FirstChr = atoi(optarg);
      break;
    }
    case 'p':{
      prefix = optarg;
      break;
    }
    case 'g':{
      genotypesfile.open(optarg);
      if (!genotypesfile.is_open()) {
        cout << "Error: cannot open genotypesfile\n";
        exit(1);
      }
      cout << "Writing genotypes to " << optarg << endl;
      break;
    }
    case 'l':{
      locusfile.open(optarg);
      if (!locusfile.is_open()) {
        cout << "Error: cannot open locusfile\n";
        exit(1);
      }
      cout << "Writing locusfile to " << optarg << endl;
      break;
    }
    case 'i':{
      incasecontrolfilename = optarg;
      break;
    }
    case 'o':{
      outcasecontrolfilename = optarg;
      break;
    }

    case 'n':{                 //number of loci
      const unsigned temp = atoi(optarg);
      if (temp > 0){
        userloci = temp - 1;
        LimitLoci = true;
      }
      break;
    }
    case 'f':{//flanking region
      //NB distance is given in Kb
      flankLength = 1000 * atof(optarg);
      break;
    }
    case 'v':{                 //verbose mode
      beVerbose = true;
      break;
    }
    default:{
      cout << "Invalid args specified: " << ich << endl;
      exit(1);
    }

    }
  }
  if (!genotypesfile.is_open())
    genotypesfile.open("genotypes.txt");
  if (!locusfile.is_open())
    locusfile.open("loci.txt");
  //set default ouput ccgenotypesfilename if none specified
  if(incasecontrolfilename && !outcasecontrolfilename)
    outcasecontrolfilename = "CaseControlGenotypes.txt";
  //ofstream outcomefile("outcome_halfmissing.txt");
}

void FPHDOptions::PrintHelpText(){
  cout << "----------------------------------------------------------------------------"
       << endl
       << "This program converts phased HapMap data to HAPMIXMAP format" << endl
       << "Copyright (c) David O'Donnell 2007" << endl
       << "All parts of this program are freely distributable" << endl << endl
       << "Usage: " << PROGNAME << " -c<chr> [options]" << endl << endl
       << "Options: " << endl
       << "-c<chr>            - chromosome number. -c0 converts all chromosomes." << endl
       << "-v                 - be verbose" << endl
       << "-p <prefix>        - prefix where HapMap files are located, defaults to '.'" << endl 
       << "-g <genotypesfile> - output genotypes file, defaults to 'genotypes.txt'" << endl
       << "-l <locusfile>     - output locus file, defaults to loci.txt" << endl
       << "-n <numloci>       - maximum number of loci per chromosome. Ignored with -i" << endl
       << "-i <input case-control file> "<< endl
       << "                   - raw case-control genotypes file. This file will be" << endl
       << "                     formatted and used to restrict the range of HapMap loci" << endl 
       << "                     output." << endl
       << "-o <output case-control-file>  " << endl
       << "                   - valid only with -i, defaults to 'CaseControlGenotypes.txt'." << endl
       << "-f<length-in-Kb>   - length in Kb of flanking region outside of typed region." << endl
       << "                     Valid only with -i. Defaults to 10." << endl << endl
       << "Note that HapMap input files must be named 'chrN_phased.txt', 'chrN_sample.txt'" << endl
       << "and 'chrN_legend.txt', where N is an integer from 1 to 22"
       << endl << endl;
}

bool FPHDOptions::Verbose()const{
  return beVerbose;
}
bool FPHDOptions::LimitedLoci()const{
  return LimitLoci;
}

const std::string& FPHDOptions::getPrefix()const{
  return prefix;
}
unsigned FPHDOptions::getFirstChr()const{
  return FirstChr;
}
unsigned FPHDOptions::getLastChr()const{
  return lastchr;
}
unsigned long FPHDOptions::getNumUserLoci()const{
  return userloci;
}

bool FPHDOptions::WriteCCFile()const{
  return (incasecontrolfilename != 0);
}
const char* FPHDOptions::getInCCFilename()const{
  return incasecontrolfilename;
}
const char* FPHDOptions::getOutCCFilename()const{
  return outcasecontrolfilename;
}
float FPHDOptions::getFlankLength()const{
  return flankLength;
}
