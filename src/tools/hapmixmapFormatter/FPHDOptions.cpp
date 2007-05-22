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
//#include "bcppcl/OptionReader.h"

#define MAXCHROMOSOMES 22//maximum number of chromosomes
#define PROGNAME "FPHD"

using namespace::std;

FPHDOptions::FPHDOptions(int argc, char** argv, std::ofstream& genotypesfile, std::ofstream& locusfile){
  //set defaults
  beVerbose = false;
  prefix = ".";
  incasecontrolfilename;
  outcasecontrolfilename;
  userloci = 1000000000;  //some number > number of HapMap loci
  LimitLoci = false;
  Chr = 0;
  flankLength = 10;//10Kb
  MinOverlap = 5000.0;
  MissingChar = 'N';

  if(argc == 1){//no args specified - print usage and exit
    PrintHelpText();
    exit(1);
  }
  ParseOptions(argc, argv, genotypesfile, locusfile);
  if (Chr <= 0 || Chr > MAXCHROMOSOMES){
    cerr << "ERROR: Invalid chromosome number: " << Chr << endl;
    exit(1);
  }

}

void FPHDOptions::ParseOptions(int argc, char** argv, std::ofstream& genotypesfile, std::ofstream& locusfile){
  bcppcl::OptionReader opt;
  opt.setVerbose(true); /* print warnings about unknown options */

  string genotypesfilename, locusfilename;
  unsigned long locuslimit = 0;
  opt.setUserOption("genotypesfile", "genotypes.txt");
  opt.setUserOption("locusfile", "loci.txt");
  DefineOptions(opt, &genotypesfilename, &locusfilename, &locuslimit);
  opt.ReadCommandLineArgs(argc, argv);
  if(opt.getFlag( "help" ) || opt.getFlag( 'h' )){
    PrintHelpText();
    exit(1);
  }
  if(!opt.SetOptions() | !opt.CheckRequiredOptions())
    exit(1);

  genotypesfile.open(genotypesfilename.c_str());
  if(!genotypesfile.is_open()){
    cout << "Error: cannot open genotypesfile\n";
    exit(1);
  }
  locusfile.open(locusfilename.c_str());
  if (!locusfile.is_open()) {
    cout << "Error: cannot open locusfile\n";
    exit(1);
  }
  //check for upper limit on loci
  if (locuslimit > 0){
    userloci = locuslimit - 1;
    LimitLoci = true;
  }
  //convert flankLength in kb to bp
  flankLength *= 1000;

  //set verbose flag
  beVerbose = opt.getFlag("verbose");

  //set default ouput ccgenotypesfilename if none specified
  if(incasecontrolfilename.size() && !outcasecontrolfilename.size())
    outcasecontrolfilename = "CaseControlGenotypes.txt";

  if(beVerbose){
    cout << "Writing genotypes to " << genotypesfilename << endl;
    cout << "Writing locusfile to " << locusfilename << endl;
    if(incasecontrolfilename.size())
      cout << "Writing case-control genotypes to " << outcasecontrolfilename << endl; 
  }
}

void FPHDOptions::DefineOptions(bcppcl::OptionReader& opt, string* genotypesfilename, string* locusfilename, unsigned long* locuslimit){
  opt.addFlag('h', "help");
  opt.addFlag('v', "verbose");
  opt.addOption('c', "chromosome", bcppcl::intOption, &Chr, true);
  opt.addOption('p', "prefix", bcppcl::stringOption, &prefix);
  opt.addOption('g', "genotypesfile", bcppcl::stringOption, genotypesfilename);
  opt.addOption('l', "locusfile", bcppcl::stringOption, locusfilename);
  opt.addOption('n', "numloci", bcppcl::longOption, &locuslimit);
  opt.addOption("minoverlap", bcppcl::floatOption, &MinOverlap);
  opt.addOption('i', "inputfile", bcppcl::stringOption, &incasecontrolfilename);
  opt.addOption('o', "outputfile", bcppcl::stringOption, &outcasecontrolfilename);
  opt.addOption('f', "flank", bcppcl::floatOption, &flankLength);
  opt.addOption('m', "missing", bcppcl::charOption, &MissingChar);
}


void FPHDOptions::PrintHelpText(){
  cout << "----------------------------------------------------------------------------"
       << endl
       << "This program converts phased HapMap data to HAPMIXMAP format" << endl
       << "Copyright (c) David O'Donnell 2007" << endl
       << "All parts of this program are freely distributable" << endl << endl
       << "Usage: " << PROGNAME << " -c<chr> [options]" << endl << endl
       << "Options ('-' signs are optional): " << endl
       << "-h   -help            print this help message and exit" << endl
       << "-c   -chromosome      chromosome number" << endl
       << "-v   -verbose         be verbose" << endl
       << "-p   -prefix          prefix where HapMap files are located, defaults to '.'" << endl 
       << "-g   -genotypesfile   output genotypes file, defaults to 'genotypes.txt'" << endl
       << "-l   -locusfile       output locus file, defaults to loci.txt" << endl
       << "-n   -numloci         maximum number of loci per (sub)chromosome. " << endl
       << "     -minoverlap      minimum overlap between sub-chromosomes in kb." << endl
       << "-i   -inputfile       raw case-control genotypes file. This file will be" << endl
       << "                      formatted and used to restrict the range of HapMap loci" << endl 
       << "                      output." << endl
       << "-o   -outputfile      output case-control-file  " << endl
       << "                      valid only with -i, defaults to 'CaseControlGenotypes.txt'." << endl
       << "-m   -missing         missing-value character, defaults to 'N'" << endl
       << "-f   -flank           length in Kb of flanking region outside of typed region." << endl
       << "                      Valid only with -i. Defaults to 10." << endl << endl
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
unsigned FPHDOptions::getChrNum()const{
  return Chr;
}

unsigned long FPHDOptions::getNumUserLoci()const{
  return userloci;
}

bool FPHDOptions::WriteCCFile()const{
  return (incasecontrolfilename.size()>0);
}
const char* FPHDOptions::getInCCFilename()const{
  return incasecontrolfilename.c_str();
}
const char* FPHDOptions::getOutCCFilename()const{
  return outcasecontrolfilename.c_str();
}
float FPHDOptions::getFlankLength()const{
  return flankLength;
}
char FPHDOptions::getMissingChar()const{
  return MissingChar;
}
float FPHDOptions::getMinOverlap()const{
  return MinOverlap;
}
