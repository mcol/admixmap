/**
  \file FPHDOptions.cpp
  Options for FPHD.
  This file is part of FPHD

  This program is free software distributed WITHOUT ANY WARRANTY. 
  You can redistribute it and/or modify it under the terms of the GNU General Public License, 
  version 2 or later, as published by the Free Software Foundation. 
  See the file COPYING for details.

  Copyright (c) David O'Donnell 2007
*/
#include "FPHDOptions.h"
#include <iostream>

#define MAXCHROMOSOMES 22//maximum number of chromosomes
#define PROGNAME "FPHD"

using namespace::std;

FPHDOptions::FPHDOptions(int argc, char** argv){
  //set defaults
  beVerbose = false;
  backup = true;
  //prefix = ".";
  inobsgenofilename;
  outobsgenofilename;
  MaxLoci = 1000000000;  //some number > number of HapMap loci
  //LimitLoci = false;
  flankLength = 10;//10Kb
  MinOverlap_kb = 5000;//5Mb
  MinOverlap_bp = (unsigned) (MinOverlap_kb * 1000.0);
  MissingChar = 'N';

  if(argc == 1){//no args specified - print usage and exit
    PrintHelpText();
    exit(1);
  }
  ParseOptions(argc, argv);

}

void FPHDOptions::ParseOptions(int argc, char** argv){
  bclib::OptionReader opt;
  opt.setVerbose(true); /* print warnings about unknown options */

  opt.setUserOption("genotypesfile", "genotypes.txt");
  opt.setUserOption("locusfile", "loci.txt");
  DefineOptions(opt);
  if(!opt.ReadCommandLineArgs(argc, argv))
    exit(1);
  if(opt.getFlag( "help" ) || opt.getFlag( 'h' )){
    PrintHelpText();
    exit(1);
  }
  if(!opt.SetOptions() | !opt.CheckRequiredOptions())
    exit(1);

  //check for upper limit on loci
 
  //convert flankLength in kb to bp
  flankLength *= 1000;

  //set verbose flag
  beVerbose = opt.getFlag("verbose");

  //set backup flag
  backup = opt.getFlag("backup");

  //set default ouput testgenotypesfilename if none specified
  if(inobsgenofilename.size() && !outobsgenofilename.size())
    outobsgenofilename = "ObservedGenotypes.txt";

  //set min overlap when dividing chromosomes
  MinOverlap_bp = (unsigned) (MinOverlap_kb * 1000.0);

  if(beVerbose){
    cout << "Writing genotypes to " << genotypesfilename << endl;
    cout << "Writing locusfile to " << locusfilename << endl;
    if(inobsgenofilename.size())
      cout << "Writing observed genotypes to " << outobsgenofilename << endl; 
  }

  StripSuffix(locusfilename);
  StripSuffix(genotypesfilename);
  StripSuffix(outobsgenofilename);
}

void FPHDOptions::StripSuffix(std::string& filename){
  if(filename.size() ==0) return;
  filename.erase(filename.find_first_of("."));
}

void FPHDOptions::DefineOptions(bclib::OptionReader& opt){
  opt.addFlag('h', "help");
  opt.addFlag('v', "verbose");
  opt.addOption('c', "chromosome", bclib::stringOption, &ChrLabel);
  opt.addOption('p', "prefix", bclib::stringOption, &prefix, true);
  opt.addOption('g', "genotypesfile", bclib::stringOption, &genotypesfilename);
  opt.addOption('l', "locusfile", bclib::stringOption, &locusfilename);
  opt.addFlag('b', "backup");
  //  opt.addOption('n', "numloci", bclib::longOption, &locuslimit);
  opt.addOption('M', "maxloci", bclib::intOption, &MaxLoci);
  opt.addOption('O', "minoverlap", bclib::floatOption, &MinOverlap_kb);
  opt.addOption('i', "inputfile", bclib::stringOption, &inobsgenofilename);
  opt.addOption('o', "outputfile", bclib::stringOption, &outobsgenofilename);
  opt.addOption('f', "flank", bclib::floatOption, &flankLength);
  opt.addOption('m', "missing", bclib::charOption, &MissingChar);
  opt.addOption("initialmixturepropsfile", bclib::stringOption, &InitialMixturePropsFilename);
  opt.addOption("initialarrivalratefile", bclib::stringOption, &InitialArrivalRateFilename);
  opt.addOption("initialallelefreqfile", bclib::stringOption, &InitialAlleleFreqFilename);
  opt.addOption("initialfreqpriorfile", bclib::stringOption, &InitialFreqPriorFilename);
}


void FPHDOptions::PrintHelpText(){
  //       |<-             This is 80 characters' width                                 ->|
  cout << "----------------------------------------------------------------------------"
       << endl
       << "This program prepares data for use with HAPMIXMAP." << endl
       << "It formats phased HapMap data as well as user-supplied case-control genotypes" << endl
       << "and prepares initial value files, breaking into separate files as necessary." << endl
       << "Copyright (c) David O'Donnell 2007" << endl
       << "All parts of this program are freely distributable" << endl << endl
       << "Usage: " << PROGNAME << "[option=value ...]" << endl << endl
       << "Options ('-' signs are optional), = denotes a default: " << endl
       << "-h   -help                  Print this help message and exit" << endl
       << "-c   -chrlabel              label to put in 4th col of locusfile" << endl
       << "-v   -verbose               Be verbose" << endl
       << "-p   -prefix = .            Prefix where HapMap files are located." << endl
       << "     Note that HapMap input files must be named '[prefix]_phased.txt', " << endl
       << "     '[prefix]_sample.txt' and '[prefix]_legend.txt'"
       << endl 
       << "-l   -locusfile = loci      Output locus file prefix" << endl
       << "-g   -genotypesfile = genotypes.txt" << endl
       << "                            Output genotypes file prefix" << endl
       << "-b   -backup                Back up original data if monomorphic loci found" << endl
       << "-i   -inputfile             Raw observed genotypes file. The genotypes will be" << endl
       << "                            encoded consistently with the HapMap genotypes and" << endl 
       << "                            the range of HapMap loci output will be restricted." << endl 
       << "-o   -outputfile = ObservedGenotypes.txt" << endl
       << "                            output testgenotypesfile prefix (valid only with -i)" << endl
       << "-M   -maxloci               Maximum number of loci per sub-chromosome/file." << endl
       << "                            If not specified, all available loci will be used" << endl
       << "-O   -minoverlap = 5000     Minimum overlap between sub-chromosomes in kb." << endl
       << "-m   -missing = N           missing-value character" << endl
       << "-f   -flank = 10            length in Kb of flanking region outside of" << endl
       << "                            typed region. Valid only with -i." << endl << endl
       << endl
       << "Initial Value Files - these will be broken up like the other output files" << endl
       << "     -initialallelefreqfile" << endl
       << "                            File with initial values of allele frequencies" << endl
       << "     -initialfreqpriorfile" << endl
       << "                            File with initial values of allele frequency prior" << endl
       << "                            parameters." << endl
       << "     -initialarrivalratefile" << endl
       << "                            File with initial values of arrival rates." << endl
       << "     -initialmixturepropsfile" << endl
       << "                            File with initial values of mixture proportions." << endl
       << endl << endl;
}

bool FPHDOptions::Verbose()const{
  return beVerbose;
}
// bool FPHDOptions::LimitedLoci()const{
//   return LimitLoci;
// }

const std::string& FPHDOptions::getPrefix()const{
  return prefix;
}
const string& FPHDOptions::getChrLabel()const{
  return ChrLabel;
}

bool FPHDOptions::Backup()const{
  return backup;
}

unsigned long FPHDOptions::getMaxLoci()const{
  return MaxLoci;
}

bool FPHDOptions::WriteObsGenoFile()const{
  return (inobsgenofilename.size()>0);
}
const char* FPHDOptions::getLocusFilename()const{
  return locusfilename.c_str();
}
const char* FPHDOptions::getGenotypesFilename()const{
  return genotypesfilename.c_str();
}
const char* FPHDOptions::getInObsGenoFilename()const{
  return inobsgenofilename.c_str();
}
const char* FPHDOptions::getOutObsGenoFilename()const{
  return outobsgenofilename.c_str();
}
float FPHDOptions::getFlankLength()const{
  return flankLength;
}
char FPHDOptions::getMissingChar()const{
  return MissingChar;
}
unsigned FPHDOptions::getMinOverlap()const{
  return MinOverlap_bp;
}
void FPHDOptions::setMaxLoci(unsigned L){
  MaxLoci = L < MaxLoci ? L : MaxLoci;
}
const char* FPHDOptions::getInitialMixturePropsFilename()const{
  return InitialMixturePropsFilename.c_str();
}
const char* FPHDOptions::getInitialArrivalRateFilename()const{
  return InitialArrivalRateFilename.c_str();
}
const char* FPHDOptions::getInitialAlleleFreqFilename()const{
  return InitialAlleleFreqFilename.c_str();
}
const char* FPHDOptions::getInitialFreqPriorFilename()const{
  return InitialFreqPriorFilename.c_str();
}
