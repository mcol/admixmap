// *-*-C++-*-*
/**
  \file FPHDOptions.h
  Options for FPHD.
  This file is part of FPHD

  This program is free software distributed WITHOUT ANY WARRANTY. 
  You can redistribute it and/or modify it under the terms of the GNU General Public License, 
  version 2 or later, as published by the Free Software Foundation. 
  See the file COPYING for details.

  Copyright (c) David O'Donnell 2007
*/

#include <fstream>
#include <string>
#include "bclib/OptionReader.h"

class FPHDOptions {
public:
  FPHDOptions(int argc, char** argv);
  ~FPHDOptions(){};
  static void PrintHelpText();
  bool Verbose()const;
  bool LimitedLoci()const;
  bool Backup()const;
  const std::string& getPrefix()const;
  //unsigned getChrNum()const;
  unsigned long getMaxLoci()const;
  bool WriteObsGenoFile()const;
  const char* getLocusFilename()const;
  const char* getGenotypesFilename()const;
  const char* getInObsGenoFilename()const;
  const char* getOutObsGenoFilename()const;
  const char* getInitialMixturePropsFilename()const;
  const char* getInitialArrivalRateFilename()const;
  const char* getInitialAlleleFreqFilename()const;
  const char* getInitialFreqPriorFilename()const;
  float getFlankLength()const;
  char getMissingChar()const;
  unsigned getMinOverlap()const;
  void setMaxLoci(unsigned);
private:
  bool beVerbose, backup;
  std::string prefix;
  //bool LimitLoci;
  //  string Chr;//chromosome number

  std::string genotypesfilename, locusfilename;
  std::string inobsgenofilename;
  std::string outobsgenofilename;
  unsigned long MaxLoci;//max loci per subchromosome
  unsigned MinOverlap_bp;//minimum overlap between subchromosomes
  float MinOverlap_kb;
  float flankLength;//length in Kb of flanking region if using CCgenotypesfile 
  char MissingChar;

  //initial value files
  std::string InitialMixturePropsFilename;
  std::string InitialArrivalRateFilename;
  std::string InitialAlleleFreqFilename;
  std::string InitialFreqPriorFilename;

  FPHDOptions();
  void DefineOptions(bclib::OptionReader& opt);
  void ParseOptions(int argc, char** argv);
  void StripSuffix(std::string& filename);
};
