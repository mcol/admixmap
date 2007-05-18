// *-*-C++-*-*
/*
  FPHDOptions.h
  Options for FPHD
  This file is part of FPHD

  This program is free software distributed WITHOUT ANY WARRANTY. 
  You can redistribute it and/or modify it under the terms of the GNU General Public License, 
  version 2 or later, as published by the Free Software Foundation. 
  See the file COPYING for details.

  Copyright (c) David O'Donnell 2007
*/

#include <fstream>
#include <string>

class FPHDOptions{
 public:
  FPHDOptions(int argc, char*const* argv, std::ofstream&, std::ofstream&);
  ~FPHDOptions(){};
  static void PrintHelpText();
  bool Verbose()const;
  bool LimitedLoci()const;
  const std::string& getPrefix()const;
  unsigned getChrNum()const;
  unsigned long getNumUserLoci()const;
  bool WriteCCFile()const;
  const char* getInCCFilename()const;
  const char* getOutCCFilename()const;
  float getFlankLength()const;
  char getMissingChar()const;

 private:
  bool beVerbose;
  std::string prefix;
  bool LimitLoci;
  unsigned Chr;//chromosome number

  char* incasecontrolfilename;
  char* outcasecontrolfilename;
  unsigned long userloci;
  float flankLength;//length in Kb of flanking region if using CCgenotypesfile 
  char MissingChar;

  FPHDOptions();
  void ParseOptions(int argc, char*const* argv, std::ofstream&, std::ofstream&);
};
