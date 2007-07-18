// *-*-C++-*-*
/** 
 *   HAMIXMAP
 *   HapMixAllelicAssocTest.h 
 *   wrapper for the hapmix model allelic assoc test
 *   Extension of the CopyNumberAssocTest   
 *   
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef HAPMIXALLELICASSOCTEST_H
#define HAPMIXALLELICASSOCTEST_H

#include "CopyNumberAssocTest.h"
class InputHapMixData;
class HapMixIndividualCollection;
class Regression;
class Genome;

namespace bclib{
  class Regression;
  class LogWriter;
}

///test using conditional distribution of copies of allele2, rather than sampled values
class HapMixAllelicAssocTest : public CopyNumberAssocTest{
public:
  HapMixAllelicAssocTest();
  ~HapMixAllelicAssocTest();
  void Initialise(const std::string& ResultsDir, const int NumLoci);
  void Update(const HapMixIndividualCollection* const IC, const bclib::Regression* const R, const Genome& Loci);
  void Output(const Genome& Loci);
  void WriteFinalTable(const std::string& ResultsDir, const Genome& Loci, const InputHapMixData& data, bclib::LogWriter& Log);

private:
  void PrintAverageInfo(bclib::LogWriter& Log, const InputHapMixData& data, const char* filename);
  //override base functions, making them private
  void Initialise(const char* filename, const int NumPopulations, const int NumLoci);
  void Update(int locus, const double* Covariates, double phi, double YMinusEY, double DInvLink, 
	      const std::vector<std::vector<double> > Probs) ;
};

#endif

