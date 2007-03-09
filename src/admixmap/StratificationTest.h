// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   StratificationTest.h 
 *   Class to implement a test for residual population stratification
 *   Copyright (c) 2005, 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef STRATIFICATIONTEST
#define STRATIFICATIONTEST 1

#include <vector>
#include <iostream>
#include "IndividualCollection.h"
#include "Genome.h"
#include "utils/LogWriter.h"

class FreqArray;

///Class to implement a test for residual population stratification
class StratificationTest
{
public:
  StratificationTest();
  
  void Initialize( AdmixOptions* const options, const Genome &Loci,  
		   const IndividualCollection* const IC, LogWriter &Log);
  void OpenOutputFile( const char * , LogWriter &);

  void calculate( const IndividualCollection* const individuals, const FreqArray& AlleleFreqs,
		  const std::vector<std::vector<int> > ChrmAndLocus, int Populations );

  void Output(LogWriter &);

private:
  int T;
  int count;
  int NumberOfTestLoci;
  std::vector<unsigned int> TestLoci;
  bool ModelIndicator;
  std::ofstream outputstream;

   std::vector<double>
  GenerateExpectedGenotype( const Individual* const, const double*, const int  );

  std::vector<unsigned short>
  SimGenotypeConditionalOnAdmixture( const std::vector<double> ProbAllele1 );

  std::vector<unsigned short>
  SimGenotypeConditionalOnAncestry( const double*, const int ancestry[2] );

  std::vector<unsigned short> SampleHeterozygotePhase( const double*, const int ancestry[2] );
};

#endif /* !defined STRATIFICATIONTEST_H */
