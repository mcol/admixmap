// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   AdmixOptions.h 
 *   header file for MisSpecAlleleFreqTest class
 *   Copyright (c) 2005, 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef MISSPECALLELEFREQTEST_H
#define MISSPECALLELEFREQTEST_H 1

#include "common.h"
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "IndividualCollection.h"
#include "Individual.h"
#include "AlleleFreqs.h"
#include "Genome.h"
#include "AdmixOptions.h"
#include "LogWriter.h"
#include "ScoreTestBase.h"

///a scalar test for SNPs only
class MisSpecifiedAlleleFreqTest : public ScoreTestBase{
public:
  MisSpecifiedAlleleFreqTest();
  ~MisSpecifiedAlleleFreqTest();

  void Initialise(const AdmixOptions* const options, const Genome* const Loci, LogWriter &Log );
  void Update(const IndividualCollection* const individuals, const AlleleFreqs* const A, const Genome* const Loci);
  void Output(int iteration, const Genome* const Loci, const Vector_s& PopLabels);
  void Reset();
  
private:
  double **Score;
  double **Info;
  double **SumScore;
  double **SumScoreSq;
  double **SumInfo;
  int NumTestLoci;     //number of comp loci with a single locus ie those used in scalar score test for misspecified allelefreqs
   int NumCompLoci; //number of composite loci
  int Populations;           //number of populations
  
  void UpdateLocus(int j, const double* const* phi, int NumCopiesAllele1,
		   const double* const AlleleFreqs, int NumStates);
  
};

///a vector test, not fully tested
class MisSpecifiedAlleleFreqTest2: public ScoreTestBase{
public:
  MisSpecifiedAlleleFreqTest2();
  ~MisSpecifiedAlleleFreqTest2();
  
  void Initialise(const AdmixOptions* const options, const Genome* const Loci, LogWriter &Log );
  void Update(const IndividualCollection* const individuals, const AlleleFreqs* const A, const Genome* const Loci);
  void Output(int iteration, const Genome* const Loci, const Vector_s& PopLabels);
  void Reset(){};

private:
  int NumCompLoci; //number of composite loci
  int Populations;           //number of populations

  double ***SumScore;
  double ***SumScoreSq;
  double ***SumInfo;

  void UpdateScoreForMisSpecOfAlleleFreqs2(const int j, const int NumberOfStates, const double* const AlleleFreqs, 
					   const IndividualCollection* const individuals);

  void OutputTestsForMisSpecifiedAlleleFreqs2( int samples, const Genome* const Loci, const std::string* const PopLabels);

};

/// a container class for the tests for misspecified allele frequencies
class MisSpecAlleleFreqTest{

public:
  MisSpecAlleleFreqTest();
  ~MisSpecAlleleFreqTest();

  void Initialise(const AdmixOptions* const options, const Genome* const Loci, LogWriter &Log );
  void Update(const IndividualCollection* const individuals, const AlleleFreqs* const A, const Genome* const Loci);
  void Output(int iteration, const Genome* const Loci, const Vector_s& PopLabels);

private:
  bool doTest1, doTest2;//indicators for the two tests

  MisSpecifiedAlleleFreqTest Test1;
  MisSpecifiedAlleleFreqTest2 Test2;

};

#endif

