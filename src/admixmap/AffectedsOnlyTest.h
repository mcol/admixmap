// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   AffectedsOnlyTest.h 
 *   header file for AffectedsOnly class
 *   Copyright (c) 2006, 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef AFFECTEDSONLYTEST_H
#define AFFECTEDSONLYTEST_H 1

#include "ScoreTestBase.h"
#include "bcppcl/RObjectWriter.h"

class Genome;
class LogWriter;

/**
   Class to implement affecteds-only score test for linkage with locus ancestry
 */
class AffectedsOnlyTest : public ScoreTestBase{

public:
  AffectedsOnlyTest();
  ~AffectedsOnlyTest();

  void Initialise(const char* filename, const int NumPopulations, const int NumLoci);

  void Reset();
  void Update(unsigned int locus, int k0, const double* const Theta, 
	      bool RandomMatingModel, bool diploid, const std::vector<std::vector<double> > AProbs);

  void Accumulate();

  void Output(const Vector_s& PopLabels, const Genome& Loci);
  void WriteFinalTable(const char* filename, const Vector_s& PopLabels, 
		       const Genome& Loci, LogWriter& Log);
  void OutputLikRatios(const char* const filename, const Vector_s& PopLabels, const Genome& Loci);

private:
  unsigned K, L;
  unsigned firstpoplabel;
  double* SumAffectedsScore2;
  double* SumAffectedsVarScore;
  double* SumAffectedsScore;
  double* SumAffectedsInfo;

  //score test objects, static so they can accumulate sums over individuals
  double *AffectedsScore;
  double *AffectedsVarScore;
  double *AffectedsInfo;

  double *LikRatio1;
  double *LikRatio2;

  RObjectWriter R;

  void OutputAffectedsOnlyTest(FileWriter& outfile, const Vector_s& PopLabels, 
			       const Genome& Loci, const std::string& sep, bool final);
};






#endif /* !defined SCORETESTS_H */
