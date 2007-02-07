// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   AffectedsOnlyTest.h 
 *   header file for AffectedsOnly class
 *   Copyright (c) 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
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

class Genome;
/**
   Class to implement affecteds-only score test for linkage with locus ancestry
 */
class AffectedsOnlyTest : public ScoreTestBase{

public:
  AffectedsOnlyTest();
  ~AffectedsOnlyTest();

  void Initialise(const char* filename, const int NumPopulations, const int NumLoci, LogWriter &Log);

  void Reset();
  void Update(unsigned int locus, int k0, const double* const Theta, 
	      bool RandomMatingModel, bool diploid, const std::vector<std::vector<double> > AProbs);

  void Accumulate(unsigned j);

  void Output(int iterations, const Vector_s& PopLabels, const Genome& Loci, bool final = false, const char* filename = 0);
  void ROutput(const int numIterations);
  void OutputLikRatios(const char* const filename, int iterations, const Vector_s& PopLabels, const Genome& Loci);

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

};






#endif /* !defined SCORETESTS_H */
