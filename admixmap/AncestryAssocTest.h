// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   AncestryAssocTest.h 
 *   header file for AncestryAssocTest class
 *   Copyright (c) 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef ANCESTRYASSOCTEST_H
#define ANCESTRYASSOCTEST_H 1

#include "ScoreTestBase.h"

class Genome;
/**
   Class to implement score test for linkage with locus ancestry
 */
class AncestryAssocTest : public ScoreTestBase{

public:
  AncestryAssocTest();
  ~AncestryAssocTest();
  void Initialise(const char* filename, const int NumPopulations, const int NumLoci, LogWriter &Log);

  void Reset();
  void Update(int locus, const double* admixtureCovars, double phi, double YMinusEY, double DInvLink, 
	      const std::vector<std::vector<double> > AProbs) ;
  void UpdateB(double DInvLink, double dispersion, const double* admixtureCovars);

  void Accumulate(unsigned j);

  void Output(int iterations, const Vector_s& PopLabels, const Genome& Loci, bool final = false, const char* filename = 0);
  void ROutput(const int numIterations);

private:
  unsigned K, L;
  unsigned firstpoplabel;

  double **AncestryScore;
  double **AncestryInfo;
  double **AncestryVarScore;
  double **AncestryInfoCorrection;
  double *B;//used for ancestry score test
  double *PrevB;//holds B for previous iteration while B accumulates for this iteration
  double *Xcov; //column matrix of covariates used to calculate B and for score test, 
                       //static only for convenience since it is reused each time
  double* SumAncestryScore;
  double* SumAncestryInfo;
  double* SumAncestryVarScore;
  double* SumAncestryScore2;
};






#endif /* !defined ANCESTRYASSOCTEST_H */
