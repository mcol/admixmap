// *-*-C++-*-*
/** 
 *   CopyNumberAssocTest.h 
 *   header file for CopyNumberAssocTest class
 *   Copyright (c) 2006, 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef COPYNUMBERASSOCTEST_H
#define COPYNUMBERASSOCTEST_H 1

#include "ScoreTestBase.h"

class Genome;

/**
   Class originally designed to implement score test for linkage with locus ancestry
 */
class CopyNumberAssocTest : public ScoreTestBase{

public:
  CopyNumberAssocTest();
  virtual ~CopyNumberAssocTest();
  virtual void Initialise(const char* filename, const int NumPopulations, const int NumLoci, LogWriter &Log);

  void Reset();
  virtual void Update(int locus, const double* Covariates, double phi, double YMinusEY, double DInvLink, 
	      const std::vector<std::vector<double> > Probs) ;
  void UpdateB(double DInvLink, double dispersion, const double* Covariates);

  void Accumulate();

  void Output(const Vector_s& PopLabels, const Genome& Loci, bool final = false, const char* filename = 0);
  virtual void ROutput() = 0;

protected:
  virtual void OpenOutputFile(LogWriter &Log, const char* filename) = 0;
  bool useprevb;
  unsigned K, L;
private:
  unsigned firstpoplabel;

  double **Score;
  double **Info;
  double **VarScore;
  double **InfoCorrection;
  double *B;//used to correct info
  double *PrevB;//holds B for previous iteration while B accumulates for this iteration
  double *Xcov; //column matrix of covariates used to calculate B and for score test, 
                       //static only for convenience since it is reused each time
  double* SumScore;
  double* SumInfo;
  double* SumVarScore;
  double* SumScore2;

  void Accumulate(unsigned j);
};






#endif /* !defined COPYNUMBERASSOCTEST_H */
