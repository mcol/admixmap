// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   AdmixtureAssocTest.h 
 *   Copyright (c) 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef ADMIXTUREASSOCTEST_H
#define ADMIXTUREASSOCTEST_H 1

#include "ScoreTestBase.h"
#include <fstream>

class IndividualCollection;
namespace bclib{
  class LogWriter;
}

/**
   Class to implement score test for admixture association (admixturescoretest)
 */
class AdmixtureAssocTest : public ScoreTestBase{

public:
  AdmixtureAssocTest();
  ~AdmixtureAssocTest();

  void Initialise(const unsigned K, const unsigned NumOutcomes, const char* filename,
		  const Vector_s& PLabels, bclib::LogWriter &);

  void Reset();

  void Output( );

  void Accumulate();
  void UpdateIndividualScore( const double* const Theta, double YMinusEY,double phi, double DInvLink, bool RandomMatingModel);

private:
  double* Score; 
  double* Info; 
  double* SumScore; 
  double* SumScore2; 
  double* SumInfo;
  unsigned NumPopulations, NumOutcomeVars;
  std::ofstream outputfile;

  void InitialiseAssocScoreFile(const Vector_s&);
};

#endif /* !defined ADMIXTUREASSOCTEST_H */
