//=============================================================================
//
// Copyright (C) 2006  David O'Donnell, Clive Hoggart and Paul McKeigue
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file AdmixtureAssocTest.h
/// Definition of the AdmixtureAssocTest class.
//=============================================================================

#ifndef ADMIXTUREASSOCTEST_H
#define ADMIXTUREASSOCTEST_H 1

#include "ScoreTestBase.h"
#include "common.h"   // for Vector_s
#include <fstream>

class IndividualCollection;
namespace bclib{
  class LogWriter;
}


/** \addtogroup admixmap
 * @{ */


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


/** @} */


#endif /* !defined ADMIXTUREASSOCTEST_H */
