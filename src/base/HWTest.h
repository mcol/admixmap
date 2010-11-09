//=============================================================================
//
// Copyright (C) 2002-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file HWTest.h
/// Definition of the HWTest class.
//=============================================================================

#ifndef HWTEST_H
#define HWTEST_H 1


#include "ScoreTestBase.h"
#include "Genome.h"
#include "IndividualCollection.h"

namespace bclib{
  class LogWriter;
}



/** \addtogroup base
 * @{ */


class Options;

/**
  Class to implement a score test for deviation from Hardy-Weinberg equilibrium
  in order to test for genotyping errors.
  This version evaluates for each single locus and sums over individuals.
  Code for evaluation for each individual (rather than summing) is commented out.
*/
class HWTest : public ScoreTestBase{

public:
  HWTest();

  void Initialise(int nloci);

  void Output(const char *filename, const Vector_s& LocusLabels,
              bclib::LogWriter& Log);

  void Update(const IndividualCollection* const IC, const Genome* const Loci);

  ~HWTest();

private:
  int NumInd, NumLoci;
  //dmatrix score;
  //dmatrix sumscore;
  //dmatrix sumscore2;
  //dmatrix suminfo;
  double *score;
  double *sumscore;
  double *sumscore2;
  double *suminfo;

  void Reset();//required function

  /**
   *  UNIMPLEMENTED: to avoid undesired copying.
   */
  HWTest(const HWTest&);
  HWTest& operator=(const HWTest&);
};



/** @} */



#endif /* !defined HWTEST_H */
