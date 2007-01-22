// *-*-C++-*-*
/** 
 *   HWTest.h 
 *   header file for HWTest class
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef HWTEST_H
#define HWTEST_H 1

#include "ScoreTestBase.h"
#include "common.h"
#include "AdmixOptions.h"
#include "Genome.h"
#include "IndividualCollection.h"

class LogWriter;
//class IndividualCollection;
/**
  Class to implement a score test for deviation from Hardy-Weinberg equilibrium
  in order to test for genotyping errors.
  This version evaluates for each single locus and sums over individuals.
  Code for evaluation for each individual (rather than summing) is commented out.
*/
class HWTest : public ScoreTestBase{

public:
  HWTest();

  void Initialise(const AdmixOptions* const options, int nloci, LogWriter &Log);

  void Output(const Vector_s locuslabels);

  void Update(const IndividualCollection* const IC, const Genome* const Loci);

  ~HWTest();

private:
  int NumInd, NumLoci, samples;
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
  void operator=(const HWTest&);
};






#endif /* !defined HWTEST_H */
