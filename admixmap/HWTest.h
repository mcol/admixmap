// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   HWTest.h 
 *   header file for HWTest class
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
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

#include "common.h"
#include "AdmixOptions.h"
#include "Genome.h"
#include "IndividualCollection.h"

class LogWriter;

class HWTest{

public:
  HWTest();

  void Initialise(AdmixOptions *options, int nloci, LogWriter *Log);

  void Output(Matrix_s locusdata);

  void Update(IndividualCollection *IC, Chromosome **C, Genome *Loci);

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

  std::ofstream outputfile;

  string double2R( double x );
  /**
   *  UNIMPLEMENTED: to avoid undesired copying.
   */    
  HWTest(const HWTest&);
  void operator=(const HWTest&);
};






#endif /* !defined HWTEST_H */
