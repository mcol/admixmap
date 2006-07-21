// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   DispersionTest.h 
 *   header file for DispersionTest class
 *   Copyright (c) 2005 LSHTM
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
#ifndef DISPTEST_H
#define DISPTEST_H 1

#include "IndividualCollection.h"
#include "Genome.h"
#include "AlleleFreqs.h"
#include "LogWriter.h"

///Class to implement test for dispersion of allele frequencies between unadmixed populations sampled and the corresponding ancestry-specific allele frequencies in the admixed population under study
class DispersionTest{
 public:
  DispersionTest();
  ~DispersionTest();
  void Initialise(const AdmixOptions* const, LogWriter &, int);
  void Output(int, const Genome &, const Vector_s& PopLabels);
  void TestForDivergentAlleleFrequencies(const AlleleFreqs* const, const IndividualCollection* const IC);

 private:
  int NumberOfCompositeLoci;
  int NumberOfPopulations;
  const AdmixOptions *options;
  std::ofstream dispersionoutputstream;
  int **divergentallelefreqstest;

};











#endif /* !defined DISPTEST_H */
