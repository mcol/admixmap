//=============================================================================
//
// Copyright (C) 2005-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file DispersionTest.h
/// Definition of the DispersionTest class.
//=============================================================================


#ifndef DISPTEST_H
#define DISPTEST_H 1

#include "IndividualCollection.h"
#include "Genome.h"
#include "AlleleFreqs.h"
namespace bclib{
  class LogWriter;
}


/** \addtogroup admixmap
 * @{ */


/// Class to implement test for dispersion of allele frequencies between
/// unadmixed populations sampled and the corresponding ancestry-specific
/// allele frequencies in the admixed population under study.
/// This is evaluated for each subpopulation at each locus, and as a global
/// test over all loci. Valid only if option priorallelefreqfile is specified.
/// The results are "Bayesian p-values".
class DispersionTest{
public:
  DispersionTest();
  ~DispersionTest();
  void Initialise(const std::string& resultsDir, bclib::LogWriter &Log, int NumLoci, int NumPopulations);
  void Output(int, const Genome &, const Vector_s& PopLabels);
  void TestForDivergentAlleleFrequencies(const AlleleFreqs* const, const IndividualCollection* const IC);
  
private:
  int NumberOfCompositeLoci;
  int NumberOfPopulations;
  std::ofstream dispersionoutputstream;
  int **divergentallelefreqstest;

};


/** @} */


#endif /* !defined DISPTEST_H */
