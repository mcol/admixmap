//=============================================================================
//
// Copyright (C) 2006, 2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file ResidualLDTest.h
/// Definition of the ResidualLDTest class.
//=============================================================================

#ifndef RESIDUALLDTEST_H
#define RESIDUALLDTEST_H 1


#include "ScoreTestBase.h"
#include "Options.h"
#include "bclib/RObjectWriter.h"


/** \addtogroup base
 * @{ */


class Chromosome;
class Genome;
class IndividualCollection;
class FreqArray;

/**
   Class to implement score tests for residual allelic association between adjacent pairs of linked loci
 */
class ResidualLDTest : public ScoreTestBase{

public:
  ResidualLDTest();

  void Initialise(const char* filename , const IndividualCollection* const, const Genome* const);

  void Output(const std::vector<std::string>& LocusLabels);
  void WriteFinalTable(const char* filename, const std::vector<std::string>& LocusLabels, bclib::LogWriter& Log);

  void Update(double);
  void Update(const FreqArray& Allelefreqs, bool ishapmixmodel);
  void Reset();

  ~ResidualLDTest();

private:
  bclib::RObjectWriter R;
  std::vector<std::vector<std::vector<double> > > Score;
  std::vector<std::vector<std::vector<double> > > Info;
  std::vector<std::vector<std::vector<double> > > SumScore;
  std::vector<std::vector<std::vector<double> > > SumScore2;
  std::vector<std::vector<std::vector<double> > > SumInfo;

  const IndividualCollection *individuals;
  const Genome* Lociptr;//Pointer to Loci
  const Chromosome* const* chrm;//Copy of pointer to array of chromosomes
  //std::vector<unsigned> Tcount;
  unsigned NumIntervals;
  
  //OUTPUT
  void OutputTestsForResidualAllelicAssociation(bclib::DelimitedFileWriter& outputstream, bool final, 
						const std::vector<std::string>& LocusLabels);
  
  void UpdateScoresForResidualAllelicAssociation(int c, int locus, const double* const AlleleFreqsA, const double* const AlleleFreqsB);
  void UpdateScoresForResidualAllelicAssociation2(int c, int locus, const double* const AlleleFreqsA, const double* const AlleleFreqsB);

};


/** @} */

#endif /* !defined SCORETESTS_H */
