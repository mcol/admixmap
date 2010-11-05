//=============================================================================
//
// Copyright (C) 2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file AllelicAssocTest.h
/// Definition of the AllelicAssocTest class.
//=============================================================================

#ifndef ALLELICASSOCTEST_H
#define ALLELICASSOCTEST_H 1

#include "ScoreTestBase.h"
#include "AdmixOptions.h"
#include "IndividualCollection.h"
#include "Genome.h"
#include "AllelicAssocSubTest.h"
#include "bclib/RObjectWriter.h"


/** \addtogroup admixmap
 * @{ */


/**
 *  Class to implement 
 *   (1) Score test for allelic association
 *   (2) Score test for within-haplotype association
 */
class AllelicAssocTest : public ScoreTestBase{

public:
  AllelicAssocTest();
  ~AllelicAssocTest();

  void Initialise( AdmixOptions* , const IndividualCollection* const, const Genome* const ,
		  bclib::LogWriter &);

  void Reset();

  void Output();
  void WriteFinalTables(bclib::LogWriter& Log);

  void MergeRareHaplotypes(const std::vector<double> &alpha0);

  void Update( const PedBase * , double, double, double, bool);

  void Accumulate();


private:
  std::vector<AllelicAssocSubTest*> SubTests;

  std::vector<HaplotypeTest*> HaplotypeAssocTests; 

  int *locusObsIndicator;

  bclib::RObjectWriter AllelicAssocRObject;
  bclib::RObjectWriter HaplotypeAssocRObject;
  //outputfile = allelicAssocStream
  //std::ofstream allelicAssocScoreStream;

  const AdmixOptions *options;
  const IndividualCollection *individuals;
  const Genome* Lociptr;
  unsigned NumCompositeLoci;
  std::vector<unsigned> NumLoci;//number of loci within each comp locus
  std::vector<unsigned> NumMergedHaplotypes;
  const Chromosome* const* chrm;//Copy of pointer to array of chromosomes
  unsigned NumOutputs;//counts calls to output function for dimensions of R objects
  void ROutput();
};


/** @} */


#endif /* !defined ALLELICASSOCTEST_H */
