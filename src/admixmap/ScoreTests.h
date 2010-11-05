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
/// \file ScoreTests.h
/// Definition of the ScoreTests class.
//=============================================================================

#ifndef SCORETESTS_H
#define SCORETESTS_H 1

#include "IndividualCollection.h"
#include "common.h"
#include "AdmixtureAssocTest.h"
#include "AllelicAssocTest.h"
#include "AffectedsOnlyTest.h" 
#include "AncestryAssocTest.h"

class FreqArray;
namespace bclib {
  class LogWriter;
}


/** \addtogroup admixmap
 * @{ */


/**
 * Class acts as a container for the following tests:
 *   (1) Score test for allelic association
 *   (2) Score test for within-halpotype association
 *   (3) Score test for admixture association (admixturescoretest)
 *   (4) Score test for linkage with locus ancestry
 *   (5) Affecteds-only score test for linkage with locus ancestry
 */
class ScoreTests{

public:
  ScoreTests();

  void Initialise(AdmixOptions* , const IndividualCollection* const, const Genome* const ,
		  const Vector_s&, bclib::LogWriter &);

  void Output(const Vector_s& PLabels);
  void WriteFinalTables(const Vector_s& PLabels, bclib::LogWriter& Log);

  void MergeRareHaplotypes(const std::vector<double> &alpha0);
  void Update(const vector<bclib::Regression* >& R);
  void UpdateScoresForResidualAllelicAssociation(const FreqArray& Allelefreqs);

  ~ScoreTests();

  AffectedsOnlyTest& getAffectedsOnlyTest()  { return AffectedsOnlyScoreTest; }
  CopyNumberAssocTest& getAncestryAssocTest(){ return AncestryAssocScoreTest; }
  void OutputLikelihoodRatios(const Vector_s& PopLabels);

private:
  const AdmixOptions *options;
  const IndividualCollection *individuals;
  const Genome* Lociptr;//Pointer to Loci

  void Reset();

  AllelicAssocTest AllelicAssociationTest;
  AdmixtureAssocTest AdmixtureAssocScoreTest;
  AffectedsOnlyTest AffectedsOnlyScoreTest;
  AncestryAssocTest AncestryAssocScoreTest;

};


/** @} */


#endif /* !defined SCORETESTS_H */
