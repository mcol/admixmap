// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   ScoreTests.h 
 *   header file for Scoretests class
 *   Copyright (c) 2005-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef SCORETESTS_H
#define SCORETESTS_H 1

#include <sstream>
#include "IndividualCollection.h"
#include "bclib/LogWriter.h"
#include "common.h"
#include "AdmixtureAssocTest.h"
#include "AllelicAssocTest.h"
#include "AffectedsOnlyTest.h" 
#include "AncestryAssocTest.h"

class FreqArray;


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
