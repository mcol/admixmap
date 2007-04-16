// *-*-C++-*-*
/** 
 *   ScoreTests.h 
 *   header file for Scoretests class
 *   Copyright (c) 2005, 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
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
#include "bcppcl/LogWriter.h"
#include "common.h"
#include "AdmixtureAssocTest.h"
#include "AllelicAssocTest.h"
#include "AffectedsOnlyTest.h" 
#include "CopyNumberAssocTest.h"
#include "ResidualLDTest.h"

class FreqArray;

/**
   Class to implement various score tests. 
 *   Class implements the following score tests:

 *   (1) Score test for allelic association
 *   (2) Score test for within-halpotype association
 *
 * Also acts as a wrapper for the following tests:
 *   (3) Score test for admixture association (admixturescoretest)
 *   (4) Score test for linkage with locus ancestry
 *   (5) Affecteds-only score test for linkage with locus ancestry
 *   (6) Score test for residual allelic association between adjacent pairs of linked loci
 *   (7) hapmix allelic assoc test (using CopyNumberAssocTest)
 */
class ScoreTests{

public:
  ScoreTests();

  void Initialise(Options* , const IndividualCollection* const, const Genome* const ,
		  const Vector_s&, LogWriter &);

  void Output(int iterations, const Vector_s& PLabels, const Vector_s& LocusLabels, bool final);

  void ROutput();
  void SetAllelicAssociationTest(const std::vector<double> &alpha0);
  void Update(const vector<Regression* >& R);
  void UpdateScoresForResidualAllelicAssociation(const FreqArray& Allelefreqs);

  ~ScoreTests();

  AffectedsOnlyTest& getAffectedsOnlyTest();
  CopyNumberAssocTest& getAncestryAssocTest();
  void OutputLikelihoodRatios(const char* const filename, int iterations, const Vector_s& PopLabels);

private:
  const Options *options;
  const IndividualCollection *individuals;
  const Genome* Lociptr;//Pointer to Loci
  const Chromosome* const* chrm;//Copy of pointer to array of chromosomes
  int rank, worker_rank, NumWorkers;

//OUTPUT

  void Reset();

  AllelicAssocTest AllelicAssociationTest;
  ResidualLDTest ResidualAllelicAssocScoreTest;//here temporarily, until rest of scoretests classes have been created
  AdmixtureAssocTest AdmixtureAssocScoreTest;
  AffectedsOnlyTest AffectedsOnlyScoreTest;
  CopyNumberAssocTest AncestryAssocScoreTest;

  CopyNumberAssocTest HapMixAllelicAssocTest;//test using conditional distribution of copies of allele2, rather than sampled values
};






#endif /* !defined SCORETESTS_H */
