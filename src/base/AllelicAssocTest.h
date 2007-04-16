// *-*-C++-*-*
/** 
 *   AllelicAssocTest.h 
 *   header file for AllelicAssocTest class
 *   Copyright (c) 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef ALLELICASSOCTEST_H
#define ALLELICASSOCTEST_H 1

#include "ScoreTestBase.h"
#include "Options.h"
#include "IndividualCollection.h"
#include "Genome.h"

#include "AllelicAssocSubTest.h"

/**
 *  Class to implement 
 *   (1) Score test for allelic association
 *   (2) Score test for within-haplotype association
 */
class AllelicAssocTest : public ScoreTestBase{

public:
  AllelicAssocTest();
  ~AllelicAssocTest();

  void Initialise( Options* , const IndividualCollection* const, const Genome* const ,
		  LogWriter &);

  void Reset();

  void Output(const Vector_s& LocusLabels, bool final);

  void ROutput();

  void SetAllelicAssociationTest(const std::vector<double> &alpha0);

  void Update( const Individual* const , double, double, double, bool);

  void Accumulate();


private:
  std::vector<AllelicAssocSubTest*> SubTests;

  std::vector<HaplotypeTest*> HaplotypeAssocTests; 

  bool onFirstLineAllelicAssoc;

  int *locusObsIndicator;

  bool onFirstLineHapAssoc;

  std::ofstream HaplotypeAssocScoreStream;
  //outputfile = allelicAssocStream
  //std::ofstream allelicAssocScoreStream;

  const Options *options;
  const IndividualCollection *individuals;
  const Genome* Lociptr;//Pointer to Loci
  const Chromosome* const* chrm;//Copy of pointer to array of chromosomes
  int rank, worker_rank, NumWorkers;
  unsigned NumOutputs;//counts calls to output function for dimensions of R objects

};

#endif /* !defined ALLELICASSOCTEST_H */
