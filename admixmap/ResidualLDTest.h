// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   ResidualLDTest.h 
 *   header file for ResidualLDTest class
 *   Copyright (c) 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef RESIDUALLDTEST_H
#define RESIDUALLDTEST_H 1

#include <sstream>
#include "ScoreTestBase.h"
#include "AdmixOptions.h"
#include "Genome.h"
#include "AlleleFreqs.h"

/**
   Class to implement score tests for residual allelic association between adjacent pairs of linked loci
 */
class ResidualLDTest : public ScoreTestBase{

public:
  ResidualLDTest();

  void Initialise(AdmixOptions* , const IndividualCollection* const, const Genome* const ,
		  LogWriter &);

  void Output(int iterations, bool final, const std::vector<std::string>& LocusLabels);
  void ROutput();

  void Update(double);
  void Update(const array_of_allelefreqs& Allelefreqs, bool ishapmixmodel);
  void Reset();

  ~ResidualLDTest();

private:
  std::vector<std::vector<std::vector<double> > > Score;
  std::vector<std::vector<std::vector<double> > > Info;
  std::vector<std::vector<std::vector<double> > > SumScore;
  std::vector<std::vector<std::vector<double> > > SumScore2;
  std::vector<std::vector<std::vector<double> > > SumInfo;

  const AdmixOptions *options;
  const IndividualCollection *individuals;
  const Genome* Lociptr;//Pointer to Loci
  const Chromosome* const* chrm;//Copy of pointer to array of chromosomes
  int rank, worker_rank, NumWorkers;
  
  //OUTPUT
  void OutputTestsForResidualAllelicAssociation(int iterations, ofstream* outputstream, bool final, 
						const std::vector<std::string>& LocusLabels);
  
  void UpdateScoresForResidualAllelicAssociation(int c, int locus, const double* const AlleleFreqsA, const double* const AlleleFreqsB);
  void UpdateScoresForResidualAllelicAssociation2(int c, int locus, const double* const AlleleFreqsA, const double* const AlleleFreqsB);

};

#endif /* !defined SCORETESTS_H */
