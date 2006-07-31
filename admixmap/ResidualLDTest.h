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
#ifdef PARALLEL
  void SetComm(const MPI::Intracomm* c, const std::vector<std::string>* locuslabels);
#endif

  void Output(int iterations, bool final);
  void ROutput();

  void Update(double);
  void Update(const array_of_allelefreqs& Allelefreqs);
  void Reset();

  ~ResidualLDTest();

private:
  double*** Score;
  double*** Info;
  double*** SumScore;
  double*** SumScore2;
  double*** SumInfo;

#ifdef PARALLEL
  double* sendresallelescore;
  double* recvresallelescore;
  double* sendresalleleinfo;
  double* recvresalleleinfo;
  int dimresallelescore, dimresalleleinfo;
  
  const std::vector<std::string> * LocusLabels;
  const MPI::Intracomm * Comm;//pointer to the workers_and_master communicator in admixmap.cc
#endif
  
  const AdmixOptions *options;
  const IndividualCollection *individuals;
  const Genome* Lociptr;//Pointer to Loci
  const Chromosome* const* chrm;//Copy of pointer to array of chromosomes
  int rank, worker_rank, NumWorkers;
  
  //OUTPUT
  void OutputTestsForResidualAllelicAssociation(int iterations, ofstream* outputstream, bool final);
  
  void UpdateScoresForResidualAllelicAssociation(int c, int locus, 
						 const double* const AlleleFreqsA, 
						 const double* const AlleleFreqsB);
};

#endif /* !defined SCORETESTS_H */
