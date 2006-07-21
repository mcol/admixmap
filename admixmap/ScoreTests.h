// *-*-C++-*-*
/** 
 *   ADMIXMAP
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
#include "LogWriter.h"
#include "common.h"
#include "ResidualLDTest.h"

/**
   Class to implement various score tests. 
 *   Class implements the following score tests:
 *   (1) Score test for admixture association (admixturescoretest)
 *   (2) Score test for allelic association
 *   (3) Score test for within-halpotype association
 *   (4) Score test for linkage with locus ancestry
 *   (5) Affecteds-only score test for linkage with locus ancestry
 *   (6) Score test for residual allelic association between adjacent pairs of linked loci
 */
class ScoreTests{

public:
  ScoreTests();

  void Initialise(AdmixOptions* , const IndividualCollection* const, const Genome* const ,
		  const Vector_s&, LogWriter &);
#ifdef PARALLEL
  void SetComm(const MPI::Intracomm* c, const std::vector<std::string>* locuslabels);
#endif

  void InitialiseAssocScoreFile(const Vector_s&);

  void Output(int iterations, const Vector_s& PLabels, bool final);

  void ROutput();

  void SetAllelicAssociationTest(const std::vector<double> &alpha0);

  void Update(double);
  void UpdateScoresForResidualAllelicAssociation(const array_of_allelefreqs& Allelefreqs);

  ~ScoreTests();

private:
  double* SumAncestryScore;
  double* SumAncestryInfo;
  double* SumAncestryVarScore;
  double* SumAncestryScore2;

  double* SumAffectedsScore2;
  double* SumAffectedsVarScore;
  double* SumAffectedsScore;
  double* SumAffectedsInfo;

  double **LocusLinkageAlleleScore;
  double **LocusLinkageAlleleInfo;
  double **SumLocusLinkageAlleleScore2;
  double **SumLocusLinkageAlleleScore;
  double **SumLocusLinkageAlleleInfo;
  int *locusObsIndicator;
  unsigned *dim_;

  double ***ScoreWithinHaplotype;
  double ***InfoWithinHaplotype;
  double **SumScoreWithinHaplotype;
  double **SumScore2WithinHaplotype;
  double **SumInfoWithinHaplotype;

  double* AdmixtureScore; 
  double* AdmixtureInfo; 
  double* SumAdmixtureScore; 
  double* SumAdmixtureScore2; 
  double* SumAdmixtureInfo;

#ifdef PARALLEL
  double *sendallelescore;
  double *sendalleleinfo;
  double *recvallelescore;
  double *recvalleleinfo;
  double *sendhapscore;
  double *sendhapinfo;
  double *recvhapscore;
  double *recvhapinfo;
  int dimallelescore, dimalleleinfo, dimhapscore, dimhapinfo;

  const std::vector<std::string> * LocusLabels;
  const MPI::Intracomm * Comm;//pointer to the workers_and_master communicator in admixmap.cc
#endif

  std::ofstream assocscorestream;
  std::ofstream ancestryAssociationScoreStream;
  std::ofstream HaplotypeAssocScoreStream;
  std::ofstream affectedsOnlyScoreStream;
  std::ofstream allelicAssocScoreStream;
  //const Vector_s& PopLabels;

  const AdmixOptions *options;
  const IndividualCollection *individuals;
  const Genome* Lociptr;//Pointer to Loci
  const Chromosome* const* chrm;//Copy of pointer to array of chromosomes
  int rank, worker_rank, NumWorkers;

//OUTPUT
  void OpenFile(LogWriter &Log, std::ofstream* outputstream, const char* filename, std::string testname);
  void OutputTestsForLocusLinkage( int iteration, ofstream* outputstream, const Vector_s& PopLabels,
			      const double* Score, const double* VarScore,
				   const double* Score2, const double* Info, string sep );

  void OutputScoreTest( int iterations, ofstream* outputstream, unsigned dim, std::vector<std::string> labels,
			const double* score, const double* scoresq, const double* info, bool final, unsigned );

  void OutputScalarScoreTest( int iterations, ofstream* outputstream, string label,
			      const double score, const double scoresq, const double info, bool final);
  void OutputAdmixtureScoreTest( int );


  //void UpdateScoreForWithinHaplotypeAssociation( const Individual* const ind, int locus, double p,double phi, double DInvLink);
  void UpdateScoreForWithinHaplotypeAssociation( const Individual* const ind, const std::vector<int> allele2Counts, 
						 int locus, double p,double phi, double DInvLink);

  void CentreAndSum(unsigned dim, double *score, double* info, 
				double *sumscore, double* sumscoresq, double* suminfo);

  void UpdateAlleleScores( double* score, double* info, const double* admixtureProps, const vector<int> Counts, 
			   double YMinusEY, double phi, double DInvLink);

  void UpdateScoreForAllelicAssociation( const Individual* const , double, double, double, bool);

  void UpdateScoreForAdmixtureAssociation( const double* const Theta, double YMinusEY,double phi, double DInvLink);


  static int ResidualAlleleInfoIndex(int M, int N, int m1, int n1, int m2, int n2);

  static std::string double2R( double );
  static std::string double2R( double x, int precision );

  static void R_output3DarrayDimensions(std::ofstream* stream, const std::vector<int> dim, const std::vector<std::string> labels);    

  void Reset();

  ResidualLDTest ResidualAllelicAssocScoreTest;//here temporarily, until rest of scoretests classes have been created

};






#endif /* !defined SCORETESTS_H */
