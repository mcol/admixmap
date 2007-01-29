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
#include "utils/LogWriter.h"
#include "common.h"
#include "AdmixtureAssocTest.h"
#include "AffectedsOnlyTest.h" 
#include "AncestryAssocTest.h"
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

  void Initialise(Options* , const IndividualCollection* const, const Genome* const ,
		  const Vector_s&, LogWriter &);

  void Output(int iterations, const Vector_s& PLabels, const Vector_s& LocusLabels, bool final);

  void ROutput();

  void SetAllelicAssociationTest(const std::vector<double> &alpha0);

  void Update(const vector<Regression* >& R);
  void UpdateScoresForResidualAllelicAssociation(const array_of_allelefreqs& Allelefreqs);

  ~ScoreTests();

  AffectedsOnlyTest& getAffectedsOnlyTest();
  AncestryAssocTest& getAncestryAssocTest();
  void OutputLikelihoodRatios(const char* const filename, int iterations, const Vector_s& PopLabels);

private:
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

  std::ofstream HaplotypeAssocScoreStream;
  std::ofstream allelicAssocScoreStream;

  const Options *options;
  const IndividualCollection *individuals;
  const Genome* Lociptr;//Pointer to Loci
  const Chromosome* const* chrm;//Copy of pointer to array of chromosomes
  int rank, worker_rank, NumWorkers;
  unsigned NumOutputs;//counts calls to output function for dimensions of R objects

//OUTPUT
  void OpenFile(LogWriter &Log, std::ofstream* outputstream, const char* filename, std::string testname);

  void OutputScoreTest( int iterations, ofstream* outputstream, unsigned dim, std::vector<std::string> labels,
			const double* score, const double* scoresq, const double* info, bool final, unsigned );

  void OutputScalarScoreTest( int iterations, ofstream* outputstream, string label,
			      const double score, const double scoresq, const double info, bool final);

  //void UpdateScoreForWithinHaplotypeAssociation( const Individual* const ind, int locus, double p,double phi, double DInvLink);
  void UpdateScoreForWithinHaplotypeAssociation( const Individual* const ind, const std::vector<int> allele2Counts, 
						 int locus, double p,double phi, double DInvLink);

  void CentreAndSum(unsigned dim, double *score, double* info, 
				double *sumscore, double* sumscoresq, double* suminfo);

  void UpdateAlleleScores( double* score, double* info, const double* admixtureProps, const vector<int> Counts, 
			   double YMinusEY, double phi, double DInvLink);

  void UpdateScoreForAllelicAssociation( const Individual* const , double, double, double, bool);

  static std::string double2R( double );
  static std::string double2R( double x, int precision );

  static void R_output3DarrayDimensions(std::ofstream* stream, const std::vector<int> dim, const std::vector<std::string> labels);    

  void Reset();

  ResidualLDTest ResidualAllelicAssocScoreTest;//here temporarily, until rest of scoretests classes have been created
  AdmixtureAssocTest AdmixtureAssocScoreTest;
  AffectedsOnlyTest AffectedsOnlyScoreTest;
  AncestryAssocTest AncestryAssocScoreTest;

  AncestryAssocTest NewAllelicAssocTest;//test using conditional distribution of copies of allele2, rather than sampled values
};






#endif /* !defined SCORETESTS_H */
