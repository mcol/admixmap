// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   admixmap.cc 
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

class ScoreTests{

public:
  ScoreTests();

  void Initialise(AdmixOptions* , const IndividualCollection* const, const Genome* const ,
		  const std::string *, LogWriter &);

  void InitialiseAssocScoreFile(const std::string *);

  void Output(int, const std::string *);
  void WriteFinalTables();

  void ROutput();

  void SetAllelicAssociationTest(const std::vector<double> &alpha0);

  void Update(double);
  void UpdateScoresForResidualAllelicAssociation(const double* const * AlleleFreqs);

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
  bool *locusObsIndicator;
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

  double*** SumAlleleScore;
  double*** SumAlleleScore2;
  double*** SumAlleleInfo;

  std::ofstream assocscorestream;
  std::ofstream ancestryAssociationScoreStream;
  std::ofstream HaplotypeAssocScoreStream;
  std::ofstream affectedsOnlyScoreStream;
  std::ofstream allelicAssocScoreStream;
  std::ofstream ResAlleleScoreFile;
  const std::string *PopLabels;

  const AdmixOptions *options;
  const IndividualCollection *individuals;
  const Genome* Lociptr;//Pointer to Loci
  const Chromosome* const* chrm;//Copy of pointer to array of chromosomes

//OUTPUT
  void OpenFile(LogWriter &Log, std::ofstream* outputstream, const char* filename, std::string testname);
  void OutputTestsForLocusLinkage( int iteration, ofstream* outputstream,
			      const double* Score, const double* VarScore,
				   const double* Score2, const double* Info, string sep );

  void OutputTestsForAllelicAssociation( int iterations, ofstream* outputstream, int locus, unsigned dim, 
					 const double* score, const double* scoresq, const double* info, string sep);
  
  void OutputTestsForHaplotypeAssociation( int iterations, ofstream* outputstream, string sep );
  
  void OutputAdmixtureScoreTest( int );

  void OutputTestsForResidualAllelicAssociation(int iterations, ofstream* outputstream, bool final);

  //void UpdateScoreForWithinHaplotypeAssociation( const Individual* const ind, int locus, double p,double phi, double DInvLink);
  void UpdateScoreForWithinHaplotypeAssociation( const Individual* const ind, const std::vector<int> allele2Counts, 
						 int locus, double p,double phi, double DInvLink);

  void ScoreTests::CentreAndSum(unsigned dim, double *score, double* info, 
				double *sumscore, double* sumscoresq, double* suminfo);

  void UpdateAlleleScores( double* score, double* info, const double* admixtureProps, const vector<int> Counts, 
			   double YMinusEY, double phi, double DInvLink);

  void UpdateScoreForAllelicAssociation( const Individual* const , double, double, double);

  void UpdateScoreForAdmixtureAssociation( const double* const Theta, double YMinusEY,double phi, double DInvLink);


  static int ResidualAlleleInfoIndex(int M, int N, int m1, int n1, int m2, int n2);

  void UpdateScoresForResidualAllelicAssociation(int c, int locus, 
						 const double* const AlleleFreqsA, const double* const AlleleFreqsB);
  void UpdateScoresForResidualAllelicAssociation_1D(int c, int locus,  
						    const double* const AlleleFreqsA, const double* const AlleleFreqsB);
  static std::string double2R( double );
  static std::string double2R( double x, int precision );

  static void R_output3DarrayDimensions(std::ofstream* stream, const std::vector<int> dim, const std::vector<std::string> labels);    

  void Reset();


};






#endif /* !defined SCORETESTS_H */
