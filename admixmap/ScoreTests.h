// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   admixmap.cc 
 *   header file for Scoretests class
 *   Copyright (c) 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
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
		  const Chromosome* const*, const std::string *, LogWriter &);

  void InitialiseAssocScoreFile(const std::string *);

  void Output(int, const std::string *);

  void ROutput();

  void SetAllelicAssociationTest(const std::vector<double> &alpha0);

  void Update(double);

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

  std::ofstream assocscorestream;
  std::ofstream * ancestryAssociationScoreStream;
  std::ofstream * SNPsAssociationScoreStream;
  std::ofstream * affectedsOnlyScoreStream;
  std::ofstream genescorestream;
  const std::string *PopLabels;

  const AdmixOptions *options;
  const IndividualCollection *individuals;
  const Genome* Lociptr;//Pointer to Loci member of Latent
  const Chromosome* const* chrm;//Copy of pointer to array of chromosomes

//OUTPUT
  
  void OutputTestsForLocusLinkage( int iteration, ofstream* outputstream,
			      const double* Score, const double* VarScore,
			      const double* Score2, const double* Info );

  void OutputTestsForAllelicAssociation( int iteration, int locus, unsigned dim, const double* score, const double* scoresq, 
					 const double* info);
  
  void OutputTestsForSNPsInHaplotype( int );
  
  void OutputAdmixtureScoreTest( int );

  //void UpdateScoreForWithinHaplotypeAssociation( const Individual* const ind, int locus, double p,double phi, double DInvLink);
  void UpdateScoreForWithinHaplotypeAssociation( const Individual* const ind, const std::vector<int> allele2Counts, 
						 int locus, double p,double phi, double DInvLink);

  void ScoreTests::CentreAndSum(unsigned dim, double *score, double* info, 
				double *sumscore, double* sumscoresq, double* suminfo);

  void UpdateAlleleScores( double* score, double* info, const double* admixtureProps, const vector<int> Counts, 
			   double YMinusEY, double phi, double DInvLink);

  void UpdateScoreForAllelicAssociation( const Individual* const , double, double, double);

  void UpdateScoreForAdmixtureAssociation( const double* const Theta, double YMinusEY,double phi, double DInvLink);

  static std::string double2R( double );

  static void R_output3DarrayDimensions(std::ofstream* stream, const std::vector<int> dim, const std::vector<std::string> labels);    

  void Reset();


};






#endif /* !defined SCORETESTS_H */
