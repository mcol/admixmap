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

#include "vector_d.h"

class ScoreTests{

public:
  ScoreTests();

  void Initialise(AdmixOptions *,IndividualCollection *,Genome *,Chromosome **,std::string *, LogWriter *);

  void InitialiseAssocScoreFile(std::string *);

  void Output(int,std::string *);

  void ROutput();

  void SetAllelicAssociationTest(Vector_d *alpha0);

  void Update(double);

  ~ScoreTests();

private:

//   Matrix_d *LocusLinkageScore;
//   Matrix_d *LocusLinkageInfo;
//   Matrix_d *SumLocusLinkageScore;
//   Matrix_d *SumLocusLinkageScore2;
//   Matrix_d *SumLocusLinkageInfo;

  Matrix_d SumAncestryScore;
  Matrix_d SumAncestryInfo;
  Matrix_d SumAncestryVarScore;
  Matrix_d SumAncestryScore2;

  Matrix_d *LocusLinkageAlleleScore;
  Matrix_d *LocusLinkageAlleleInfo;
  Matrix_d *SumLocusLinkageAlleleScore2;
  Matrix_d *SumLocusLinkageAlleleScore;
  Matrix_d *SumLocusLinkageAlleleInfo;
  bool *locusObsIndicator;

  Matrix_d SumAffectedsScore2;
  Matrix_d SumAffectedsVarScore;
  Matrix_d SumAffectedsScore;
  Matrix_d SumAffectedsInfo;
  
  Matrix_d **ScoreWithinHaplotype;
  Matrix_d **InfoWithinHaplotype;
  Matrix_d *SumScoreWithinHaplotype;
  Matrix_d *SumScore2WithinHaplotype;
  Matrix_d *SumInfoWithinHaplotype;

  Matrix_d AdmixtureScore; 
  Matrix_d AdmixtureInfo; 
  Matrix_d SumAdmixtureScore; 
  Matrix_d SumAdmixtureScore2; 
  Matrix_d SumAdmixtureInfo;

  std::ofstream assocscorestream;
  std::ofstream * ancestryAssociationScoreStream;
  std::ofstream * SNPsAssociationScoreStream;
  std::ofstream * affectedsOnlyScoreStream;
  std::ofstream genescorestream;
  std::string *PopLabels;

  AdmixOptions *options;
  IndividualCollection *individuals;
  Genome* Lociptr;//Pointer to Loci member of Latent
  Chromosome** chrm;//Copy of pointer to array of chromosomes

  LogWriter *Logptr;
//OUTPUT
  //void OutputTestsForLocusLinkage( int, std::ofstream *, Matrix_d *, Matrix_d *, Matrix_d *);
  
  void OutputTestsForLocusLinkage2( int, std::ofstream *,Matrix_d, Matrix_d, Matrix_d, Matrix_d );

  void OutputTestsForAllelicAssociation( int );
  
  void OutputTestsForSNPsInHaplotype( int );
  
  void OutputAdmixtureScoreTest( int );

  void UpdateScoreForWithinHaplotypeAssociation( Individual *ind, int locus, double p,double phi, double DInvLink);
  Vector_i GetAlleleCountsInHaplotype(unsigned short **genotype, int NumberOfLoci);
  void SumScoreForWithinHaplotypeAssociation();

  void UpdateScoreForLinkageAffectedsOnly( Individual*);
  
  void UpdateScoreForAllelicAssociation( Individual*, double,double, double);

  //void UpdateScoreForAncestryOld( Individual* ind, double Y, int regressonindicator, double EY, double lambda0);

  void UpdateScoreForAdmixtureAssociation( double *Theta, double YMinusEY,double phi, double DInvLink);

  static std::string double2R( double );

  static void 
  R_output3DarrayDimensions(std::ofstream* stream,std::vector<int> dim,std::vector<std::string> labels);    

  void Reset();


};






#endif /* !defined SCORETESTS_H */
