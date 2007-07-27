// *-*-C++-*-*
/* 
 *   AllelicAssocSubTest.h 
 *   header file for AllelicAssocSubTest class
 *   Copyright (c) 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef ALLELICASSOCSUBTEST_H
#define ALLELICASSOCSUBTEST_H 1

#include <string>
#include <vector>
class CompositeLocus;
namespace bclib{
  class DelimitedFileWriter;
}

///base class for tests in the AllelicAssocTest class
class AllelicAssocSubTest{
public:
  AllelicAssocSubTest(unsigned d);
  virtual ~AllelicAssocSubTest();
  virtual void Reset();
  virtual void Resize(unsigned){};
  virtual void Update(const int* const happair, CompositeLocus* const Locus, const double* covariates,  
			    double YMinusEY, double phi, double DInvLink) = 0;
  virtual void Accumulate() = 0;
  virtual void Output(bclib::DelimitedFileWriter& outfile, const CompositeLocus* const Locus, 
		      bool final, unsigned numUpdates) = 0;
  static void SetNumCovars(unsigned);

  unsigned getDim()const;
protected:
  double* Score;
  double* Info;
  unsigned dim;
  static unsigned NumCovars;

  AllelicAssocSubTest();

  void Update(const std::vector<int> Counts, const double* covariates,  
	      double YMinusEY, double phi, double DInvLink);
  static void CentreAndSum(unsigned dim, double *score, double* info, 
			   double *sumscore, double* sumscoresq, double* suminfo);

  static void OutputScoreTest(bclib::DelimitedFileWriter& outputstream, unsigned dim, std::vector<std::string> labels,
			      const double* score, const double* scoresq, const double* info, 
			      bool final, unsigned, unsigned numUpdates );
  
  static void OutputScalarScoreTest( bclib::DelimitedFileWriter& outputstream, std::string label,
				     const double score, const double scoresq, const double info, 
				     bool final, unsigned numUpdates);

};

///allelic association test for SNP
class SNPTest : public AllelicAssocSubTest{
public:
  SNPTest();
  ~SNPTest();
  void Update(const int* const happair, CompositeLocus* const Locus, const double* covariates,  
	      double YMinusEY, double phi, double DInvLink);
  void Accumulate();
  void Output(bclib::DelimitedFileWriter& outfile, const CompositeLocus* const Locus, 
	      bool final, unsigned numUpdates);
private:
  double SumScore;
  double SumScore2;
  double SumInfo;

};

///allelic association test for multiallelic locus
class MultiAllelicLocusTest : public AllelicAssocSubTest{
public:
  MultiAllelicLocusTest(unsigned d);
  ~MultiAllelicLocusTest();
  virtual void Update(const int* const happair, CompositeLocus* const Locus, const double* covariates,  
		      double YMinusEY, double phi, double DInvLink);
  virtual void Accumulate();
  virtual void Output(bclib::DelimitedFileWriter& outfile, const CompositeLocus* const, 
		      bool final, unsigned numUpdates);
protected:
  double* SumScore;
  double* SumScore2;
  double* SumInfo;

private:
  MultiAllelicLocusTest();
};

///haplotype association test
class HaplotypeTest : public MultiAllelicLocusTest{
public:
  HaplotypeTest(unsigned d);
  ~HaplotypeTest();
  void Resize(unsigned d);
  void Update(const int* const happair, CompositeLocus* const Locus, const double* covariates,  
	      double YMinusEY, double phi, double DInvLink);
  void Accumulate();
  void Output(bclib::DelimitedFileWriter& outfile, const CompositeLocus* const, 
	      bool final, unsigned numUpdates);
private:
  HaplotypeTest();
};

///Test for within-haplotype association
class WithinHaplotypeTest : public AllelicAssocSubTest{
public:
  WithinHaplotypeTest(unsigned d);
  ~WithinHaplotypeTest();
  void Reset();
  void Resize(unsigned){};
  void Update(const int* const happair, CompositeLocus* const Locus, const double* covariates,  
	      double YMinusEY, double phi, double DInvLink);
  void Accumulate();
  void Output(bclib::DelimitedFileWriter& outfile, const CompositeLocus* const Locus, 
	      bool final, unsigned numUpdates);

private:
  double** WithinHaplotypeScore;
  double** WithinHaplotypeInfo;
  double* SumWithinHaplotypeScore;
  double* SumWithinHaplotypeScore2;
  double* SumWithinHaplotypeInfo;

  WithinHaplotypeTest();
  void UpdateWithinHaplotypeAssociationTest( const double* covariates, const std::vector<int> allele2Counts, 
					     double YMinusEY, double phi, double DInvLink);

};


#endif /* !defined ALLELICASSOCSUBTEST_H */
