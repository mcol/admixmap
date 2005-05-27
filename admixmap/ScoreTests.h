// *-*-C++-*-*
#ifndef SCORETESTS_H
#define SCORETESTS_H 1

#include <sstream>
#include "IndividualCollection.h"
#include "LogWriter.h"

class ScoreTests{

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

  void UpdateScoreForAdmixtureAssociation( Matrix_d Theta, double YMinusEY,double phi, double DInvLink);

  static std::string double2R( double );

  static void 
  R_output3DarrayDimensions(std::ofstream* stream,std::vector<int> dim,std::vector<std::string> labels);    

  void Reset();

public:
  ScoreTests();

  void Initialise(AdmixOptions *,IndividualCollection *,Genome *,Chromosome **,std::string *, LogWriter *);

  void InitialiseAssocScoreFile(std::string *);

  void Output(int,std::string *);

  void ROutput();

  void SetAllelicAssociationTest();

  void Update(double);

  ~ScoreTests();

};






#endif /* !defined SCORETESTS_H */
