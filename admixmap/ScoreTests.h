// *-*-C++-*-*
#ifndef SCORETESTS_H
#define SCORETESTS_H 1

#include <sstream>
#include "IndividualCollection.h"
#include "LogWriter.h"

class ScoreTests{

private:
// these objects are composite-locus specific and could be moved there
//   MatrixArray_d LocusLinkageScore;
//   MatrixArray_d LocusLinkageInfo;
//   MatrixArray_d SumLocusLinkageScore;
//   MatrixArray_d SumLocusLinkageScore2;
//   MatrixArray_d SumLocusLinkageInfo;

  Matrix_d SumAncestryScore;
  Matrix_d SumAncestryInfo;
  Matrix_d SumAncestryVarScore;
  Matrix_d SumAncestryScore2;

  MatrixArray_d LocusLinkageAlleleScore;
  MatrixArray_d LocusLinkageAlleleInfo;
  MatrixArray_d SumLocusLinkageAlleleScore2;
  MatrixArray_d SumLocusLinkageAlleleScore;
  MatrixArray_d SumLocusLinkageAlleleInfo;

  Matrix_d SumAffectedsScore2;
  Matrix_d SumAffectedsVarScore;
  Matrix_d SumAffectedsScore;
  Matrix_d SumAffectedsInfo;
  
  MatrixArray_d *ScoreWithinHaplotype;
  MatrixArray_d *InfoWithinHaplotype;
  MatrixArray_d SumScoreWithinHaplotype;
  MatrixArray_d SumScore2WithinHaplotype;
  MatrixArray_d SumInfoWithinHaplotype;

  Matrix_d U; 
  Matrix_d I; 
  Matrix_d SumU; 
  Matrix_d SumU2; 
  Matrix_d SumI;

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
  void OutputTestsForLocusLinkage( int, std::ofstream *, MatrixArray_d, MatrixArray_d, MatrixArray_d );
  
  void OutputTestsForLocusLinkage2( int, std::ofstream *,Matrix_d, Matrix_d, Matrix_d, Matrix_d );

  void OutputTestsForAllelicAssociation( int );
  
  void OutputTestsForSNPsInHaplotype( int );
  
  void OutputScoreTest( int );

  void UpdateScoreForWithinHaplotypeAssociation( Individual *ind, int locus, double p,double phi, double DInvLink);
  
  void SumScoreForWithinHaplotypeAssociation();

  void UpdateScoreForLinkageAffectedsOnly( Individual*);
  
  void UpdateScoreForAllelicAssociation( Individual*, double,double, double);

  //void UpdateScoreForAncestryOld( Individual* ind, double Y, int regressonindicator, double EY, double lambda0);

  void UpdateScoreForAssociation( Matrix_d Theta, double YMinusEY,double phi, double DInvLink);


public:
  ScoreTests();

  void Initialise(AdmixOptions *,IndividualCollection *,Genome *,Chromosome **,std::string *, LogWriter *);

  void InitialiseAssocScoreFile(std::string *);

  void Output(int,std::string *);

  void ROutput();

  static std::string  double2R( double );

  static void 
  R_output3DarrayDimensions(std::ofstream* stream,std::vector<int> dim,std::vector<std::string> labels);    

  void Reset();

  void SetAllelicAssociationTest();

  void Update(double);

  ~ScoreTests();

};






#endif /* !defined SCORETESTS_H */
