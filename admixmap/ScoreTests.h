// *-*-C++-*-*
#ifndef SCORETESTS_H
#define SCORETESTS_H 1

#include <sstream>
#include "IndividualCollection.h"
#include "LogWriter.h"

class ScoreTests{

private:
// these objects are composite-locus specific and could be moved there
  MatrixArray_d LocusLinkageScore;
  MatrixArray_d LocusLinkageInfo;
  MatrixArray_d SumLocusLinkageScore;
  MatrixArray_d SumLocusLinkageScore2;
  MatrixArray_d SumLocusLinkageInfo;
  MatrixArray_d AffectedsScore;
  MatrixArray_d AffectedsVarScore;
  MatrixArray_d AffectedsInfo;
  MatrixArray_d LocusLinkageAlleleScore;
  MatrixArray_d LocusLinkageAlleleInfo;
  MatrixArray_d SumAffectedsScore2;
  MatrixArray_d SumAffectedsVarScore;
  MatrixArray_d SumAffectedsScore;
  MatrixArray_d SumAffectedsInfo;
  MatrixArray_d SumLocusLinkageAlleleScore2;
  MatrixArray_d SumLocusLinkageAlleleScore;
  MatrixArray_d SumLocusLinkageAlleleInfo;
  MatrixArray_d *ScoreWithinHaplotype;
  MatrixArray_d *InfoWithinHaplotype;
  MatrixArray_d SumScoreWithinHaplotype;
  MatrixArray_d SumScore2WithinHaplotype;
  MatrixArray_d SumInfoWithinHaplotype;

  Vector_d pp;
  Matrix_d U; 
  Matrix_d I; 
  Matrix_d SumU; 
  Matrix_d SumU2; 
  Matrix_d SumI;


  std::ofstream allelefreqscorestream;
  std::ofstream allelefreqscorestream2;
  std::ofstream assocscorestream;
  std::ofstream * ancestryAssociationScoreStream;
  std::ofstream * SNPsAssociationScoreStream;
  std::ofstream * affectedsOnlyScoreStream;
  std::ofstream genescorestream;
  std::string *PopLabels;

  AdmixOptions *options;
  IndividualCollection *individuals;
  Genome* Lociptr;//Pointer to Loci member of Latent
  Genome* chrm; //population

  LogWriter *Logptr;
//OUTPUT
  void OutputTestsForLocusLinkage( int, std::ofstream *, MatrixArray_d, MatrixArray_d, MatrixArray_d );
  
  void OutputTestsForLocusLinkage2( int, std::ofstream *,MatrixArray_d, MatrixArray_d, MatrixArray_d, MatrixArray_d );

  void OutputTestsForAllelicAssociation( int );
  
  void OutputTestsForSNPsInHaplotype( int );
  
  void OutputTestsForMisSpecifiedAlleleFreqs( int);
  
  void OutputTestsForMisSpecifiedAlleleFreqs2( int);
  
  void OutputScoreTest( int );

  void UpdateScoreForWithinHaplotypeAssociation( int ind, int locus, int indicator, double p,double );
  
  void SumScoreForWithinHaplotypeAssociation();

  void UpdateScoreForLinkageAffectedsOnly( Individual*);
  
  void UpdateScoreForLinkage( Individual*, int, int, int, MatrixArray_d *, MatrixArray_d *,double,double);
  
  void UpdateScoreForAssociation_Binary( int, Matrix_d, Matrix_d *, Matrix_d *, double );
  
  void UpdateScoreForAssociation_Continuous( int, Matrix_d, Matrix_d *, Matrix_d *, double ,double);

  void TransformScoreStatistics( int, Matrix_d, Matrix_d, Matrix_d *, Matrix_d * );

public:
  ScoreTests();

  void Initialise(AdmixOptions *,IndividualCollection *,Genome *,Genome *,LogWriter *);

  void InitialiseAssocScoreFile(std::string *);

  void Output(int,std::string *);

  void ROutput();

  static std::string  double2R( double );

  static void 
  R_output3DarrayDimensions(std::ofstream* stream,std::vector<int> dim,std::vector<std::string> labels);    
//UPDATES
  void Reset();

  void Output2(Vector_d *, std::ofstream *);//needs a more descriptive name

  void Update(int, Matrix_d*, double);
  //void UpdateScoreIndLevel(int , Matrix_d*,double);

  //void UpdateScoreStats();

  void UpdateScoresForMisSpecOfAlleleFreqs( int i );

  ~ScoreTests();

};//end of class declaration






#endif /* !defined SCORETESTS_H */
