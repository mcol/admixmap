// *-*-C++-*-*
#ifndef COMPOSITE_LOCUS_H
#define COMPOSITE_LOCUS_H 1

//#include "AbstractCompLocus.h"
#include "rand.h"
#include "vector.h"
#include "matrix.h"
#include "vector_d.h"
#include "matrix_d.h"
#include "MatrixArray_d.h"
#include "vector_i.h"
#include "matrix_i.h"
#include "MatrixArray_i.h"
#include "VectorLoop.h"
#include "DARS.h"
#include "TuneRW.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

class Vector_i;
class Matrix_d;
class VectorLoop;
class TuneRW;

class CompositeLocus //: public AbstractCompLocus
{
private: // members
  int NumberOfLoci;
  int NumberOfStates;
  int Populations;
  int NumberOfMergedHaplotypes;
  Vector_i NumberOfAlleles;
  MatrixArray_d HaplotypeProbs;
  MatrixArray_d HaplotypeProbsMAP;
  VectorLoop HapLoop;
  VectorLoop DipLoop;
  std::string *Label;
  Vector_i pmult;
  Vector_i MergeHaplotypes;
  Matrix_i HapLabels;
  Vector_i PossHaplotypes;

   //possibly move out
  Matrix_d ScoreGene;
  Matrix_d InfoGene;
  Matrix_d SumScoreGene;
  Matrix_d SumScoreGeneSq;
  Matrix_d SumInfoGene;
  MatrixArray_d SumNewScore;
  MatrixArray_d SumNewScoreSq;
  MatrixArray_d SumNewInfo;

  Vector_i GetHaplotypes(const std::vector<unsigned int>&);
  Matrix_d GetHaploidLocusProbs(const std::vector<unsigned int>& x);
  Vector_i Haplotype( Vector_i, int );
  void SetNoMergeHaplotypes();
  void setHaplotypeProbsMAP(Matrix_d);
  double GetAlleleProbs( int x, int ancestry , Matrix_d &Freqs);

  // UNIMPLEMENTED
  // to avoid use
  CompositeLocus(const CompositeLocus&);
  CompositeLocus& operator=(const CompositeLocus&);

public:
  CompositeLocus();
   ~CompositeLocus();

  Vector_i decodeGenotype(const std::vector<unsigned int>& encoded);
  void ConstructHaplotypeProbs(Matrix_d &AlleleFreqs);
  int HapLoopGetDecimal(Vector_i x);
  Matrix_d GetGenotypeProbs(const std::vector<unsigned int>&,bool, int);

  void GetPossibleHaplotypes(std::vector<unsigned int>& genotype);
 
  void setHaplotypeProbsMAP();
  void SetNumberOfLabels();
  void SetNumberOfPopulations( int );
  void SetNumberOfStates( int );
  void InitialiseHaplotypes(Matrix_d &);
  void InitialiseMuProposal(int);

  Vector_i SampleHaplotype(const std::vector<unsigned int>& genotype, Vector_i , Matrix_d &);
  int GetNumberOfLoci();
  int GetNumberOfStates();
  int GetSize();
  std::string GetLabel(int);
  int GetMergedHaplotype( int i );
  int GetNumberOfMergedHaplotypes();
  Vector_i GetHapLabels( int );
  Vector_i GetNumberOfAlleles();
  int GetNumberOfAllelesOfLocus( int );
  Vector_i GetAlleleCountsInHaplotype(const std::vector<unsigned int>&);
  void SetDefaultMergeHaplotypes( Vector_d alpha, Matrix_d AlleleFreqs );
  void SetLabel( int, std::string );
  void SetNumberOfLoci( int );
  void SetNumberOfAllelesOfLocus( int, int );
  void AddLocus( int );
 

  //score test for misspecified allelefreqs
  void InitialiseScoreTest(int);
  void SumScoreForMisSpecOfAlleleFreqs();
  void ResetScoreForMisSpecOfAlleleFreqs();
  void UpdateScoreForMisSpecOfAlleleFreqs( Matrix_d phi, std::vector<unsigned int> x , Matrix_d );
  void UpdateScoreForMisSpecOfAlleleFreqs2(Matrix_d&, Matrix_i&);
  Matrix_d GetScore();
  Matrix_d GetInfo();
  Matrix_d GetScoreSq();
  Matrix_d GetNewScore( int );
  Matrix_d GetNewInfo( int );
  Matrix_d GetNewScoreSq( int );



};

double GetMarginalLikelihood( Vector_d PriorAlleleFreqs, Vector_d AlleleCounts );




#endif /* !COMPOSITE_LOCUS_H */
