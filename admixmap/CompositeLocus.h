// *-*-C++-*-*
#ifndef COMPOSITE_LOCUS_H
#define COMPOSITE_LOCUS_H 1

#include "AbstractCompLocus.h"
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

#include "LocusVisitor.h"

class Vector_i;
class Matrix_d;
class VectorLoop;
class TuneRW;

class CompositeLocus : public AbstractCompLocus
{
private: // members
  int NumberOfLoci;
  int NumberOfStates;
  int Populations;
  int RandomAlleleFreqs;
  int NumberOfMergedHaplotypes;
  Vector_i NumberOfAlleles;
  MatrixArray_d HaplotypeProbs;
  MatrixArray_d HaplotypeProbsMAP;
  bool Historical;
  
  VectorLoop HapLoop;
  VectorLoop DipLoop;
  std::string *Label;
  Vector_i pmult;
  Vector_i MergeHaplotypes;
  Vector_d Fst;
  Vector_d SumFst;
  Matrix_i HapLabels;

  Vector_i PossHaplotypes;

  //to move out
  Matrix_d AlleleFreqs;
  Matrix_d AlleleFreqsMAP;
  Matrix_d HistoricalAlleleFreqs;
  Matrix_i LikelihoodAlleleFreqs;
  Matrix_d HistoricLikelihoodAlleleFreqs;
  Matrix_d PriorAlleleFreqs;
  Matrix_d SumAlleleFreqs;
  Vector_d SumEta;

  //possibly move out
  Matrix_d ScoreGene;
  Matrix_d InfoGene;
  Matrix_d SumScoreGene;
  Matrix_d SumScoreGeneSq;
  Matrix_d SumInfoGene;
  MatrixArray_d SumNewScore;
  MatrixArray_d SumNewScoreSq;
  MatrixArray_d SumNewInfo;




  // HELPER METHODS
  Matrix_d GetGenotypeProbs(const std::vector<unsigned int>&,bool);
  Vector_i GetHaplotypes(const std::vector<unsigned int>&);
  void ConstructHaplotypeProbs();
  Matrix_d GetLocusProbs(const std::vector<unsigned int>& x,bool);
  Matrix_d GetHaploidLocusProbs(const std::vector<unsigned int>& x);
  double GetAlleleProbsMAP( int, int );
  Vector_i Haplotype( Vector_i, int );
  void SetNoMergeHaplotypes();
  void SamplePriorAlleleFreqs1D( Vector_d eta );
  void SamplePriorAlleleFreqsMultiDim( Vector_d eta );
//    DARS SampleMu;
   std::vector<TuneRW> MuProposal;
  
  // UNIMPLEMENTED
  // to avoid use
  CompositeLocus(const CompositeLocus&);
  CompositeLocus& operator=(const CompositeLocus&);

public:
  CompositeLocus();
   virtual ~CompositeLocus();

 static Vector_i
 decodeGenotype(const std::vector<unsigned int>& encoded);

  void GetPossibleHaplotypes(std::vector<unsigned int>& genotype);

  virtual void
  accept(LocusVisitor&);

  virtual Vector_d GetPriorAlleleFreqs( int );
  virtual int GetNumberOfLoci();
  virtual int GetNumberOfStates();
  virtual std::string GetLabel(int);
  virtual int IsRandom();
  virtual void SetDefaultMergeHaplotypes( Vector_d alpha );
  virtual int GetMergedHaplotype( int i );
  virtual int GetNumberOfMergedHaplotypes();
   virtual Vector_i GetHapLabels( int );
  virtual void ResetSumAlleleFreqs();
  virtual Vector_d GetFst();
  virtual Matrix_d GetLikelihood( std::vector<unsigned int>, bool, bool );
  virtual void UpdateLikelihoodAlleleFreqs(const std::vector<unsigned int>&, Vector_i );
  virtual void
  UpdateLikelihoodAlleleFreqs_HaploidData(const std::vector<unsigned int>&, int );
  virtual Vector_i SampleHaplotype(const std::vector<unsigned int>& genotype, Vector_i );
  virtual Vector_i GetAlleleCountsInHaplotype(const std::vector<unsigned int>&);
  virtual void SumScoreForMisSpecOfAlleleFreqs();
  virtual Matrix_d GetAlleleFreqs();
   Vector_d getAlleleFreqsMAP(int);
   void setHaplotypeProbsMAP();
   void setAlleleFreqsMAP();
  virtual void UpdateScoreForMisSpecOfAlleleFreqs( Matrix_d phi, std::vector<unsigned int> x );
  void UpdateScoreForMisSpecOfAlleleFreqs2();

  virtual int GetSize();

  virtual Matrix_d GetScore();
  virtual Matrix_d GetInfo();
  virtual Matrix_d GetScoreSq();
  virtual Matrix_d GetSumAlleleFreqs();
  virtual Matrix_i GetLikelihoodAlleleFreqs();
   Vector_i GetLikelihoodAlleleFreqs(int);
  virtual void SamplePriorAlleleFreqs( Vector_d eta );
  virtual void SampleAlleleFreqs( int );
  virtual Vector_d GetStatsForEta( int );
  virtual void UpdatePriorAlleleFreqs( int, const Vector_d& );
  virtual void UpdatePriorAlleleFreqsGlobal( int, const std::vector<Vector_d>& ){};
  virtual void ResetLikelihoodAlleleFreqs();
  virtual void UpdateFst();
  virtual void SetLabel( int, std::string );
  virtual void SetNumberOfLoci( int );
  virtual void SetNumberOfAllelesOfLocus( int, int );
  virtual void AddLocus( int );
  virtual void SetAlleleFreqs( Matrix_d );
  virtual void SetHistoricalAlleleFreqs( Matrix_d );
  virtual void SetPriorAlleleFreqs( Matrix_d, bool );
  virtual void SetDefaultAlleleFreqs( int );
  virtual Vector_i GetNumberOfAlleles();
  virtual int GetNumberOfAllelesOfLocus( int );
  virtual void ResetScoreForMisSpecOfAlleleFreqs();
  Matrix_d GetNewScore( int );
  Matrix_d GetNewInfo( int );
  Matrix_d GetNewScoreSq( int );
   void SetNumberOfLabels();
   void SetLikelihoodAlleleFreqs( int, Matrix_i );
   void SetNumberOfPopulations( int );
   void SetNumberOfStates( int );
  double GetAlleleProbs( int, int );
};

double GetMarginalLikelihood( Vector_d PriorAlleleFreqs, Vector_d LikelihoodAlleleFreqs );
double fMu( Vector_d &, MatrixArray_i &, MatrixArray_d &, double );
double dfMu( Vector_d &, MatrixArray_i &, MatrixArray_d &, double );
double ddfMu( Vector_d &, MatrixArray_i &, MatrixArray_d &, double );

#endif /* !COMPOSITE_LOCUS_H */
