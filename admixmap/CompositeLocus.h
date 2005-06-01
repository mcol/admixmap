// *-*-C++-*-*
#ifndef COMPOSITE_LOCUS_H
#define COMPOSITE_LOCUS_H 1

#include "rand.h"
#include "common.h"
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>


#include "vector.h"
#include "matrix.h"
#include "vector_d.h"
#include "matrix_d.h"
#include "vector_i.h"
#include "matrix_i.h"
#include "TuneRW.h"

class Vector_i;
class Matrix_d;
class TuneRW;

typedef struct  
{
   int haps[2];
} hapPair; 

class CompositeLocus 
{


public:
  CompositeLocus();
  ~CompositeLocus();

  void setHaplotypeProbsMAP();
  void SetNumberOfLabels();
  void SetNumberOfPopulations( int );
  void SetNumberOfStates( int );
  void SetLabel( int, std::string );
  void SetNumberOfLoci( int );
  void SetNumberOfAllelesOfLocus( int, int );
  void AddLocus( int );
  int GetNumberOfLoci();
  int GetNumberOfStates();
  std::string GetLabel(int);

  int GetNumberOfAllelesOfLocus( int );

  void getLocusAlleleProbs(double **P, int k);

  void setPossibleHaplotypePairs(unsigned short **Genotype, std::vector<hapPair> &PossibleHapPairs);
  void decodeIntAsHapAlleles(const int h, int *hapAlleles);
  void GetGenotypeProbs(double **Probs, std::vector<hapPair > &HaplotypePairs, bool fixed, int RandomAlleleFreqs);
  void SetHapPairProbs();
  void SampleHapPair(int hap[2], std::vector<hapPair > &HapPairs, Vector_i ancestry);
  void SampleHapPair(int hap[2], std::vector<hapPair > &HapPairs, int *ancestry);
  void Initialise(Matrix_d &);
  void SetAlleleProbs(Matrix_d &alleleFreqs);

  //to be rewritten
  int HapLoopGetDecimal(Vector_i x);

  //functions used for haplotype association score test 
  int GetMergedHaplotype( int i );
  int GetNumberOfMergedHaplotypes();
  Vector_i GetHapLabels( int );
  void SetDefaultMergeHaplotypes( Vector_d alpha, Matrix_d AlleleFreqs );
  void SetDefaultMergeHaplotypes( double *alpha, Matrix_d AlleleFreqs );

private: 
  int NumberOfLoci;
  int NumberOfStates;
  int Populations;
  int *NumberOfAlleles;
  dmatrix AlleleProbs;
  double ****HapPairProbs; //haplotype pair probabilities
  double ****HapPairProbsMAP; //Posterior estimates of hap pair probs
  std::string *Label;
  int *base;

  //possibly move out
  Vector_i MergeHaplotypes;
  Matrix_i HapLabels;
  int NumberOfMergedHaplotypes;

  Matrix_d GetHaploidLocusProbs(const std::vector<unsigned int>& x);
  void SetNoMergeHaplotypes();
  void setHaplotypeProbsMAP(Matrix_d);
  double GetAlleleProbs( int x, int ancestry , Matrix_d &Freqs);

  void intToBits(int n, const int length, bool *bits) ;
  void setBaseForHapCode();
  void setBaseMissing(const int *missingLoci, const int numMissingLoci, int baseMissing[][2]);
  void setMissingAlleles(const int baseMissing[][2], const int numMissingLoci, const int permMissing, int MissingAlleles[][2]) ;
  int codeHapAllelesAsInt(const int *hapAlleles);
  void codeHapAllelesPairAsIntPair(const int HapAllelesPair[][2], int *hpair);
  void permuteHetLoci(const bool *isHet, const int numHetLoci, const int permHet, 
		      unsigned short **Genotype, int HapAllelesPair[][2]);
  void permuteMissingLoci(const bool *isMissing, const int numMissingLoci, const int permMissing, 
			  const int HapAllelesPair[][2], const int baseMissing[][2], int HapAllelesPairNoMissing[][2]) ;
  // UNIMPLEMENTED
  // to avoid use
  CompositeLocus(const CompositeLocus&);
  CompositeLocus& operator=(const CompositeLocus&);
};

double GetMarginalLikelihood( Vector_d PriorAlleleFreqs, Vector_d AlleleCounts );




#endif /* !COMPOSITE_LOCUS_H */
