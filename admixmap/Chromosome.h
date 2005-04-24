// *-*-C++-*-*
#ifndef CHROMOSOME_H
#define CHROMOSOME_H 1

#include "Genome.h"
#include "AdmixOptions.h"
#include "matrix_d.h"
#include "matrix_i.h"
#include "MatrixArray_d.h"
#include "MatrixArray_i.h"
#include "vector_d.h"
#include "HMM.h"
#include "Latent.h"
#include "AlleleFreqs.h"
#include <vector>

class Individual;

class Chromosome:public Genome
{
private:
  int _startLoci;
  int populations;
  int D;
  int L;
  std::string _Label;
  HMM SampleStates;
  double *StationaryDist;
  double **GenotypeProbs;
  Matrix_d Prob;//used to construct genotypeprobs
  int *CodedStates;//used to sample hidden states from HMM
  double **Tpat, **Tmat;//paternal and maternal transition probability matrices
  double *_product1, *_product2;
  
  // UNIMPLEMENTED
  // to avoid use
 // Private default constructor
  Chromosome(const Chromosome&);
  Chromosome& operator=(const Chromosome&);

public:
  Chromosome();
  Chromosome(int size,int start, int);
  void ResetStuffForX();
  ~Chromosome();
  void SetLabel( int, std::string );
  std::string GetLabel( int );
  int GetLocus(int);
  unsigned int GetSize();
  void UpdateParameters(Individual*,AlleleFreqs *, Matrix_d&,AdmixOptions*,double *[], bool,bool);
  void SampleForLocusAncestry(Matrix_i*, bool);
  void setAncestryProbs(int); 
  void getAncestryProbs(int, double[][3]);
  void CreateAncestryProbs();
  double getLogLikelihood();

};

#endif /* !defined CHROMOSOME_H */
