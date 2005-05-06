// *-*-C++-*-*
#ifndef CHROMOSOME_H
#define CHROMOSOME_H 1

#include "Genome.h"
#include "AdmixOptions.h"
#include "matrix_d.h"
#include "matrix_i.h"
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
  //double *StationaryDist;
  //double **GenotypeProbs;
  double ***Lambda;
  //Matrix_d Prob;//used to construct genotypeprobs
  int *CodedStates;//used to sample hidden states from HMM
  // double **Tpat, **Tmat;//paternal and maternal transition probability matrices
  //double *_product1, *_product2;
  
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
  //void UpdateParameters(Individual*,AlleleFreqs *, Matrix_d&,AdmixOptions*,double *[], bool,bool);
  void NewUpdateParameters(Individual* ind, AlleleFreqs *A, Matrix_d& Admixture, AdmixOptions* options, double * f[],
		    bool fixedallelefreqs, bool diploid );
  //void SampleForLocusAncestry(Matrix_i*, bool);
  void NewSampleForLocusAncestry(Matrix_i *OrderedStates, Matrix_d &Admixture, double *f[], int Mcol,bool isdiploid);
  void getAncestryProbs(int, double[][3]);
  double getLogLikelihood();

};

#endif /* !defined CHROMOSOME_H */
