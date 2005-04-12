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
  double **Likelihood;
  Matrix_d Prob;//used to construct likelihood
  int *CodedStates;//used to sample hidden states from HMM
  // Matrix_d *AncestryProbs;//Conditional Probs of Locus Ancestry
  
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
  int GetSize();
  Vector_i SampleForHaploidLocusAncestry(Individual*, AlleleFreqs *);
  void UpdateParameters(Individual*,AlleleFreqs *, Matrix_d&,AdmixOptions*,std::vector<Vector_d>&, bool);
  void UpdateParametersHaploid(Individual*,AlleleFreqs *,Matrix_d&,AdmixOptions*,std::vector<Vector_d>&, bool);
  Matrix_i SampleForLocusAncestry();
  void setAncestryProbs(int); 
  Matrix_d getAncestryProbs(int);
  void CreateAncestryProbs();
  double getLogLikelihood();

};

#endif /* !defined CHROMOSOME_H */
