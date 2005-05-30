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
  unsigned int L;
  std::string _Label;
  HMM SampleStates;
  double ***Lambda;
  double ***StateArrivalProbs;

  // f0 and f1 are arrays of scalars of the form exp(- rho*x), where x is distance between loci
  // With a global rho model, this array is same for all individuals and calculated only once.
  // required to calculate transition matrices 
  double *f[2]; 
  int *CodedStates;//used to sample hidden states from HMM
  
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
  void SetLabel(std::string );
  std::string GetLabel( int );
  int GetLocus(int);
  unsigned int GetSize();

  void InitialiseLociCorr(const double rho);
  void SetLociCorr(const double rho);

  void UpdateParameters(Individual* ind, AlleleFreqs *A, Matrix_d& Admixture, AdmixOptions* options,  
			double *f[2], bool fixedallelefreqs, bool diploid );

  void UpdateParameters(Individual* ind, AlleleFreqs *A, Matrix_d& Admixture, AdmixOptions* options,  
			std::vector< double > _rho,  bool fixedallelefreqs, bool diploid );


  void SampleLocusAncestry(Matrix_i *OrderedStates, Matrix_d &Admixture,  bool isdiploid);
  void getAncestryProbs(int, double[][3]);
  double getLogLikelihood();
  void SampleJumpIndicators(const Matrix_i &LocusAncestry,   
				      const unsigned int gametes, std::vector< std::vector<bool> > *xi, int *sumxi, 
				      double *Sumrho0);
};

#endif /* !defined CHROMOSOME_H */
