// *-*-C++-*-*
#ifndef CHROMOSOME_H
#define CHROMOSOME_H 1

#include "Genome.h"
#include "matrix_d.h"
#include "matrix_i.h"
#include "vector_d.h"
#include "HMM.h"
#include "AlleleFreqs.h"
#include <vector>

class Individual;
class AdmixOptions;

class Chromosome:public Genome
{
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


  void SampleLocusAncestry(Matrix_i *OrderedStates, Matrix_d &Admixture, bool isdiploid);
  void getAncestryProbs(int, double[][3]);
  double getLogLikelihood();
  void SampleJumpIndicators(const Matrix_i &LocusAncestry, const unsigned int gametes, 
			    int *sumxi, double *Sumrho0, Matrix_i *SumLocusAncestry, Matrix_i *SumLocusAncestry_X, bool isX, 
			    unsigned int SumN[], unsigned int SumN_X[], bool RhoIndicator);
private:
  int _startLocus;
  int populations;
  int D; 
  std::string _Label;
  HMM SampleStates;
  double ***Lambda;
 
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
};

#endif /* !defined CHROMOSOME_H */
