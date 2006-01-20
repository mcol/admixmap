// *-*-C++-*-*
//To sample allele freqs for a single compound locus
#include "common.h"
#include <math.h>
#include <stdlib.h>
#include "functions.h"
#include "HamiltonianMonteCarlo.h"
#include "IndividualCollection.h"
#include <algorithm>
#include <gsl/gsl_linalg.h>

// typedef struct //from CompositeLocus.h 
// {
//   int haps[2];
// } hapPair; 


typedef struct{
  unsigned NumPops;// #populations
  unsigned NumStates;// #alleles/haplotypes
  unsigned locus;// current locus
  const IndividualCollection* IP;//pointer to individuals
  const double* PriorParams;//parameters of Dirichlet prior on allele freqs

}AlleleFreqArgs;

class AlleleFreqSampler{
public:
  AlleleFreqSampler();
  void SampleAlleleFreqs(double *phi, const double* Prior, IndividualCollection* IC, unsigned locus, 
					  unsigned NumStates, unsigned NumPops);



private:
  HamiltonianMonteCarlo Sampler;
  AlleleFreqArgs Args;

  static double logLikelihood(const double *phi, const int Anc[2], const std::vector<hapPair > H, unsigned NumPops);
  static double logPrior(const double* PriorParams, const double* phi, unsigned NumPops, unsigned NumStates);
  static double logJacobian(const double* a, const double z, unsigned H);
  static void logLikelihoodFirstDeriv(double *phi, const int Anc[2], const std::vector<hapPair > H, unsigned NumStates, unsigned NumPops,
			       double* FirstDeriv);

  static double getEnergy(const double * const phi, const void* const vargs);
  static void gradient(const double * const phi, const void* const vargs, double* g);
};
