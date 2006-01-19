//compile with:
//g++ AlleleFreqSampler.cc HamiltonianMonteCarlo.cc StepSizeTuner.cc functions.cc rand_serial.cc -lgsl -lgslcblas -o AlleleFreqSampler


//To sample allele freqs for a single compound locus
#include "common.h"
#include <math.h>
#include <stdlib.h>
#include "functions.h"
#include "HamiltonianMonteCarlo.h"
#include <algorithm>
#include <gsl/gsl_linalg.h>

typedef struct //from CompositeLocus.h 
{
  int haps[2];
} hapPair; 
HamiltonianMonteCarlo Sampler;

typedef struct{
  unsigned NumPops;
  unsigned NumStates;
  const int AncestryPair[2];
  const std::vector<hapPair> HapPairs;
  const double* PriorParams;

}AlleleFreqArgs;

//to begin with, just one individual, then sum over individuals

// requires: sampled ancestry pair A, PossibleHapPairs (compatible with genotype) H, current values of AlleleFreqs at this locus, phi
//, number of alleles/haplotypes NumStates, number of populations, NumPops

double ::logLikelihood(const double *phi, const int Anc[2], const std::vector<hapPair > H, unsigned NumStates, unsigned NumPops){
  double sum = 0.0, sum2 = 0.0;//sums of products and products of squares
  double phiphi;
  unsigned NumPossHapPairs = H.size();
  for(unsigned hpair = 0; hpair < NumPossHapPairs; ++hpair){
    unsigned j0 = H[hpair].haps[0];//j
    unsigned j1 = H[hpair].haps[1];//j'
    int index0 = j0 * NumPops +Anc[0];//jk
    int index1 = j1 * NumPops +Anc[1];//j'k'
    phiphi = phi[index0] * phi[index1];
    sum += phiphi;
    sum2 += phiphi*phiphi;
  }
  return log(sum2) - log(sum);
}

double ::logPrior(const double* PriorParams, const double* phi, unsigned NumPops, unsigned NumStates){
  double logprior = 0.0;
  std::vector<double> DirichletParams(NumStates);
  std::vector<double> freqs(NumStates);

  for(unsigned k = 0; k < NumPops; ++k){
    for(int s = 0; s < NumStates; ++s){
    DirichletParams[s] = PriorParams[k*NumStates + s];
    freqs[s] = phi[k*NumStates + s];
    }

    logprior += getDirichletLogDensity( DirichletParams, freqs );
  }
  return logprior;
}

double ::logJacobian(const double* a, const double z, unsigned H){
  //computes logJacobian for softmax transformation

  //construct matrix
  gsl_matrix *J = gsl_matrix_calloc(H-1, H-1);
  for(unsigned i = 0; i < H-1; ++i){
    gsl_matrix_set(J,i,i, a[i]*(z-a[i])/(z*z));//diagonal elements
    for(unsigned j = i+1; j < H-1; ++j){
      gsl_matrix_set(J, i, j, -exp(2*a[i]+a[j])/(z*z));//upper triangle
      gsl_matrix_set(J, j, i, -exp(2*a[j]+a[i])/(z*z));//lower triangle
    }
  }
  //LU decomposition
  gsl_permutation *p = gsl_permutation_alloc(H-1);
  gsl_permutation_init(p);
  int signum =1;
  
  //int status = 
  gsl_linalg_LU_decomp ( J , p, &signum);

  gsl_permutation_free(p);
  double logJ = gsl_linalg_LU_lndet(J); 
  gsl_matrix_free(J);
  return logJ; 
}

//energy function for Hamiltonian sampler
double ::getEnergy(const double * const phi, const void* const vargs){
  const AlleleFreqArgs* args = (const AlleleFreqArgs*)vargs;
  double energy = 0.0;
  energy -= logPrior(args->PriorParams, phi, args->NumPops, args->NumStates) 
    + logLikelihood(phi, args->AncestryPair, args->HapPairs, args->NumStates, args->NumPops);
  for(unsigned k = 0; k < args->NumPops; ++k){
    double freqs[args->NumStates];
    double z = 0.0;
    for(int s = 0; s < args->NumStates; ++s){
      freqs[s] = phi[k*args->NumStates + s];
      z += exp(freqs[s]);
    }   
    energy -= logJacobian(freqs, z, args->NumStates);

  }
  return energy;
}

//first derivative of log likelihood
void ::logLikelihoodFirstDeriv(double *phi, const int Anc[2], const std::vector<hapPair > H, unsigned NumStates, unsigned NumPops,
			       double* FirstDeriv){
  unsigned NumPossHapPairs = H.size();
  unsigned dim = NumStates*NumPops;

  double A[dim], B[dim], C[dim], D[dim], E[dim], F[dim];
  //fill(A, A+dim, 0.0);  fill(B, B+dim, 0.0); fill(C, C+dim, 0.0);
  //fill(D, D+dim, 0.0);  fill(E, E+dim, 0.0); fill(F, F+dim, 0.0);
  for(unsigned d = 0; d < dim; ++d){
    A[d] = B[d] = C[d] = D[d] = E[d] = F[d] = 0.0;
  }

  double sum = 0.0, sum2 = 0.0;//sums of products and products of squares
  double phiphi;
  for(unsigned hpair = 0; hpair < NumPossHapPairs; ++hpair){
    unsigned j0 = H[hpair].haps[0];//j
    unsigned j1 = H[hpair].haps[1];//j'
    int index0 = j0 * NumPops +Anc[0];//jk
    int index1 = j1 * NumPops +Anc[1];//j'k'
    phiphi = phi[index0] * phi[index1];
    sum += phiphi;
    sum2 += phiphi*phiphi;

    if(Anc[0] == Anc[1]){
      A[index0] = A[index1] = 1.0;
      D[index0] = D[index1] = 1.0;
    }
    B[index0] += phi[index1];
    B[index1] += phi[index0];
    E[index0] += phi[index1]*phi[index1];
    E[index1] += phi[index0]*phi[index0];

  }
  for(unsigned d = 0;d < dim; ++d){
    double phi2 = phi[d]*phi[d];//phi^2
    double phi3 = phi[d]*phi2;//phi^3
    double phi4 = phi2*phi2;//phi^4

    C[d] = sum - A[d]*phi2 - B[d]*phi[d];
    F[d] = sum2 - D[d]*phi4 - E[d]*phi2;

    FirstDeriv[d] = ( (4*D[d]*phi3 + 2*E[d]*phi[d]) / sum2 ) - ( (2*A[d]*phi[d] + B[d]) / sum );
  }

}

//adds derivative of log prior to derivative of log likelihood
void ::addLogPriorDeriv(double *phi, const double* PriorParams, unsigned NumStates, unsigned NumPops,
			       double* FirstDeriv){

}

void gradient(const double * const phi, const void* const vargs, double* g){
  const AlleleFreqArgs* args = (const AlleleFreqArgs*)vargs;
  logLikelihoodFirstDeriv(phi, args->AncestryPair, args->HapPairs, args->NumStates, args->NumPops, g);
  addLogPriorDeriv(phi, args->PriorParams, args->NumStates, args->NumPops, g);
  //TODO: g should be negative of this
}

void ::SampleAlleleFreqs(double *phi, const int Anc[2], const std::vector<hapPair > H, unsigned NumStates, unsigned NumPops){
  unsigned dim = NumStates*NumPops;

  //initialise Hamiltonian

  //transform phi 

  //call Sample on transformed variables 

  //reverse transformation
}

int main(){


  return 0;
}
