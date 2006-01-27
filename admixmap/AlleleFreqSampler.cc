/** 
 *   ADMIXMAP
 *   AlleleFreqSampler.cc 
 *   Class to sample allele/haplotype frequencies.
 *   Copyright (c) 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "AlleleFreqSampler.h"
#include "IndividualCollection.h"

//#define DEBUG 1

AlleleFreqSampler::AlleleFreqSampler(){
  //default do-nothing constructor
}

//requires: AlleleFreqs phi, parameters of Dirichlet prior Prior, pointer to individuals, current locus number, 
//number of alleles/haplotypes NumStates, number of populations, NumPops
void AlleleFreqSampler::SampleAlleleFreqs(double *phi, const double* Prior, IndividualCollection* IC, unsigned locus, 
					  unsigned NumStates, unsigned NumPops, double coolness){
  unsigned dim = NumStates*NumPops;

  Args.IP = IC;
  Args.NumStates = NumStates;
  Args.NumPops = NumPops;
  Args.locus = locus;
  Args.PriorParams = Prior;
  Args.coolness = coolness;

  //initialise Hamiltonian Sampler
  double step0 = 0.1;//initial step size
  double min = -100.0, max = 100.0; //min and max stepsize
  int numleapfrogsteps = 1;
  Sampler.SetDimensions(dim, step0, min, max, numleapfrogsteps, 0.44/*target acceptrate*/, getEnergy, gradient);

  //transform phi 
  double* params = new double[dim];
  double freqs[NumStates];//frequencies for one population
  for(unsigned k = 0; k < NumPops; ++k){
    freqs[NumStates-1] = 1.0;
    for(unsigned s = 0; s < NumStates - 1; ++s) {//NB supplied frequencies are missing last element
      freqs[s] = phi[s*NumPops+k];
      freqs[NumStates-1] -= freqs[s];
    }
    inv_softmax(NumStates, freqs, params+k*NumStates);
  }

  //call Sample on transformed variables 
  Sampler.Sample(params, &Args);

  //reverse transformation
  for(unsigned k = 0; k < NumPops; ++k){
    softmax(NumStates, freqs, params+k*NumStates);
    for(unsigned s = 0; s < NumStates - 1; ++s) phi[s*NumPops+k] = freqs[s];
  }

  delete[] params;
}

void AlleleFreqSampler::SampleSNPFreqs(double *phi, const double* Prior, const int* AlleleCounts, const int* hetCounts, unsigned locus, 
				  unsigned NumPops, double coolness){
  Args.IP = 0;
  Args.NumStates = 2;
  Args.NumPops = NumPops;
  Args.locus = locus;
  Args.PriorParams = Prior;
  Args.coolness = coolness;
  Args.AlleleCounts = AlleleCounts;
  Args.hetCounts = hetCounts;

  //initialise Hamiltonian Sampler
  double step0 = 0.1;//initial step size
  double min = -100.0, max = 100.0; //min and max stepsize
  int numleapfrogsteps = 20;
  Sampler.SetDimensions(2*NumPops, step0, min, max, numleapfrogsteps, 0.44, getEnergySNP, gradientSNP);

  //transform phi 
  double* params = new double[2*NumPops];
  double freqs[2];//frequencies for one population

  for(unsigned k = 0; k < NumPops; ++k){
    freqs[0] = phi[k];
    freqs[1] = 1.0-phi[k];
    inv_softmax(2, freqs, params+k*2);
  }

  //call Sample on transformed variables 
  Sampler.Sample(params, &Args);

  //reverse transformation
  for(unsigned k = 0; k < NumPops; ++k){
    softmax(2, freqs, params+k*2);
    phi[k] = freqs[0];
  }
  delete[] params;
}

// requires: sampled ancestry pair A, PossibleHapPairs (compatible with genotype) H, current values of AlleleFreqs at this locus, phi
//, number of alleles/haplotypes NumStates, number of populations, NumPops

double AlleleFreqSampler::logLikelihood(const double *phi, const int Anc[2], const std::vector<hapPair > H, 
					unsigned NumStates){
  double sum = 0.0, sum2 = 0.0;//sums of products and products of squares
  double phiphi;
  unsigned NumPossHapPairs = H.size();
  for(unsigned hpair = 0; hpair < NumPossHapPairs; ++hpair){
    unsigned j0 = H[hpair].haps[0];//j
    unsigned j1 = H[hpair].haps[1];//j'
    int index0 = Anc[0]*NumStates + j0;//jk
    int index1 = Anc[1]*NumStates + j1;//j'k'
    phiphi = phi[index0] * phi[index1];
    sum += phiphi;
    sum2 += phiphi*phiphi;
  }
  return log(sum2) - log(sum);
}

double AlleleFreqSampler::logPrior(const double* PriorParams, const double* phi, unsigned NumPops, unsigned NumStates){
  double logprior = 0.0;
  std::vector<double> DirichletParams(NumStates);

  for(unsigned k = 0; k < NumPops; ++k){
    for(unsigned s = 0; s < NumStates; ++s){
    DirichletParams[s] = PriorParams[k*NumStates + s];
    }

    logprior += getDirichletLogDensity( DirichletParams, phi+k*NumStates );
  }
  return logprior;
}

double AlleleFreqSampler::logJacobian(const double* a, const double z, unsigned H){
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
double AlleleFreqSampler::getEnergy(const double * const params, const void* const vargs){
  const AlleleFreqArgs* args = (const AlleleFreqArgs*)vargs;
  double energy = 0.0;
  unsigned States = args->NumStates;

  //transform params to freqs
  double phi[args->NumStates * args->NumPops];
  for(unsigned k = 0; k < args->NumPops; ++k)
    softmax(States, phi + k*States, params + k*States); 

  //accumulate likelihood over individuals
  for(int i = 0; i < args->IP->getSize(); ++i){
    const Individual* ind = args->IP->getIndividual(i);
    int Anc[2];
    ind->GetLocusAncestry(args->locus, Anc);
    energy -= logLikelihood(phi, Anc, ind->getPossibleHapPairs(args->locus), States);
  }
  energy *= args->coolness;//NOTE: here coolness will always be <1 as otherwise the other sampler is in use

  energy -= logPrior(args->PriorParams, phi, args->NumPops, States) ;
  //add Jacobian
  for(unsigned k = 0; k < args->NumPops; ++k){
    double z = 0.0;
    for(unsigned s = 0; s < States; ++s){
      z += exp(phi[s]);
    }   
    energy -= logJacobian(phi+k*States, z, States);
  }
  return energy;
}

void AlleleFreqSampler::gradient(const double * const params, const void* const vargs, double* g){
  const AlleleFreqArgs* args = (const AlleleFreqArgs*)vargs;
  unsigned States = args->NumStates;
  fill(g, g+States* args->NumPops, 0.0);
  //transform params to freqs
  double phi[States * args->NumPops];
  for(unsigned k = 0; k < args->NumPops; ++k)
    softmax(States, phi + k*States, params + k*States); 

  double* dEdphi = new double[States * args->NumPops];fill(dEdphi, dEdphi+States*args->NumPops, 0.0);
  for(int i = 0; i < args->IP->getSize(); ++i){
    const Individual* ind = args->IP->getIndividual(i);
    int Anc[2];
    ind->GetLocusAncestry(args->locus, Anc);
    //first compute gradient wrt phi
    logLikelihoodFirstDeriv(phi, Anc, ind->getPossibleHapPairs(args->locus), States, args->NumPops, dEdphi);
  }
  //subtract derivative of log prior
  for(unsigned s = 0; s < States; ++s)
    for(unsigned k = 0; k < args->NumPops; ++k){
      if(args->PriorParams[s*args->NumPops +k] > 0.0)
	dEdphi[k*States+s] -= (args->PriorParams[s*args->NumPops +k] - 1.0) / phi[k*States+s]; 
    }
  
   //now use chain rule to obtain gradient wrt args
  for(unsigned k = 0; k < args->NumPops; ++k){
    double sum = 0.0;
    for(unsigned s = 0; s < States; ++s){
      for(unsigned s1 = 0; s1 < States; ++s1)sum += exp(params[k*States+s1]);
      for(unsigned s1 = 0; s1 < States; ++s1)
	g[k*States+s] -= dEdphi[k*States+s1] * exp(params[k*States+s1])*exp(params[k*States+s]) / (sum*sum); 

    }
  }
  delete[] dEdphi;
}

//first derivative of  -log likelihood
void AlleleFreqSampler::logLikelihoodFirstDeriv(double *phi, const int Anc[2], const std::vector<hapPair > H, 
						unsigned NumStates, unsigned NumPops, double* FirstDeriv){
  unsigned NumPossHapPairs = H.size();
  unsigned dim = NumStates*NumPops;

  double A[dim], B[dim], /*C[dim],*/ D[dim], E[dim]/*, F[dim]*/;
  //fill(A, A+dim, 0.0);  fill(B, B+dim, 0.0); fill(C, C+dim, 0.0);
  //fill(D, D+dim, 0.0);  fill(E, E+dim, 0.0); fill(F, F+dim, 0.0);
  for(unsigned d = 0; d < dim; ++d){
    A[d] = B[d] = /*C[d] =*/ D[d] = E[d] = /*F[d] =*/ 0.0;
  }

  double sum = 0.0, sum2 = 0.0;//sums of products and products of squares
  double phiphi;
  for(unsigned hpair = 0; hpair < NumPossHapPairs; ++hpair){
    unsigned j0 = H[hpair].haps[0];//j
    unsigned j1 = H[hpair].haps[1];//j'
    int index0 = Anc[0]*NumStates + j0;//jk
    int index1 = Anc[1]*NumStates + j1;//j'k'
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
    //double phi4 = phi2*phi2;//phi^4

    //C[d] = sum - A[d]*phi2 - B[d]*phi[d];
    //F[d] = sum2 - D[d]*phi4 - E[d]*phi2;

    FirstDeriv[d] -= ( (4*D[d]*phi3 + 2*E[d]*phi[d]) / sum2 ) - ( (2*A[d]*phi[d] + B[d]) / sum );
  }

}

//energy function for Hamiltonian sampler, case of SNP
double AlleleFreqSampler::getEnergySNP(const double * const params, const void* const vargs){
  const AlleleFreqArgs* args = (const AlleleFreqArgs*)vargs;
  double energy = 0.0;
  unsigned Pops = args->NumPops;
  //transform params to freqs
  double phi[2 * Pops];
  for(unsigned k = 0; k < Pops; ++k)
    softmax(2, phi + k*2, params + k*2);

  //get loglikelihood
  for(unsigned k = 0; k < Pops; ++k){
    energy -= args->AlleleCounts[k] * log(phi[k*2])/*1k*/ + args->AlleleCounts[Pops + k] * log(phi[k*2+1])/*2k*/;
    for(unsigned k1 = k+1; k1< Pops; ++k1)
      energy -= (args->hetCounts[k*Pops+k1] + args->hetCounts[k1*Pops+k]) * 
	( log(phi[k*2]*phi[k*2]*phi[k1*2+1]*phi[k1*2+1] + phi[k*2+1]*phi[k*2+1]*phi[k1*2]*phi[k1*2]) 
	  - log(phi[k*2]*phi[k1*2+1] + phi[k*2+1]*phi[k1*2]) );
  }

  energy *= args->coolness;//NOTE: here coolness will always be <1 as otherwise the other sampler is in use
  
  energy -= logPrior(args->PriorParams, phi, Pops, 2) ;
  //add Jacobian
  for(unsigned k = 0; k < Pops; ++k){
    double z = 0.0;
    for(unsigned s = 0; s < 2; ++s){
      z += exp(phi[s]);
    }   
    energy -= logJacobian(phi+k*2, z, 2);
  }
  return energy;
}

void AlleleFreqSampler::gradientSNP(const double * const params, const void* const vargs, double* g){
  const AlleleFreqArgs* args = (const AlleleFreqArgs*)vargs;
  fill(g, g+ 2* args->NumPops, 0.0);
  unsigned Pops = args->NumPops;
  //transform params to freqs
  double phi[2 * Pops];
  for(unsigned k = 0; k < Pops; ++k)
    softmax(2, phi + k*2, params + k*2);

  double* dEdphi = new double[2 * Pops];fill(dEdphi, dEdphi+ 2*Pops, 0.0);
  //derivative of log likelihood
  for(unsigned k  = 0; k < Pops; ++k){
    dEdphi[2*k] -= args->AlleleCounts[k] / phi[k*2];//allele1
    dEdphi[2*k+1] -= args->AlleleCounts[Pops + k] / phi[k*2+1];//allele2
    for(unsigned k1 = 0; k1 < Pops; ++k1)if(k!=k1){
      double denom1 = (phi[k*2]*phi[k*2]*phi[k1*2+1]*phi[k1*2+1] + phi[k*2+1]*phi[k*2+1]*phi[k1*2]*phi[k1*2]);
      double denom2 = (phi[k*2]*phi[k1*2+1] + phi[k*2+1]*phi[k1*2]);

      dEdphi[2*k] -= (args->hetCounts[k*Pops+k1] + args->hetCounts[k1*Pops+k]) * 
	((2*phi[k1*2+1]*phi[k1*2+1]*phi[k*2]) / denom1  - (phi[k1*2+1]) / denom2);
      dEdphi[2*k+1] -= (args->hetCounts[k*Pops+k1] + args->hetCounts[k1*Pops+k]) * 	
	((2*phi[k1*2]*phi[k1*2]*phi[k*2+1]) / denom1 - (phi[k1*2]) / denom2);
    }
  }

  //subtract derivative of log prior
  for(unsigned s = 0; s < 2; ++s)
    for(unsigned k = 0; k < Pops; ++k){
      if(args->PriorParams[s*Pops +k] > 0.0)
	dEdphi[k*2+s] -= (args->PriorParams[s*Pops +k] - 1.0) / phi[k*2+s]; 
    }
  
   //now use chain rule to obtain gradient wrt args
  for(unsigned k = 0; k < Pops; ++k){
    double sum = 0.0;
    for(unsigned s = 0; s < 2; ++s){
      for(unsigned s1 = 0; s1 < 2; ++s1)sum += exp(params[k*2+s1]);
      for(unsigned s1 = 0; s1 < 2; ++s1)
	g[k*2+s] -= dEdphi[k*2+s1] * exp(params[k*2+s1])*exp(params[k*2+s]) / (sum*sum); 

    }
  }
  delete[] dEdphi;
}


