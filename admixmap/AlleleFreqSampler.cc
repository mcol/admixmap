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
#include "functions.h"

//#define DEBUG 1

bool AlleleFreqSampler::ishapmixmodel;

///default do-nothing constructor
AlleleFreqSampler::AlleleFreqSampler(){
  params = 0;
}

AlleleFreqSampler::AlleleFreqSampler(unsigned NumStates, unsigned NumPops, 
				     const double* const Prior, bool hapmixmodel = false) {
  unsigned dim = NumStates*NumPops;
  //initialise Hamiltonian Sampler
  double step0 = 0.02;//initial step size
  double min = -100.0, max = 100.0; //min and max stepsize
  int numleapfrogsteps = 5; // 10;
  Args.PriorParams = Prior;
  ishapmixmodel = hapmixmodel;

  if(NumStates == 2){//case of SNP
    params = new double[NumPops];
    step0 = 0.03;//initial step size
    numleapfrogsteps = 20;
    Sampler.SetDimensions(NumPops, step0, min, max, numleapfrogsteps, 0.7, getEnergySNP, 
			  gradientSNP);
  } else {
    params = new double[dim];
    Sampler.SetDimensions(dim, step0, min, max, numleapfrogsteps, 0.7/*target acceptrate*/, 
			  getEnergy, gradient);
  }
}

AlleleFreqSampler::~AlleleFreqSampler(){
  delete[] params;
}

///Samples allele frequencies.
///requires: AlleleFreqs phi, parameters of Dirichlet prior Prior, pointer to individuals, current locus number, 
///number of alleles/haplotypes NumStates, number of populations, NumPops
void AlleleFreqSampler::SampleAlleleFreqs(double *phi,  IndividualCollection* IC, unsigned locus, 
					  unsigned NumStates, unsigned NumPops, double coolness){
  Args.IP = IC;
  Args.NumStates = NumStates;
  Args.NumPops = NumPops;
  Args.locus = locus;
  Args.coolness = coolness;

  //transform phi 
  for(unsigned k = 0; k < NumPops; ++k){
    inv_softmax(NumStates, phi+k*NumStates, params+k*NumStates);
  }

  try { // call Sample on transformed variables 
    Sampler.Sample(params, &Args);
  } catch(string s) {
    throw string("Error sampling allele freqs:\n" + s);
  }

  //reverse transformation
  for(unsigned k = 0; k < NumPops; ++k) {
    softmax(NumStates, phi+k*NumStates, params+k*NumStates);
  }
}

void AlleleFreqSampler::SampleSNPFreqs(double *phi, const int* AlleleCounts, 
				       const int* hetCounts, unsigned locus, unsigned NumPops,
				       double coolness) {
  Args.IP = 0;
  Args.NumStates = 2;
  Args.NumPops = NumPops;
  Args.locus = locus;
  Args.coolness = coolness;
  Args.AlleleCounts = AlleleCounts;
  Args.hetCounts = hetCounts;

  //transform phi 
   //set params = logit(phi)
  for(unsigned k = 0; k < NumPops; ++k){
    params[k] = log(phi[k*2] / (1.0 - phi[k*2]));
  }
  try{
    //call Sample on transformed variables 
    Sampler.Sample(params, &Args);
  }
  catch(string s){
    throw string("Error sampling SNP allele freqs:\n" + s);
  }
  //reverse transformation
  for(unsigned k = 0; k < NumPops; ++k){
    phi[k*2] = exp(params[k]) / (1.0 + exp(params[k]));//allele 1
    phi[k*2+1] = 1.0 - phi[k*2];//allele 2
  }
}

// log normalized prior density - required for updating prior params
double AlleleFreqSampler::logPrior(const double* PriorParams, const double* phi, unsigned NumPops, unsigned NumStates){
  double logprior = 0.0;
  std::vector<double> DirichletParams(NumStates);
  if(ishapmixmodel)
    logprior = NumPops * NumStates * (* PriorParams);
  else{
    for(unsigned k = 0; k < NumPops; ++k){
      for(unsigned s = 0; s < NumStates; ++s){
	DirichletParams[s] = PriorParams[k*NumStates + s];
      }
      logprior += getDirichletLogDensity( DirichletParams, phi+k*NumStates );
    }
  }
  return logprior;
}

///computes logJacobian for softmax transformation
// generic function - should be in functions file
double AlleleFreqSampler::logJacobian(const double* a, const double z, unsigned H){
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
  int signum = 1;
  gsl_linalg_LU_decomp( J , p, &signum);
  gsl_permutation_free(p);
  double logJ = gsl_linalg_LU_lndet(J); 
  gsl_matrix_free(J);
  return logJ; 
}

///energy function for Hamiltonian sampler
double AlleleFreqSampler::getEnergy(const double * const params, const void* const vargs){
  const AlleleFreqArgs* args = (const AlleleFreqArgs*)vargs;
  double energy = 0.0;
  unsigned States = args->NumStates;

  //transform params to freqs
  double *phi = new double[args->NumStates * args->NumPops];
  for(unsigned k = 0; k < args->NumPops; ++k)
     softmax(States, phi + k*States, params + k*States); 

  //accumulate likelihood over individuals
  for(int i = 0; i < args->IP->getSize(); ++i){
    const Individual* ind = args->IP->getIndividual(i);
    int Anc[2];
    ind->GetLocusAncestry(args->locus, Anc);
    energy -= logLikelihood(phi, Anc, ind->getPossibleHapPairs(args->locus), States);
  }
  energy *= args->coolness;

  // subtract log prior density in softmax basis
  if(ishapmixmodel)
    energy -= args->NumPops * States * *(args->PriorParams );
  else{
    for(unsigned k = 0; k < args->NumPops; ++k){
      for( unsigned i = 0; i < States; ++i ) {
	if( args->PriorParams[i] > 0.0 ) {
	  energy -= *(args->PriorParams+k*States+i) * mylog( *(phi+k*States+i) );
	}
      }
    }
  }
  delete[] phi;
  return energy;
}

void AlleleFreqSampler::gradient(const double* const params, const void* const vargs, double* g){
  const AlleleFreqArgs* args = (const AlleleFreqArgs*)vargs;
  unsigned States = args->NumStates;
  fill(g, g+States* args->NumPops, 0.0);
  //transform params to freqs
  double* phi = new double[States * args->NumPops];
  for(unsigned k = 0; k < args->NumPops; ++k) {
    softmax(States, phi + k*States, params + k*States); 
  }
  double* dE_dphi = new double[States * args->NumPops]; 
  fill(dE_dphi, dE_dphi+States*args->NumPops, 0.0);
  for(int i = 0; i < args->IP->getSize(); ++i){ // compute gradient of minus log-likelihood wrt phi
    const Individual* ind = args->IP->getIndividual(i);
    int Anc[2];
    ind->GetLocusAncestry(args->locus, Anc);
    logLikelihoodFirstDeriv(phi, Anc, ind->getPossibleHapPairs(args->locus), States, args->NumPops, dE_dphi);
  }
  //subtract derivative wrt phi of log prior density in softmax basis
  for(unsigned s = 0; s < States; ++s) {
    for(unsigned k = 0; k < args->NumPops; ++k) {
      dE_dphi[k*States+s]*=args->coolness;
    }
  }
  for(unsigned k = 0; k < args->NumPops; ++k) {
    for(unsigned s = 0; s < States; ++s) {
      if(ishapmixmodel) {
	dE_dphi[k*States+s] -= *(args->PriorParams) / phi[k*States+s];
      } else {
	if(args->PriorParams[s*args->NumPops +k] > 0.0)
	  dE_dphi[k*States+s] -= args->PriorParams[s*args->NumPops +k]  / phi[k*States+s]; 
      }
    }
  }
  
  //now use chain rule to obtain gradient wrt args: dE/dbeta_j = sum_i{dE/dphi_i dphi_i/dbeta_j}
  for(unsigned k = 0; k < args->NumPops; ++k) {
    double sum = 0.0; 
    for(unsigned s = 0; s < States; ++s) { 
      sum += exp(params[k*States+s]);
    }
    double sumsq = sum * sum;
    for(unsigned s = 0; s < States; ++s) { // loop over j
      for(unsigned s1 = 0; s1 < States; ++s1) { // loop over i to evaluate dE/dbeta_j as a sum
	if(s1 != s) {
	  g[k*States+s] -= dE_dphi[k*States+s1] * exp(params[k*States+s1]) * exp(params[k*States+s]) / sumsq; 
	} else {
	  g[k*States+s] += dE_dphi[k*States+s1] * (sum - exp(params[k*States+s])) * exp(params[k*States+s]) / sumsq; 
	}
      }
    }
  }
  delete[] phi;
  delete[] dE_dphi;
}

// requires: sampled ancestry pair A, PossibleHapPairs (compatible with genotype) H, 
// current values of AlleleFreqs at this locus, phi
//, number of alleles/haplotypes NumStates, number of populations, NumPops
double AlleleFreqSampler::logLikelihood(const double *phi, const int Anc[2], const std::vector<hapPair > H, 
					unsigned NumStates){
  double sum = 0.0;
  unsigned NumPossHapPairs = H.size();
  for(unsigned hpair = 0; hpair < NumPossHapPairs; ++hpair){
    unsigned j0 = H[hpair].haps[0];//j
    unsigned j1 = H[hpair].haps[1];//j'
    if( (Anc[0]==Anc[1]) || (j0==j1) ) {
      int index0 = Anc[0]*NumStates + j0;//jk
      int index1 = Anc[1]*NumStates + j1;//j'k'
      sum += log(phi[index0]*phi[index1]);
    } else {
      sum += log(phi[Anc[0]*NumStates + j1] * phi[Anc[1]*NumStates + j0] + 
		 phi[Anc[0]*NumStates + j0] * phi[Anc[1]*NumStates + j1] );
    }
  }
  return sum;
}

///first derivative of minus log likelihood wrt phi
void AlleleFreqSampler::logLikelihoodFirstDeriv(const double *phi, const int Anc[2], const std::vector<hapPair > H, 
						unsigned NumStates, unsigned, double* FirstDeriv){
  unsigned NumPossHapPairs = H.size();
  for(unsigned hpair = 0; hpair < NumPossHapPairs; ++hpair){
    unsigned j0 = H[hpair].haps[0];//j
    unsigned j1 = H[hpair].haps[1];//j'
    int index0 = Anc[0]*NumStates + j0;//jk
    int index1 = Anc[1]*NumStates + j1;//j'k'
    if( (Anc[0] == Anc[1]) || (index0 == index1) ) {
      FirstDeriv[index0] -= 1.0 / phi[index1];
      FirstDeriv[index1] -= 1.0 / phi[index0];
    } else {
      int index2 = Anc[0]*NumStates + j1;//jk'
      int index3 = Anc[1]*NumStates + j0;//j'k
      double denom = phi[index0] * phi[index1] + phi[index2] * phi[index3];
      FirstDeriv[index0] -= phi[index1] / denom;
      FirstDeriv[index1] -= phi[index0] / denom;
      FirstDeriv[index2] -= phi[index3] / denom;
      FirstDeriv[index3] -= phi[index2] / denom;
    }
  }
}

///energy function for Hamiltonian sampler, case of SNP
double AlleleFreqSampler::getEnergySNP(const double * const params, const void* const vargs){
  const AlleleFreqArgs* args = (const AlleleFreqArgs*)vargs;
  double energy = 0.0;
  unsigned Pops = args->NumPops;
  //transform params to freqs
  double* phi = new double[Pops];//phi[k] is freq of allele 1 in population k
  for(unsigned k = 0; k < Pops; ++k)
     phi[k] = 1.0 / (1.0+exp(-params[k]));//inverse-logit transformation
  //calculate minus loglikelihood
  for(unsigned k = 0; k < Pops; ++k){
    energy -= args->AlleleCounts[k]        * mylog(phi[k]);     /*allele1*/
    energy -= args->AlleleCounts[Pops + k] * mylog(1.0-phi[k]); /*allele2*/
    for(unsigned k1 = k+1; k1 < Pops; ++k1){ // loop over upper triangle of array hetCounts
      energy -= (args->hetCounts[k*Pops+k1] + args->hetCounts[k1*Pops+k]) * 
	mylog(phi[k] + phi[k1] - 2.0*phi[k]*phi[k1]);
    }
  }
  energy *= args->coolness;
  // subtract log prior density in logit basis
  for(unsigned k = 0; k < Pops; ++k){
    if(ishapmixmodel) { // prior is the same across populations and alleles
      energy -= Pops * *(args->PriorParams) * ( mylog(phi[k]) + mylog(1 - phi[k]) );
    } else {
      if(args->PriorParams[k] > 0.0)
	energy -= args->PriorParams[2*k]*mylog(phi[k]) + args->PriorParams[2*k+1]*mylog(1 - phi[k]);
    }
  }
  delete[] phi;
  return energy;
}

void AlleleFreqSampler::gradientSNP(const double * const params, const void* const vargs, double* g){
  const AlleleFreqArgs* args = (const AlleleFreqArgs*)vargs;
  unsigned Pops = args->NumPops;
  double* phi = new double[Pops];//phi[k] is freq of allele 1 in population k
  double* dE_dphi = new double[Pops];// derivative of energy wrt phi
  for(unsigned k = 0; k < Pops; ++k) {// transform params to freqs
    phi[k] = 1.0 / ( 1.0 + exp(-params[k]) );// inverse-logit transformation
  }  
  // calculate derivative of minus log-likelihood wrt phi
  for(unsigned k  = 0; k < Pops; ++k){
    dE_dphi[k] = -args->AlleleCounts[k] / phi[k];//allele1
    dE_dphi[k] += args->AlleleCounts[Pops + k] / (1.0 - phi[k]);//allele2
    for(unsigned k1 = k+1; k1 < Pops; ++k1) { // loop over upper triangle of array hetCounts
      double denom = phi[k] + phi[k1] - 2.0*phi[k]*phi[k1];
      dE_dphi[k]  -= (args->hetCounts[k*Pops+k1] + args->hetCounts[k1*Pops+k]) * (1.0 - 2.0*phi[k1])/denom;
      dE_dphi[k1] -= (args->hetCounts[k*Pops+k1] + args->hetCounts[k1*Pops+k]) * (1.0 - 2.0*phi[k]) /denom;
    }
  }
  for(unsigned k = 0; k < Pops; ++k){ // anneal likelihood 
    dE_dphi[k] *= args->coolness;
  }
  //subtract derivative of log prior density in logit basis wrt phi
  for(unsigned k = 0; k < Pops; ++k){
    dE_dphi[k] *= args->coolness;
    if(ishapmixmodel)
      dE_dphi[k] -= *(args->PriorParams) *( 1.0/phi[k] - 1.0/(1 - phi[k]) ); 
    else{
      if(args->PriorParams[k] > 0.0)
	dE_dphi[k] -= args->PriorParams[2*k]/phi[k] - args->PriorParams[2*k+1]/(1 - phi[k]); 
    }
  }
  //now use chain rule to obtain gradient wrt logit(phi): dphi/dtheta = phi * (1 - phi)
  for(unsigned k = 0; k < Pops; ++k){
    g[k] = dE_dphi[k] * phi[k] * (1.0 - phi[k]);
  }
  delete[] phi;
  delete[] dE_dphi;
}


