// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   AlleleFreqSampler.h 
 *   Class to sample allele/haplotype frequencies at a given composite locus.
 *   Copyright (c) 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef ALLELEFREQSAMPLER_H
#define ALLELEFREQSAMPLER_H
#include "common.h"
#include <math.h>
#include <stdlib.h>
#include "functions.h"
#include "HamiltonianMonteCarlo.h"
#include "CompositeLocus.h"
#include <algorithm>
#include <gsl/gsl_linalg.h>

class IndividualCollection;
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
  const int* AlleleCounts;
  const int* hetCounts;
  double coolness;
}AlleleFreqArgs;

class AlleleFreqSampler{
public:
  AlleleFreqSampler();
  void SampleAlleleFreqs(double *phi, const double* Prior, IndividualCollection* IC, unsigned locus, 
			 unsigned NumStates, unsigned NumPops, double coolness);
  void SampleSNPFreqs(double *phi, const double* Prior, const int* AlleleCounts, const int* hetCounts, unsigned locus, 
			 unsigned NumPops, double coolness);


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
  static double getEnergySNP(const double * const phi, const void* const vargs);
  static void gradientSNP(const double * const phi, const void* const vargs, double* g);
};
#endif
