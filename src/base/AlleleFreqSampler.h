//=============================================================================
//
// Copyright (C) 2006  David O'Donnell, Clive Hoggart and Paul McKeigue
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file AlleleFreqSampler.h
/// Definition of the AlleleFreqSampler class.
//=============================================================================

#ifndef ALLELEFREQSAMPLER_H
#define ALLELEFREQSAMPLER_H


#include <gsl/gsl_linalg.h>
#include "bclib/HamiltonianMonteCarlo.h"
#include "HapPair.h"


/** \addtogroup base
  @{ */



class IndividualCollection;

/// Arguments for sampling allele freqs
struct AlleleFreqArgs
  {
  unsigned			NumPops	    ; ///< #populations
  unsigned			NumStates   ; ///< #alleles/haplotypes
  unsigned			locus	    ; ///< current locus
  const IndividualCollection *	IP	    ; ///< pointer to individuals
  const double *		PriorParams ; ///< parameters of Dirichlet prior on allele freqs
  const int *			AlleleCounts;
  const int *			hetCounts   ;
  double			coolness    ;
  };



///Class to sample allele/haplotype frequencies at a given composite locus.
class AlleleFreqSampler{
public:
  AlleleFreqSampler();
  ~AlleleFreqSampler();
  AlleleFreqSampler(unsigned NumStates, unsigned NumPops, const double* const Prior, bool);
  void SampleAlleleFreqs(double *phi, IndividualCollection* IC, unsigned locus, 
			 unsigned NumStates, unsigned NumPops, double coolness);
  void SampleSNPFreqs(double *phi, const int* AlleleCounts, const int* hetCounts, unsigned locus, 
			 unsigned NumPops, double coolness);
  void resetStepSizeApproximator(int k);
  double getStepSize()const;
  double getAcceptanceRate()const;

private:
  bool samplerInitialized; ///< Flag to tell us if (NumStates==2)&&(NumPops<=1), i.e. sampler is uninitialized.
  bclib::HamiltonianMonteCarlo Sampler;
  AlleleFreqArgs Args;
  double* params;
  static bool ishapmixmodel;

  static double logLikelihood(const double *phi, const int Anc[2], const std::vector<hapPair> & H, const unsigned NumPops);
  static double logPrior(const double* PriorParams, const double* phi, const unsigned NumPops, const unsigned NumStates);
  static double logJacobian(const double* a, const double z, const unsigned H);
  static void minusLogLikelihoodFirstDeriv(const double *phi, const int Anc[2], const std::vector<hapPair > & H,
					   const unsigned NumStates, double* FirstDeriv);
  static double getEnergy(const double * const phi, const void* const vargs);
  static void gradient(const double * const phi, const void* const vargs, double* g);
  static double getEnergySNP(const double * const phi, const void* const vargs);
  static void gradientSNP(const double * const phi, const void* const vargs, double* g);

  AlleleFreqSampler(const AlleleFreqSampler&);
  void operator=(const AlleleFreqSampler&);
};


/**@}*/


#endif
