// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   CompositeLocus.h 
 *   header file for CompositeLocus class
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef COMPOSITE_LOCUS_H
#define COMPOSITE_LOCUS_H 1

#include "rand.h"
#include "common.h"
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

typedef struct  
{
   int haps[2];
} hapPair; 

class CompositeLocus 
{


public:
  CompositeLocus();
  ~CompositeLocus();

  void setHaplotypeProbsMAP();
  void SetNumberOfLabels();
  void SetNumberOfPopulations( int );
  void SetRandomAlleleFreqs(bool);
  void SetNumberOfStates( int );
  void SetLabel( int, std::string );
  void SetNumberOfLoci( int );
  void SetNumberOfAllelesOfLocus( int, int );
  void AddLocus( int );
  void SetHapPairProbs();
  void SetHapPairProbs(const double* const* alleleProbs);
  void Initialise(const double* const allelefreqs);
  void SetAlleleProbs(const double* const alleleFreqs, bool);

  int GetNumberOfLoci()const;
  int GetNumberOfStates()const;
  const std::string GetLabel(int)const;
  int GetNumberOfAllelesOfLocus( int )const;
  void getLocusAlleleProbs(double **P, int k)const;

  void setPossibleHaplotypePairs(const unsigned short* const* Genotype, std::vector<hapPair> &PossibleHapPairs);
  void decodeIntAsHapAlleles(const int h, int *hapAlleles)const;
  void GetGenotypeProbs(double *Probs, const std::vector<hapPair > &HaplotypePairs, bool chibindicator)const;
  void SetHapPairProbsToPosteriorMeans(int iterations);
  void SampleHapPair(int hap[2], const std::vector<hapPair > &HapPairs, const int ancestry[2])const;

  //functions used for haplotype association score test 
  int GetMergedHaplotype( int i )const;
  int GetNumberOfMergedHaplotypes()const;
  const int *GetHapLabels( int ) const;
  void SetDefaultMergeHaplotypes( const double* const alpha);

private: 
  int NumberOfLoci;
  int NumberOfStates;
  int Populations;
  int *NumberOfAlleles;
  double **AlleleProbs;
  double **SumAlleleProbs;//sums of alleleprobs for a single population, used to compute loglikelihood at posterior means
  double *HapPairProbs; //haplotype pair probabilities
  double *HapPairProbsMAP; //Posterior estimates of hap pair probs
  std::string *Label;
  int *base;
  bool RandomAlleleFreqs;

  //possibly move out
  int *MergeHaplotypes;
  int *HapLabels;
  int NumberOfMergedHaplotypes;

  void SetNoMergeHaplotypes();

  void intToBits(int n, const int length, bool *bits) ;
  void setBaseForHapCode();
  void setBaseMissing(const std::vector<int> missingLoci, const int numMissingLoci, std::vector<int> baseMissing[2]);
  void setMissingAlleles(const std::vector<int> baseMissing[2], int numMissingLoci, int permMissing, std::vector<int> MissingAlleles[2]) ;
  int codeHapAllelesAsInt(const int *hapAlleles);
  void codeHapAllelesPairAsIntPair(const std::vector<int> HapAllelesPair[2], int *hpair);
  void permuteHetLoci(const std::vector<bool> isHet, const int numHetLoci, const int permHet, 
		      const unsigned short* const* Genotype, std::vector<int> HapAllelesPair[2]);
  void permuteMissingLoci(const std::vector<bool> isMissing, const int numMissingLoci, const int permMissing, 
			  const std::vector<int> HapAllelesPair[2], const std::vector<int> baseMissing[2], std::vector<int> HapAllelesPairNoMissing[2]) ;
  // UNIMPLEMENTED
  // to avoid use
  CompositeLocus(const CompositeLocus&);
  CompositeLocus& operator=(const CompositeLocus&);
};

double GetMarginalLikelihood( const std::vector<double> PriorAlleleFreqs, const std::vector<int> AlleleCounts );

inline void CompositeLocus::GetGenotypeProbs(double *Probs, const std::vector<hapPair > &HapPairs, bool chibindicator)const{
  for(int k0 = 0; k0 < Populations * Populations; ++k0){
    Probs[k0] = 0.0;
    for(unsigned int h = 0; h < HapPairs.size() ; ++h)
      if(RandomAlleleFreqs && chibindicator )
	Probs[k0] += HapPairProbsMAP[HapPairs[h].haps[0] * NumberOfStates * Populations * Populations +
				     HapPairs[h].haps[1] * Populations * Populations +
				     k0];
      else
	Probs[k0] += HapPairProbs[HapPairs[h].haps[0] * NumberOfStates * Populations * Populations +
				  HapPairs[h].haps[1] * Populations * Populations +
				  k0];
  }
}


#endif /* !COMPOSITE_LOCUS_H */
