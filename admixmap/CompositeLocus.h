// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   CompositeLocus.h 
 *   header file for CompositeLocus class
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef COMPOSITE_LOCUS_H
#define COMPOSITE_LOCUS_H 1

#include "common.h"

typedef struct
 {
   int haps[2];
}hapPair  ; 

class CompositeLocus 
{

public:
  CompositeLocus();
  ~CompositeLocus();

  void setHaplotypeProbsMAP();
  void SetNumberOfPopulations( int );
  void SetRandomAlleleFreqs(bool);
  void SetNumberOfStates( int );
  void SetLabel( int, std::string );
  void SetNumberOfLoci( int );
  void SetNumberOfAllelesOfLocus( int, int );
  void AddLocus( int, std::string );
  void SetHapPairProbs();
  void SetHapPairProbs(const double* alleleProbs);
  void InitialiseHapPairProbs(const double* const allelefreqs);
  void InitialiseHapPairProbsMAP();
  void AccumulateAlleleProbs();

  int GetNumberOfLoci()const;
  int GetNumberOfStates()const;
  const std::string GetLabel(int)const;
  int GetNumberOfAllelesOfLocus( int )const;
  void getLocusAlleleProbs(double **P, int k)const;

  const std::vector<int> getAlleleCounts(int a, const int* happair)const;
  const std::vector<int> getHaplotypeCounts(const int* happair);
  void setPossibleHaplotypePairs(const std::vector<std::vector<unsigned short> > Genotype, std::vector<hapPair> &PossibleHapPairs);
  void decodeIntAsHapAlleles(const int h, int *hapAlleles)const;
  void GetGenotypeProbs(double *Probs, const std::vector<hapPair > &HaplotypePairs, bool chibindicator)const;
  void SetHapPairProbsToPosteriorMeans(int iterations);
  void SampleHapPair(hapPair*, const std::vector<hapPair > &HapPairs, const int ancestry[2])const;

  //functions used for haplotype association score test 
  int GetMergedHaplotype( int i )const;
  int GetNumberOfMergedHaplotypes()const;
  const int *GetHapLabels( int ) const;
  void SetDefaultMergeHaplotypes( const double* const alpha);

private: 
  int NumberOfLoci;
  int NumberOfStates;
  int Populations;
  std::vector<int> NumberOfAlleles;
  const double *AlleleProbs;//pointer to allele frequencies held in AlleleFreqs
  double **SumAlleleProbs;//sums of alleleprobs for a single population, used to compute loglikelihood at posterior means
#ifndef PARALLEL
  double *HapPairProbs; //haplotype pair probabilities
  double *HapPairProbsMAP; //Posterior estimates of hap pair probs
#endif
  std::vector<std::string> Label;
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
		      const std::vector<std::vector<unsigned short> > Genotype, std::vector<int> HapAllelesPair[2]);
  void permuteMissingLoci(const std::vector<bool> isMissing, const int numMissingLoci, const int permMissing, 
			  const std::vector<int> HapAllelesPair[2], const std::vector<int> baseMissing[2], std::vector<int> HapAllelesPairNoMissing[2]) ;
  // UNIMPLEMENTED
  // to avoid use
  CompositeLocus(const CompositeLocus&);
  CompositeLocus& operator=(const CompositeLocus&);
};

double GetMarginalLikelihood( const std::vector<double> PriorAlleleFreqs, const std::vector<int> AlleleCounts );

typedef std::vector<hapPair>::const_iterator happairiter;

#ifdef PARALLEL
inline void CompositeLocus::GetGenotypeProbs(double *Probs, const std::vector<hapPair > &HapPairs, bool chibindicator) const {
  //TODO: happairprobsMAP
  double *q = Probs;

  happairiter end = HapPairs.end();
  for(int k0 = 0; k0 < Populations; ++k0)
    for(int k1 = 0; k1 < Populations; ++k1) {
    *q = 0.0;
    happairiter h = HapPairs.begin();
    for( ; h != end ; ++h) {
      //in parallel version, the happairprobs are not stored so we calculate them from allele probs
      *q += AlleleProbs[k0*NumberOfStates+h->haps[0]] * AlleleProbs[k1*NumberOfStates+h->haps[1]];
    }
    q++;
  }
}
#else
inline void CompositeLocus::GetGenotypeProbs(double *Probs, const std::vector<hapPair > &HapPairs, bool chibindicator) const {
  int Ksq = Populations*Populations;
  double *q = Probs;
  double *p;
  if(!chibindicator || !RandomAlleleFreqs) 
    p = HapPairProbs;
  else 
    p = HapPairProbsMAP;

  happairiter end = HapPairs.end();
  for(int k0 = 0; k0 < Ksq; ++k0) {
    *q = 0.0;
    happairiter h = HapPairs.begin();
    for( ; h != end ; ++h) {
      *q += *(p + (h->haps[0] * NumberOfStates + h->haps[1]) * Ksq);

    }
    p++;
    q++;
  }
}
#endif

#endif /* !COMPOSITE_LOCUS_H */
