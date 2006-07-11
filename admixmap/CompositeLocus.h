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

///struct to hold a pair of haplotypes, coded as integers
typedef struct
 {
   int haps[2];
}hapPair  ; 

void operator<<(std::ostream& os, const hapPair &h);

///   Class to represent a composite locus
class CompositeLocus 
{

public:
  CompositeLocus();
  ~CompositeLocus();

  static void SetNumberOfPopulations( int );
  static void SetRandomAlleleFreqs(bool);
  void SetNumberOfStates( int );
  void SetLabel( int, std::string );
  void SetNumberOfLoci( int );
  void SetNumberOfAllelesOfLocus( int, int );
  void AddLocus( int, std::string );
  void SetHapPairProbs();
  void SetHapPairProbs(const double* alleleProbs);
  void InitialiseHapPairProbs(const double* const allelefreqs);
  void InitialiseHapPairProbsMAP();
  void setHaplotypeProbsMAP();
  void setAlleleProbsMAP(const double* const Freqs);
  void AccumulateAlleleProbs();

  int GetNumberOfLoci()const;
  int GetNumberOfStates()const;
  const std::string GetLabel(int)const;
  int GetNumberOfAllelesOfLocus( int )const;
  void getLocusAlleleProbs(double **P, int k)const;
  const double* getAlleleProbs()const{return AlleleProbs;}

  const std::vector<int> getAlleleCounts(int a, const int* happair)const;
  const std::vector<int> getHaplotypeCounts(const int* happair);
  void setPossibleHaplotypePairs(const std::vector<std::vector<unsigned short> > Genotype, std::vector<hapPair> &PossibleHapPairs);
  void setPossibleXHaplotypes(const std::vector<std::vector<unsigned short> > Genotype, std::vector<hapPair> &PossibleHapPairs);
  void decodeIntAsHapAlleles(const int h, int *hapAlleles)const;
  void GetGenotypeProbs(double *Probs, const std::vector<hapPair > &HaplotypePairs, bool chibindicator)const;
  void GetHaploidGenotypeProbs(double *Probs, const std::vector<hapPair > &HapPairs, bool chibindicator) const; 
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
  static int Populations;
  std::vector<int> NumberOfAlleles;
  const double *AlleleProbs;//< pointer to allele frequencies held in AlleleFreqs
  const double* AlleleProbsMAP;//< pointer to AlleleFreqsMAP held in AlleleFreqs
  double **SumAlleleProbs;//< sums of alleleprobs for a single population, used to compute loglikelihood at posterior means
#ifndef PARALLEL
  double *HapPairProbs; //< haplotype pair probabilities
  double *HapPairProbsMAP; //< Posterior estimates of hap pair probs
#endif
  std::vector<std::string> Label;
  int *base;
  static bool RandomAlleleFreqs;

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

  void setMissingAlleles(int numMissingLoci, int permMissing,  std::vector<int>& MissingAlleles, 
			 const std::vector<int>& MissingLoci); 
  void permuteMissingLoci(const std::vector<bool>& isMissing, const int numMissingLoci, const int permMissing, 
			  const std::vector<int>& HapAlleles,  const std::vector<int>& MissingLoci, 
			  std::vector<int>& HapAllelesNoMissing) ;

  // UNIMPLEMENTED
  // to avoid use
  CompositeLocus(const CompositeLocus&);
  CompositeLocus& operator=(const CompositeLocus&);
};

double GetMarginalLikelihood( const std::vector<double> PriorAlleleFreqs, const std::vector<int> AlleleCounts );

typedef std::vector<hapPair>::const_iterator happairiter;

#ifdef PARALLEL
inline void CompositeLocus::GetGenotypeProbs(double *Probs, const std::vector<hapPair > &HapPairs, bool /*chibindicator*/) const {
  double *q = Probs;
  const double* const p;
  if(!chibindicator || !RandomAlleleFreqs) 
    p = AlleleProbs;
  else 
    p = AlleleProbsMAP;

  happairiter end = HapPairs.end();
  for(int k0 = 0; k0 < Populations; ++k0)
    for(int k1 = 0; k1 < Populations; ++k1) {
    *q = 0.0;
    happairiter h = HapPairs.begin();
    for( ; h != end ; ++h) {
      //in parallel version, the happairprobs are not stored so we calculate them from allele probs
      *q += (p + k0*NumberOfStates+h->haps[0]) * (p + k1*NumberOfStates+h->haps[1]);
    }
    q++;
  }
}
#else
inline void CompositeLocus::GetGenotypeProbs(double *Probs, const std::vector<hapPair > &HapPairs, bool chibindicator) const {
  int Ksq = Populations*Populations;
  double *q = Probs;
  const double *p;
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
inline void CompositeLocus::GetHaploidGenotypeProbs(double *Probs, const std::vector<hapPair > &HapPairs, bool chibindicator) const {
  //TODO: check if this workd in parallel version
  double *q = Probs;
  const double *p;
  if(!chibindicator || !RandomAlleleFreqs) 
    p = AlleleProbs;
  else 
    p = AlleleProbsMAP;

  happairiter end = HapPairs.end();
  for(int k0 = 0; k0 < Populations; ++k0) {
    *q = 0.0;
    happairiter h = HapPairs.begin();
    for( ; h != end ; ++h) {
      *q += *(p + (h->haps[0] ) * Populations);
    }
    p++;
    q++;
  }
}

#endif /* !COMPOSITE_LOCUS_H */
