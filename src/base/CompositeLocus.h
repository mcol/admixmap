// *-*-C++-*-*
/** 
 *   CompositeLocus.h 
 *   header file for CompositeLocus class
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
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
#include "HapPair.h"
#include "Haplotype.h"
#include "bclib/pvector.h"

///   Class to represent a composite locus
class CompositeLocus
{

public:
  CompositeLocus();
  virtual ~CompositeLocus();

  static void SetNumberOfPopulations( int );
  static void SetRandomAlleleFreqs(bool);
  void SetNumberOfStates( int );
  void SetLabel( int, std::string );
  void SetNumberOfLoci( int );
  void SetNumberOfAllelesOfLocus( int, int );
  void AddLocus( int, std::string );
  void SetHapPairProbs();
  void InitialiseHapPairProbs(const double* const allelefreqs, bool AllHaploid);
  void InitialiseHapPairProbsMAP();
  void SetHapPairProbsMAP();
  void setAlleleProbsMAP(const double* const Freqs);
  void AccumulateAlleleProbs();

  void getConditionalHapPairProbs(bclib::pvector<double>& Probs, const std::vector<hapPair > &PossibleHapPairs, const int ancestry[2])const;
  /** Similar to getConditionalHapPairProbs, but only the first and last
   * elements of Probs are being set.
   * 
   * It's a shortcut to speed it up in a special case when a genotype
   * is missing. 
   */
  double getFirstConditionalHapPairProbs(const int ancestry[2]) const;
  double getLastConditionalHapPairProbs(const int ancestry[2]) const;
  void SampleHapPair(hapPair*, const std::vector<hapPair > &PossibleHapPairs, const int ancestry[2])const;


  int GetNumberOfLoci()const;
  int GetNumberOfStates()const;
  const std::string GetLabel(int)const;
  int GetNumberOfAllelesOfLocus( int )const;
  void getLocusAlleleProbs(double **P, int k)const;
  const double* getAlleleProbs()const{return AlleleProbs;}

  const std::vector<int> getAlleleCounts(int a, const int* happair)const;
  const std::vector<int> getHaplotypeCounts(const int* happair);

  void GetGenotypeProbs(double *Probs, const std::vector<hapPair > &HaplotypePairs, bool chibindicator)const;
  void GetHaploidGenotypeProbs(double *Probs, const std::vector<hapPair > &HapPairs, bool chibindicator) const; 
  void SetHapPairProbsToPosteriorMeans(int iterations);
  const int *GetHapLabels( int ) const;

  //functions used for haplotype association score test 
  int GetMergedHaplotype( int i )const;
  int GetNumberOfMergedHaplotypes()const;
  void SetDefaultMergeHaplotypes( const double* const alpha);

  Haplotype HaplotypeSetter;

protected:
  int NumberOfStates;
  static int Populations;
  static int PopulationsSquared;
  static int PopulationsSquared_x_3; //< Caching a value for efficiency
  double *HapPairProbs; //< haplotype pair probabilities calculated using AlleleProbs

private:
  int NumberOfLoci;
  /// Squared number of populations, stored for efficiency reasons.
  std::vector<int> NumberOfAlleles;
  const double *AlleleProbs;//< pointer to allele frequencies held in AlleleFreqs
  const double *AlleleProbsMAP;//< pointer to AlleleFreqsMAP held in AlleleFreqs
  double *SumAlleleProbs;//< sums of alleleprobs for a single population, used to compute loglikelihood at posterior means
  double *HapPairProbsMAP; //< hap pair probs calculated using AlleleProbsMAP
  std::vector<std::string> Label;
  int *base;
  static bool RandomAlleleFreqs;

  //possibly move out
  int *MergeHaplotypes;
  int *HapLabels;
  int NumberOfMergedHaplotypes;

  void SetNoMergeHaplotypes();


  void SetHapPairProbs(const double* alleleProbs, double* happairprobs);

  // UNIMPLEMENTED
  // to avoid use
  CompositeLocus(const CompositeLocus&);
  CompositeLocus& operator=(const CompositeLocus&);
};

double GetMarginalLikelihood( const std::vector<double> PriorAlleleFreqs, const std::vector<int> AlleleCounts );

typedef std::vector<hapPair>::const_iterator happairiter;

inline void CompositeLocus::GetGenotypeProbs(double *Probs, const std::vector<hapPair > &HapPairs, 
					     bool chibindicator) const {
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

inline void CompositeLocus::GetHaploidGenotypeProbs(double *Probs, const std::vector<hapPair > &HapPairs, bool chibindicator) const {
  double *q = Probs;
  const double *p;
  if(!chibindicator || !RandomAlleleFreqs) 
    p = AlleleProbs;
  else 
    p = AlleleProbsMAP;

  happairiter end = HapPairs.end();
  for(int k = 0; k < Populations; ++k) {
    *q = 0.0;
    happairiter h = HapPairs.begin();
    for( ; h != end ; ++h) {
//Probs[k] += p[k*NumberOfStates + (h->haps[0] )];
	*q += p[k*NumberOfStates + (h->haps[0] )];
    }
    q++;
  }
}



#endif /* !COMPOSITE_LOCUS_H */
