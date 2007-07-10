// *-*-C++-*-*
/** 
 *   HAPMIXMAP
 *   HapMixIndividual.h 
 *   header file for HapMixIndividual class
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef HAPMIXINDIVIDUAL_H
#define HAPMIXINDIVIDUAL_H 1
#include "Individual.h"
#include "HapMixGenome.hh"
#include "GPI.h"
#include <map>
#include "bcppcl/pvector.h"

using namespace std;
class InputHapMixData;

///Class to represent an individual in a hapmixmodel
class HapMixIndividual : public Individual
{
public:
  HapMixIndividual();
  HapMixIndividual(int number, const Options* const options, InputHapMixData* const Data, const double* GlobalTheta, bool isMasked);
  ~HapMixIndividual();
  static void SetGenotypeProbs(HapMixGenome* const G, const FreqArray& haploidGenotypeProbs, const FreqArray& diploidGenotypeProbs);
  void SetMissingGenotypes();
  bool GenotypeIsMissing(unsigned int locus)const;
  bool simpleGenotypeIsMissing(unsigned locus)const;
  void SampleJumpIndicators(int* SumArrivalCounts);

  void AccumulateConcordanceCounts(int* ConcordanceCounts)const;
  void calculateUnorderedGenotypeProbs(const Options& options);
  const bcppcl::pvector<double>& getStateProbs(const bool, int, int) const;
  void calculateUnorderedGenotypeProbs(unsigned);
  const std::vector<std::vector<double> >& getUnorderedProbs(const unsigned int) const;

private:
  vector<unsigned short> hgenotypes;//< genotypes coded as single integers, assuming all typed loci are SNPs

  /// Temporary vector used by calculateUnorderedProbs
  bcppcl::pvector<double> orderedGenotypeProbs;

  //shared pointer to HapMixGenome object
  static HapMixGenome* pG;

  static const FreqArray* HaploidGenotypeProbs;
  static const FreqArray* DiploidGenotypeProbs;
  static GenotypeProbIterator GPI;
  
  std::vector<std::vector<std::vector<double> > > UnorderedProbs;
  static map<int, int> ord2unord; //< Map ordered to unordered probabilities

  void UpdateHMMInputs(unsigned int j, const Options& options, 
                     const double* const , const vector<double> );

  //  void SetPossibleHaplotypePairs(unsigned locus, vector<unsigned short>::const_iterator g, vector<hapPair> &PossibleHapPairs );
  void SetPossibleHaplotypePairs(const vector<unsigned short>::const_iterator GI, vector<hapPair> &PossibleHapPairs, ploidy p);
};

#endif /* HAMPIXINDIVIDUAL_H */

