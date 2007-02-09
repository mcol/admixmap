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

using namespace::std;

///Class to represent an individual in a hapmixmodel
class HapMixIndividual : public Individual
{
public:
  HapMixIndividual();
  HapMixIndividual(int number, const Options* const options, const InputData* const Data);
  ~HapMixIndividual();
  static void SetStaticMembers(Genome* const pLoci, const Options* const options, const array_of_allelefreqs* const haploidGenotypeProbs, const double* const diploidGenotypeProbs);
  bool GenotypeIsMissing(unsigned int locus)const;
  bool simpleGenotypeIsMissing(unsigned locus)const;

private:

  vector<unsigned short> genotypes;//genotypes coded as single integers, assuming all typed loci are SNPs

  static const array_of_allelefreqs* HaploidGenotypeProbs;
  static const double* DiploidGenotypeProbs;

  void UpdateHMMInputs(unsigned int j, const Options* const options, 
                       const double* const , const vector<double> );

  void SetPossibleHaplotypePairs(unsigned locus, vector<unsigned short>::const_iterator g, vector<hapPair> &PossibleHapPairs );
};

#endif /* HAMPIXINDIVIDUAL_H */

