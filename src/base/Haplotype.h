// *-*-C++-*-*
/**
 *   CompositeLocus.h
 *   header file for CompositeLocus class
 */

/*
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *
 * This program is free software distributed WITHOUT ANY WARRANTY.
 * You can redistribute it and/or modify it under the terms of the GNU General Public License,
 * version 2 or later, as published by the Free Software Foundation.
 * See the file COPYING for details.
 *
 */

#ifndef HAPLOTYPE_H
#define HAPLOTYPE_H 1


#include <vector>
#include "HapPair.h"
class GenotypeIterator;


using namespace::std;


/** \addtogroup base
 * @{ */


class Haplotype{
public:
  Haplotype();
  ~Haplotype();
  void AddLocus(int nalleles);

  void setPossibleHaplotypePairs(const GenotypeIterator* G, vector<hapPair> &PossibleHapPairs);

  void setPossibleHaplotypes(const GenotypeIterator* G, vector<hapPair> &PossibleHapPairs);

  void decodeIntAsHapAlleles(const int h, int *hapAlleles)const;

private:
  int NumberOfLoci;
  vector<int> NumberOfAlleles;
  int *base;

  void intToBits(int n, const int length, bool *bits);
  void setBaseForHapCode();
  void setBaseMissing(const vector<int> missingLoci, const int numMissingLoci, vector<int> baseMissing[2]);
  void setMissingAlleles(const vector<int> baseMissing[2], int numMissingLoci, int permMissing,
			 vector<int> MissingAlleles[2]);
  void setMissingAlleles(int numMissingLoci, int permMissing, vector<int>& MissingAlleles,
			 const vector<int>& MissingLoci);
  int codeHapAllelesAsInt(const int *hapAlleles);

  void codeHapAllelesPairAsIntPair(const vector<int> HapAllelesPair[2], int *hpair);

  void permuteHetLoci(const vector<bool> isHet, const int numHetLoci, const int permHet,
                      const GenotypeIterator* G, vector<int> HapAllelesPair[2]);
  void permuteMissingLoci(const vector<bool> isMissing, const int numMissingLoci, const int permMissing,
			  const vector<int> HapAllelesPair[2], const vector<int> baseMissing[2],
			  vector<int> HapAllelesPairNoMissing[2]);
  void permuteMissingLoci(const vector<bool>& isMissing, const int numMissingLoci, const int permMissing,
			  const vector<int>& HapAlleles, const vector<int>& MissingLoci,
			  vector<int>& HapAllelesNoMissing);

  void setPossibleHaplotypes(int numMissingLoci, int numPermsMissing, const vector<bool>& isMissing,
			     const vector<int>& HapAlleles, vector<int>& HapAllelesNoMissing,
			     vector<hapPair> &PossibleHapPairs);

};



/** @} */



#endif /* !HAPLOTYPE_H */
