//=============================================================================
//
// Copyright (C) 2002-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
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
/// \file Haplotype.h
/// Definition of the Haplotype class.
//=============================================================================

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
  void setBaseMissing(const std::vector<int>& missingLoci,
                      int numMissingLoci, std::vector<int> baseMissing[2]);
  void setMissingAlleles(const vector<int> baseMissing[2], int numMissingLoci, int permMissing,
			 vector<int> MissingAlleles[2]);
  void setMissingAlleles(int numMissingLoci, int permMissing, vector<int>& MissingAlleles,
			 const vector<int>& MissingLoci);
  int codeHapAllelesAsInt(const int *hapAlleles);

  void codeHapAllelesPairAsIntPair(const vector<int> HapAllelesPair[2], int *hpair);

  void permuteHetLoci(const std::vector<bool>& isHet, int numHetLoci,
                      int permHet, bool *permbits, const GenotypeIterator* G,
                      std::vector<int> HapAllelesPair[2]);

  void permuteMissingLoci(const std::vector<bool>& isMissing,
                          int numMissingLoci, int permMissing,
                          const std::vector<int> HapAllelesPair[2],
                          const std::vector<int> baseMissing[2],
                          std::vector<int> HapAllelesPairNoMissing[2],
                          std::vector<int> MissingAlleles[2]);

  void permuteMissingLoci(const std::vector<bool>& isMissing,
                          int numMissingLoci, int permMissing,
                          const std::vector<int>& HapAlleles,
                          const std::vector<int>& MissingLoci,
                          std::vector<int>& HapAllelesNoMissing,
                          std::vector<int>& MissingAlleles);

  void setPossibleHaplotypes(int numMissingLoci, int numPermsMissing,
                             const std::vector<bool>& isMissing,
                             const std::vector<int>& HapAlleles,
                             std::vector<hapPair>& PossibleHapPairs);

};



/** @} */



#endif /* !HAPLOTYPE_H */
