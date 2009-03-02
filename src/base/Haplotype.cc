/** 
 *   Haplotypes.cc 
 *   Class to set Possible Haplotypes/Haplotype Pairs
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "Haplotype.h"

#include "GenotypeIterator.h" // Used forward reference in header

Haplotype::Haplotype(){
  base = NULL;
  NumberOfLoci = 0;
}
Haplotype::~Haplotype(){
  delete[] base;
}

// void Haplotype::Initialise(const vector<int>& NumAlleles){
//   NumberOfLoci = NumAlleles.size();
//   NumberOfAlleles = NumAlleles;
// }
void Haplotype::AddLocus(int nalleles){
  ++NumberOfLoci;
  NumberOfAlleles.push_back(nalleles);
}

/*
********************************************************
Diploid functions for setting PossibleHaplotypePairs
********************************************************
*/


/// updates: PossibleHapPairs, an stl vector of arrays of 2 integers.
/// arguments: genotype is 2D array of alleles with 2 rows and NumberOfLoci cols. 
/// call once for each individual at start of program 
void Haplotype::setPossibleHaplotypePairs(const GenotypeIterator* G, vector<hapPair> &PossibleHapPairs){
  setBaseForHapCode();
  PossibleHapPairs.clear();
  if(!G->isDiploid()){//haploid genotypes
    setPossibleHaplotypes(G, PossibleHapPairs);
    return;
  }

  int numHetLoci = 0;
  int numMissingLoci = 0;
  int numPermsHet = 1;
  int numPermsMissing = 1;
  vector<bool> isHet(NumberOfLoci);
  vector<bool> isMissing(NumberOfLoci);
  //loop over loci to determine which are heterozygous or missing
  for( int i = 0; i < NumberOfLoci; i++ ) {
    isHet[i] = false;
    isMissing[i] = false;
    if( (G->get(i, 0) == 0) || (G->get(i,1) == 0)  ) { // missing genotype
      isMissing[i] = true;
      numMissingLoci ++;
      numPermsMissing *= NumberOfAlleles[i] * NumberOfAlleles[i];
    } else {
      if( G->get(i, 0) !=  G->get(i, 1) ) { // heterozygous
	isHet[i] = true; 
	numHetLoci++;
	numPermsHet *= 2;
      }
    }
  }
  //   cout << numMissingLoci << " missing genotypes with " << numPermsMissing << " permutations\n";
  //   cout << numHetLoci << " heterozygous genotypes with " << numPermsHet << " permutations\n";
  //  setPossibleHaplotypePairs(numHetLoci, numMissingLoci, numPermsHet, numPermsMissing, isHet, isMissing);
  //}

  //void Haplotype::setPossibleHaplotypePairs(int numHetLoci, int numMissingLoci, int numPermsHet, int numPermsMissing, 
  //				       const vector<bool>& isHet, const vector<bool>& isMissing){
  // loop over loci to calculate vectors hetLoci and missingLoci
  vector<int> hetLoci(numHetLoci);
  vector<int> missingLoci(numMissingLoci);
  int offsetHetLoci = 0;
  int offsetMissingLoci = 0;
  for( int i = 0; i < NumberOfLoci; i++ ) {
    if( isHet[i] ) {
      hetLoci[offsetHetLoci] = i;
      offsetHetLoci++;
    }
    if( isMissing[i] ) {
      missingLoci[offsetMissingLoci] = i;
      offsetMissingLoci++;
    }
  }
  
  vector<int> HapAllelesPair[2] = {vector<int>(NumberOfLoci), vector<int>(NumberOfLoci)};
  vector<int> HapAllelesPairNoMissing[2]= {vector<int>(NumberOfLoci), vector<int>(NumberOfLoci)};;
  hapPair hpair;
  // loop over all possible ordered haplotype pairs compatible with unphased genotype
  // loop over permsHet to permute alleles at heterozygous loci
  for (int permHet = 0; permHet < numPermsHet; permHet++ ) {
    permuteHetLoci(isHet, numHetLoci, permHet, G, HapAllelesPair);
    if( numMissingLoci == 0 ) {
      codeHapAllelesPairAsIntPair(HapAllelesPair, hpair.haps);
      PossibleHapPairs.push_back(hpair);
    } else {
      vector<int> baseMissing[2] = {vector<int>(numMissingLoci), vector<int>(numMissingLoci)};
      setBaseMissing(missingLoci, numMissingLoci, baseMissing);
      // loop over permsMissing to permute alleles at missing loci
      for ( int permMissing = 0; permMissing < numPermsMissing; permMissing++ ) {
	permuteMissingLoci(isMissing, numMissingLoci, permMissing, HapAllelesPair, baseMissing, 
			   HapAllelesPairNoMissing);
	codeHapAllelesPairAsIntPair(HapAllelesPairNoMissing, hpair.haps);
	PossibleHapPairs.push_back(hpair);
      }
    }
  }
}

/*
********************************************************
Haploid functions for setting PossibleHaplotypePairs
********************************************************
*/

void Haplotype::setPossibleHaplotypes(const GenotypeIterator* G, vector<hapPair> &PossibleHapPairs){

  //setBaseForHapCode();

  int numMissingLoci  = 0;
  int numPermsMissing = 1;
  vector<bool>	isMissing	    ( NumberOfLoci );
  vector<int>	HapAlleles	    ( NumberOfLoci );
  vector<int>	HapAllelesNoMissing ( NumberOfLoci );

  for ( int i = 0; i < NumberOfLoci; i++ ) {
    unsigned short val = G->get(i, 0);
    HapAlleles[i] = val;
    if ( isMissing[i] = (val == 0) ) { // missing genotype
	numMissingLoci ++;
	numPermsMissing *= NumberOfAlleles[i];
    }
  }

  setPossibleHaplotypes(numMissingLoci, numPermsMissing, isMissing, HapAlleles, HapAllelesNoMissing, PossibleHapPairs);
}

void Haplotype::setPossibleHaplotypes(int numMissingLoci, int numPermsMissing, const vector<bool>& isMissing, 
				      const vector<int>& HapAlleles, vector<int>& HapAllelesNoMissing, 
				      vector<hapPair> &PossibleHapPairs){
  vector<int> missingLoci(numMissingLoci);
  int offsetMissingLoci = 0;
  for( int i = 0; i < NumberOfLoci; i++ ) {
    if( isMissing[i] ) {
      missingLoci[offsetMissingLoci] = i;
      offsetMissingLoci++;
    }
  }

  hapPair hpair;
  hpair.haps[1] = -1;//using -1 to denote *haplotype* rather than happair (must be negative as alleles are counted from 0)

  if( numMissingLoci == 0 ) {
    hpair.haps[0] = codeHapAllelesAsInt(&(HapAlleles[0]));
    PossibleHapPairs.push_back(hpair);
  } else {
    // loop over permsMissing to permute alleles at missing loci
    for ( int permMissing = 0; permMissing < numPermsMissing; permMissing++ ) {
      permuteMissingLoci(isMissing, numMissingLoci, permMissing, HapAlleles, missingLoci,
			 HapAllelesNoMissing);
      hpair.haps[0] = codeHapAllelesAsInt(&(HapAllelesNoMissing[0]));
      PossibleHapPairs.push_back(hpair); 
    }
    
  }

}

/*
*****************************************************************
Auxiliary functions follow
****************************************************************
*/

// arguments: integer, length of bit array
// returns: 1D array of bits representing integer
void Haplotype::intToBits(int n, const int length, bool *bits){
  for( int i = length - 1; i >= 0; i-- ) {
    bits[i] = n % 2;
    // right shift n by 1 to get next bit 
    n >>= 1; 
  }
}

// updates 1D array of counting bases used to increment hap code
// should call once to initialize CompositeLocus object 
void Haplotype::setBaseForHapCode(){
  delete[] base;
  base = new int[NumberOfLoci];
  base[NumberOfLoci - 1] = 1;
  for( int i = NumberOfLoci - 2; i >= 0; i-- ) {
    base[i] = base[i + 1] * NumberOfAlleles[i + 1];
  }
}

/// updates 2D array of counting bases used to increment permMissing
void Haplotype::setBaseMissing(const vector<int> missingLoci, const int numMissingLoci, vector<int> baseMissing[2]){
  baseMissing[1][numMissingLoci - 1] = 1;
  baseMissing[0][numMissingLoci - 1] = NumberOfAlleles[ missingLoci[numMissingLoci - 1] ];
  if( numMissingLoci > 1 ) {
    for( int i = numMissingLoci - 2; i >= 0; i-- ) {
      baseMissing[1][i] = baseMissing[0][i + 1] * NumberOfAlleles[ missingLoci[i + 1] ]; 
      baseMissing[0][i] = baseMissing[1][i] * NumberOfAlleles[ missingLoci[i] ]; 
    }
  }
}

/// updates 2D array MissingAlleles with alleles specified by permMissing
void Haplotype::setMissingAlleles(const vector<int> baseMissing[2], int numMissingLoci, int permMissing, 
				       vector<int> MissingAlleles[2]){
  int remainder = permMissing; 
  for( int i = 0; i < numMissingLoci; i++ ) { // loop over loci and gametes to extract alleles by integer division
    MissingAlleles[0][i] = 1 + remainder / baseMissing[0][i]; 
    remainder = remainder % baseMissing[0][i]; 
    MissingAlleles[1][i] = 1 + remainder / baseMissing[1][i]; 
    remainder = remainder % baseMissing[1][i]; 
  }
  // uncomment to display haplotype alleles  
  //   cout << "\npermMissing " << permMissing; 
  //   for( int col = 0; col < 2; col++ ) {
  //     cout << "\n";
  //     for( int row = 0; row < numMissingLoci; row++ ) {
  //       cout << MissingAlleles[col][row] << " ";
  //     }
  //   }
  //   cout << "\n";
}

void Haplotype::setMissingAlleles(int numMissingLoci, int permMissing, vector<int>& MissingAlleles, 
				       const vector<int>& MissingLoci){
  int remainder = permMissing; 
  vector<int> base(numMissingLoci);
  base[numMissingLoci-1] = 1;
  for( int i = numMissingLoci-2; i >=0 ; --i ){
    base[i] = NumberOfAlleles[MissingLoci[i+1]];
  }

  for( int i = 0; i < numMissingLoci; i++ ) { // loop over loci and gametes to extract alleles by integer division
    MissingAlleles[i] = 1 + remainder / base[i]; 
    remainder = remainder % NumberOfAlleles[MissingLoci[i]]; 
  }
}
/// codes haplotypes starting at 0 by incrementing counter from right to left starting at 0. 
/// arguments: NumberOfAlleles, base, and hapAlleles are 1D arrays of length NumberOfLoci
/// returns: hap code as integer
int Haplotype::codeHapAllelesAsInt(const int *hapAlleles){
  int h = hapAlleles[NumberOfLoci - 1] - 1;
  for(int i = NumberOfLoci - 2; i >= 0; i-- ) {
    h += (hapAlleles[i] - 1) * base[i];
  }
  return( h );
}

/// decode hap code to array of integers
/// updates: 1D array hapAlleles
void Haplotype::decodeIntAsHapAlleles(const int h, int *hapAlleles)const{
  int remainder = h;
  // loop over loci to extract alleles by integer division
  for(int i = 0; i < NumberOfLoci; i++ ) {
    hapAlleles[i] = 1 + remainder / base[i];
    remainder = remainder % base[i];
  }
}

/// updates 1D array of 2 integers coding haplotypes  
/// arguments: HapAllelesPair is 2D array with NumberOfLoci rows and 2 cols
void Haplotype::codeHapAllelesPairAsIntPair(const vector<int> HapAllelesPair[2], int *hpair){
  hpair[0] = HapAllelesPair[0][NumberOfLoci - 1]- 1;
  hpair[1] = HapAllelesPair[1][NumberOfLoci - 1] - 1;
  for( int i = NumberOfLoci - 2; i >= 0; i-- ) {
    hpair[0] += (HapAllelesPair[0][i] - 1) * base[i];
    hpair[1] += (HapAllelesPair[1][i] - 1) * base[i];
  }
}

void Haplotype::permuteHetLoci(const vector<bool> isHet, const int numHetLoci, const int permHet, 
			       const GenotypeIterator* G, vector<int> HapAllelesPair[2]){
  //recode permHet as array of bits, with length equal to NumHetLoci
  bool* permbits = new bool[numHetLoci];
  intToBits(permHet, numHetLoci, permbits);
  
  int hetLocusOffset = 0; // used to loop over heterozygous loci 
  //loop over loci to assign HapAllelesPair and swap alleles where element of permbits is 1;
  for(int locus = 0; locus < NumberOfLoci; locus++ ) {
    if( !isHet[locus] || !permbits[hetLocusOffset] ) { // assign alleles without swapping
      HapAllelesPair[0][locus] = G->get(locus, 0);
      HapAllelesPair[1][locus] = G->get(locus, 1);
    } else {
      if( permbits[hetLocusOffset] ) { // swap alleles
	HapAllelesPair[0][locus] = G->get(locus, 1);
	HapAllelesPair[1][locus] = G->get(locus, 0);
      }
    }
    hetLocusOffset += isHet[locus]; 
  } 
  delete[] permbits;
}

void Haplotype::permuteMissingLoci(const vector<bool> isMissing, const int numMissingLoci, const int permMissing, 
			const vector<int> HapAllelesPair[2], const vector<int> baseMissing[2], vector<int> HapAllelesPairNoMissing[2]) {
  vector<int> MissingAlleles[2] = {vector<int>(numMissingLoci), vector<int>(numMissingLoci)};
  setMissingAlleles(baseMissing, numMissingLoci, permMissing, MissingAlleles);
  int missingLocusOffset = 0; // used to loop over missing loci 
  //loop over loci to assign HapAllelesPairNoMissing
  for(int locus = 0; locus < NumberOfLoci; locus++ ) {
    if( ! isMissing[locus] ) { // assign alleles from HapAllelesPair
      HapAllelesPairNoMissing[0][locus] = HapAllelesPair[0][locus];
      HapAllelesPairNoMissing[1][locus] = HapAllelesPair[1][locus];
    } else { // assign alleles from MissingAlleles
      HapAllelesPairNoMissing[0][locus] = MissingAlleles[0][missingLocusOffset];
      HapAllelesPairNoMissing[1][locus] = MissingAlleles[1][missingLocusOffset];
    }
    missingLocusOffset += isMissing[locus]; 
  } 
}

void Haplotype::permuteMissingLoci(const vector<bool>& isMissing, const int numMissingLoci, const int permMissing, 
					const vector<int>& HapAlleles,  const vector<int>& MissingLoci, 
					vector<int>& HapAllelesNoMissing){
  vector<int> MissingAlleles = vector<int>(numMissingLoci);
  setMissingAlleles(numMissingLoci, permMissing, MissingAlleles, MissingLoci);

  int missingLocusOffset = 0; // used to loop over missing loci 
  //loop over loci to assign HapAllelesNoMissing
  for(int locus = 0; locus < NumberOfLoci; locus++ ) {
    if( ! isMissing[locus] ) { // assign alleles from HapAlleles
      HapAllelesNoMissing[locus] = HapAlleles[locus];
    } else { // assign alleles from MissingAlleles
      HapAllelesNoMissing[locus] = MissingAlleles[missingLocusOffset];
    }
    missingLocusOffset += isMissing[locus]; 
  } 
}

