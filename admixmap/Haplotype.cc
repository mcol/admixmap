#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>

using namespace std;

typedef struct  
{
   int haps[2];
} hapPair; 

// should initialize these when CompositeLocus object is initialized
const int numLoci = 3;
const int numAlleles[3] = {2, 3, 2};
int base[3]; // set with function setBaseForHapCode;


// arguments: integer, length of bit array
// returns: 1D array of bits representing integer
void ::intToBits(int n, const int length, bool *bits) 
{
  for( int i = length - 1; i >= 0; i-- ) {
    bits[i] = n % 2;
    // right shift n by 1 to get next bit 
    n >>= 1; 
  }
}

// updates 1D array of counting bases used to increment hap code
// should call once to initialize CompositeLocus object 
void ::setBaseForHapCode(int *base)
{
  base[numLoci - 1] = 1;
  for( int i = numLoci - 2; i >= 0; i-- ) {
    base[i] = base[i + 1] * numAlleles[i + 1];
  }
}

// updates 2D array of counting bases used to increment permMissing
void ::setBaseMissing(const int *missingLoci, const int numMissingLoci, int baseMissing[][2])
{
  baseMissing[numMissingLoci - 1][1] = 1;
  baseMissing[numMissingLoci - 1][0] = numAlleles[ missingLoci[numMissingLoci - 1] ];
  if( numMissingLoci > 1 ) {
    for( int i = numMissingLoci - 2; i >= 0; i-- ) {
      baseMissing[i][1] = baseMissing[i + 1][0] * numAlleles[ missingLoci[i + 1] ]; 
      baseMissing[i][0] = baseMissing[i][1] * numAlleles[ missingLoci[i] ]; 
    }
  }
}

// updates 2D array MissingAlleles with alleles specified by permMissing
void ::setMissingAlleles(const int baseMissing[][2], const int numMissingLoci, const int permMissing, int MissingAlleles[][2]) 
{
  int remainder = permMissing; 
  for( int i = 0; i < numMissingLoci; i++ ) { // loop over loci and gametes to extract alleles by integer division
    MissingAlleles[i][0] = 1 + remainder / baseMissing[i][0]; 
    remainder = remainder % baseMissing[i][0]; 
    MissingAlleles[i][1] = 1 + remainder / baseMissing[i][1]; 
    remainder = remainder % baseMissing[i][1]; 
  }
  // uncomment to display haplotype alleles  
  //   cout << "\npermMissing " << permMissing; 
  //   for( int col = 0; col < 2; col++ ) {
  //     cout << "\n";
  //     for( int row = 0; row < numMissingLoci; row++ ) {
  //       cout << MissingAlleles[row][col] << " ";
  //     }
  //   }
  //   cout << "\n";
}

// codes haplotypes starting at 0 by incrementing counter from right to left starting at 0 
// arguments: numAlleles, base, and hapAlleles are 1D arrays of length NumLoci
// returns: hap code as integer
int ::codeHapAllelesAsInt(const int *hapAlleles)
{
  int h = hapAlleles[numLoci - 1] - 1;
  for(int i = numLoci - 2; i >= 0; i-- ) {
    h += (hapAlleles[i] - 1) * base[i];
  }
  return( h );
}

// decode hap code to array of integers
// updates: 1D array hapAlleles
void ::decodeIntAsHapAlleles(const int h, int *hapAlleles)
{
  int remainder = h;
  // loop over loci to extract alleles by integer division
  for(int i = 0; i < numLoci; i++ ) {
    hapAlleles[i] = 1 + remainder / base[i];
    remainder = remainder % base[i];
  }
}


// arguments: HapAllelesPair is 2D array with NumLoci rows and 2 cols
// updates: 1D array of 2 integers coding haplotypes  
void ::codeHapAllelesPairAsIntPair(const int HapAllelesPair[][2], int *hpair)
{
  hpair[0] = HapAllelesPair[numLoci - 1][0] - 1;
  hpair[1] = HapAllelesPair[numLoci - 1][1] - 1;
  for( int i = numLoci - 2; i >= 0; i-- ) {
    hpair[0] += (HapAllelesPair[i][0] - 1) * base[i];
    hpair[1] += (HapAllelesPair[i][1] - 1) * base[i];
  }
}

// arguments: genotype as 2D array, hetLoci as array of col nums of het loci, isHet as array of length equal to numLoci
// updates: haplotype pair array with permHet th permutation of alleles at heterozygous loci  
void ::permuteHetLoci(const int *hetLoci, const bool *isHet, const int numHetLoci, const int permHet, 
				   const int Genotype[][2], int HapAllelesPair[][2])
{
  //recode permHet as array of bits, with length equal to NumHetLoci
  bool *permbits = new bool[numHetLoci];
  intToBits(permHet, numHetLoci, permbits);
  
  int hetLocusOffset = 0; // used to loop over heterozygous loci 
  //loop over loci to assign HapAllelesPair and swap alleles where element of permbits is 1;
  for(int locus = 0; locus < numLoci; locus++ ) {
    if( !isHet[locus] | !permbits[hetLocusOffset] ) { // assign alleles without swapping
      HapAllelesPair[locus][0] = Genotype[locus][0];
      HapAllelesPair[locus][1] = Genotype[locus][1];
    } else {
      if( permbits[hetLocusOffset] ) { // swap alleles
	HapAllelesPair[locus][0] = Genotype[locus][1];
	HapAllelesPair[locus][1] = Genotype[locus][0];
      }
    }
    hetLocusOffset += isHet[locus]; 
  } 
  delete [] permbits;
}
  
void ::permuteMissingLoci(const int *missingLoci, const bool *isMissing, const int numMissingLoci, const int permMissing, 
			const int HapAllelesPair[][2], const int baseMissing[][2], int HapAllelesPairNoMissing[][2]) 
{
  int (*MissingAlleles)[2] = new int[numMissingLoci][2];
  setMissingAlleles(baseMissing, numMissingLoci, permMissing, MissingAlleles);
  int missingLocusOffset = 0; // used to loop over missing loci 
  //loop over loci to assign HapAllelesPairNoMissing
  for(int locus = 0; locus < numLoci; locus++ ) {
    if( ! isMissing[locus] ) { // assign alleles from HapAllelesPair
      HapAllelesPairNoMissing[locus][0] = HapAllelesPair[locus][0];
      HapAllelesPairNoMissing[locus][1] = HapAllelesPair[locus][1];
    } else { // assign alleles from MissingAlleles
      HapAllelesPairNoMissing[locus][0] = MissingAlleles[missingLocusOffset][0];
      HapAllelesPairNoMissing[locus][1] = MissingAlleles[missingLocusOffset][1];
    }
    missingLocusOffset += isMissing[locus]; 
  } 
  delete [] MissingAlleles;
}

// arguments: genotype is 2D array of alleles with 2 rows and NumLoci cols
// updates: PossibleHapPairs, as stl vector of arrays of 2 integers
// call once for each individual at start of program 
void setPossibleHaplotypePairs(const int Genotype[][2], vector<hapPair> & PossibleHapPairs)
{
  int numHetLoci = 0;
  int numMissingLoci = 0;
  int numPerms;
  int numPermsHet = 1;
  int numPermsMissing = 1;
  bool *isHet = new bool[numLoci];
  bool *isMissing = new bool[numLoci];
  //loop over loci to determine which are heterozygous or missing
  for( int i = 0; i < numLoci; i++ ) {
    isHet[i] = false;
    isMissing[i] = false;
    if( (Genotype[i][0] == 0) | (Genotype[i][1] == 0)  ) { // missing genotype
      isMissing[i] = true;
      numMissingLoci ++;
      numPermsMissing *= numAlleles[i] * numAlleles[i];
    } else {
      if( Genotype[i][0] !=  Genotype[i][1] ) { // heterozygous
	isHet[i] = true; 
	numHetLoci++;
	numPermsHet *= 2;
      }
    }
  }
  //   cout << numMissingLoci << " missing genotypes with " << numPermsMissing << " permutations\n";
  //   cout << numHetLoci << " heterozygous genotypes with " << numPermsHet << " permutations\n";

  // loop over loci to calculate vectors hetLoci and missingLoci
  int *hetLoci = new int[numHetLoci];
  int *missingLoci = new int[numMissingLoci];
  int offsetHetLoci = 0;
  int offsetMissingLoci = 0;
  for( int i = 0; i < numLoci; i++ ) {
    if( isHet[i] ) {
      hetLoci[offsetHetLoci] = i;
      offsetHetLoci++;
    }
    if( isMissing[i] ) {
      missingLoci[offsetMissingLoci] = i;
      offsetMissingLoci++;
    }
  }
  
  int (*HapAllelesPair)[2] = new int[numLoci][2];
  int (*HapAllelesPairNoMissing)[2] = new int[numLoci][2];
  hapPair hpair;
  // loop over all possible ordered haplotype pairs compatible with unphased genotype
  // loop over permsHet to permute alleles at heterozygous loci
  for (int permHet = 0; permHet < numPermsHet; permHet++ ) {
    permuteHetLoci(hetLoci, isHet, numHetLoci, permHet, Genotype, HapAllelesPair);
    if( numMissingLoci == 0 ) {
      codeHapAllelesPairAsIntPair(HapAllelesPair, hpair.haps);
      PossibleHapPairs.push_back(hpair);
    } else {
      int (*baseMissing)[2] = new int[numMissingLoci][2];
      setBaseMissing(missingLoci, numMissingLoci, baseMissing);
      // loop over permsMissing to permute alleles at missing loci
      for ( int permMissing = 0; permMissing < numPermsMissing; permMissing++ ) {
	permuteMissingLoci(missingLoci, isMissing, numMissingLoci, permMissing, HapAllelesPair, baseMissing, 
			   HapAllelesPairNoMissing);
	codeHapAllelesPairAsIntPair(HapAllelesPairNoMissing, hpair.haps);
	PossibleHapPairs.push_back(hpair);
      }
      delete [] baseMissing;
    }
  }
  delete [] isHet;
  delete [] isMissing;
  delete [] hetLoci;
  delete [] missingLoci;
  delete [] HapAllelesPair;
  delete [] HapAllelesPairNoMissing;
}


int main() // tests these methods
{
  int hapAlleles[3] = {2, 3, 2}; 
  int Genotype[3][2] = { {1, 1}, {1, 2}, {0, 0} };
  vector<hapPair> PossibleHapPairs;
  
  setBaseForHapCode(base);
  
  int h = codeHapAllelesAsInt(hapAlleles);
  cout << "hap as integer " << h << "\n";
  
  int *hapDecoded = new int[numLoci];
  decodeIntAsHapAlleles(11, hapDecoded);
  cout << "hap as alleles ";
  for( int i = 0; i< numLoci; i++ ) {
    cout << hapDecoded[i] << " ";
  }
  cout << "\n";

  setPossibleHaplotypePairs(Genotype, PossibleHapPairs);
  cout << PossibleHapPairs.size() << " possible hap pairs\n";
  for( int perm = 0; perm < PossibleHapPairs.size(); perm++ ) {
    cout << PossibleHapPairs[perm].haps[0] << " " << PossibleHapPairs[perm].haps[1] << "\n";
  }
}  

