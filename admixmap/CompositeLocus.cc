/** 
 *   ADMIXMAP
 *   CompositeLocus.cc 
 *   Class to represent a composite locus
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

#include "CompositeLocus.h"
#include "functions.h"

using namespace std;

/**
 * No-argument constructor for a CompositeLocus object, which represents a
 * composite locus. By default, constructs a simple diallelic locus with fixed 
 * allele frequencies.
 */
CompositeLocus::CompositeLocus()
{
  NumberOfLoci = 1;
  NumberOfAlleles = new int[NumberOfLoci];
  NumberOfAlleles[0] = 2;
  NumberOfStates = 2;
  base = 0;
  
  Populations = 0;
  RandomAlleleFreqs = false;
  NumberOfMergedHaplotypes = 0;
  NumberOfAlleles = 0;
  AlleleProbs = 0;
  HapPairProbs = 0;
  HapPairProbsMAP = 0;

  MergeHaplotypes = 0;
  HapLabels = 0;
}

CompositeLocus::~CompositeLocus()
{
  delete [] Label;
  if(base)
    delete[] base;
  delete[] HapPairProbs;
  delete[] HapPairProbsMAP;
  delete[] NumberOfAlleles;
  delete[] MergeHaplotypes;
  delete[] HapLabels;
}

/**
 * Changes number of loci in this composite locus, then sets each
 * locus to be diallelic (as for a SNP).
 *
 * NewNumberOfLoci - the number of loci to be represented by this
 * object.
 */
void CompositeLocus::SetNumberOfLoci( int NewNumberOfLoci )
{
   NumberOfLoci = NewNumberOfLoci;
   NumberOfAlleles = new int[NumberOfLoci];
   NumberOfAlleles[0] = 2;
   NumberOfStates = (int)pow( 2.0, NumberOfLoci );
   base = new int[NumberOfLoci];
   setBaseForHapCode();
}

/**
 * Gets the number of loci in this composite locus.
 * 
 * Returns:
 * the number of loci
 */
int CompositeLocus::GetNumberOfLoci()
{
   return( NumberOfLoci );
}

/**
 * Returns the number of states which can exist within this composite 
 * locus. If this composite locus is composed of a single locus, the 
 * number of states will be equal to the number of alleles. If there 
 * are more than one locus in this composite, the number of states 
 * will be equal to the number of possible haplotypes.
 *
 * Returns:
 * the number of states
 */
int CompositeLocus::GetNumberOfStates()
{
   return( NumberOfStates );
}

void CompositeLocus::SetNumberOfStates( int newNumberOfStates )
{
   NumberOfStates= newNumberOfStates;
}

void CompositeLocus::SetNumberOfPopulations(int pops){
  Populations = pops;
}
void CompositeLocus::SetRandomAlleleFreqs(bool b){
  RandomAlleleFreqs = b;
}
/**
 * Sets the number of alleles at a given locus.
 * Exits with an error if the locus doesn't exist.
 *
 * alleles - the number of alleles that exist at a given locus
 * locus - the given locus
 */
void CompositeLocus::SetNumberOfAllelesOfLocus( int locus, int alleles )
{
   if( locus > NumberOfLoci - 1 ){
      cout << "Input to SetNumberOfAllelesOfLocus > NumberOfLoci\n";
      exit(1);
   }

   NumberOfStates /= NumberOfAlleles[ locus ];
   NumberOfAlleles[ locus ] = alleles;
   NumberOfStates *= alleles;
}

/**
 * Gets the number of alleles at a given locus.
 * Exits with an error if the locus does not exist.
 *
 * locus - the locus
 * 
 * returns:
 * the number of alleles at the given locus
 */
int CompositeLocus::GetNumberOfAllelesOfLocus( int locus )
{
   if( locus > NumberOfLoci - 1 ){
      cout << "Input to GetNumberOfAllelesOfLocus > NumberOfLoci\n";
      exit(0);}

   return( NumberOfAlleles[ locus ] );
}

/**
 * Extends the composite locus by adding one locus (containing a given
 * number of alleles) to the end of the composite locus.
 *
 * alleles - the number of alleles in the locus to be added
 */

void CompositeLocus::AddLocus( int alleles )
{ 
  int temp[NumberOfLoci];
  for(int i=0;i<NumberOfLoci;++i)temp[i] = NumberOfAlleles[i];
  NumberOfLoci++;

  delete[] NumberOfAlleles;
  NumberOfAlleles = new int[NumberOfLoci];
  NumberOfAlleles[NumberOfLoci-1] = alleles;
  for(int i=0; i< NumberOfLoci-1; ++i)NumberOfAlleles[i] = temp[i]; 

  NumberOfStates *= alleles;
}

void CompositeLocus::Initialise(Matrix_d &AFreqs){
  AlleleProbs = alloc2D_d(NumberOfStates, Populations);
  //set size of array of haplotype pair probs
  HapPairProbs = new double[NumberOfStates * NumberOfStates * Populations * Populations];
  HapPairProbsMAP = new double[NumberOfStates * NumberOfStates * Populations * Populations];
  SetAlleleProbs(AFreqs);
  SetHapPairProbs();
  SetNoMergeHaplotypes();
  //Initialise HapPairProbsMAP to values in HapPairProbs
  for(int h0 = 0; h0 < NumberOfStates * NumberOfStates * Populations * Populations; ++h0)
    HapPairProbsMAP[h0] = HapPairProbs[h0]; 
}

void CompositeLocus::SetAlleleProbs(Matrix_d &AFreqs){
  try{
  //check dimensions
    if(!AlleleProbs)throw("AlleleProbs not allocated");
    if(AFreqs.GetNumberOfRows() != NumberOfStates - 1)throw("Wrong Number Of rows");
    if(AFreqs.GetNumberOfCols() != Populations)throw("Wrong number of cols");
  }
  catch(char * str){
    cout<<"Error in CompositeLocus::SetAlleleProbs: "<<str<<endl;
  }

  // put ones in last row
  for (int k = 0; k < Populations; k++ ) {
    AlleleProbs[NumberOfStates - 1][k] = 1.0; 
  }
  for( int a = 0; a < NumberOfStates - 1; a++ ) {
    for(int k = 0; k < Populations; ++k){
      // set allele probs in all but last row
      AlleleProbs[a][k] = AFreqs(a,k); 
      // accumulate subtraction from 1 in last row 
      AlleleProbs[NumberOfStates-1][k] -= AFreqs(a,k);
    }
  }

}

void CompositeLocus::SetNumberOfLabels()
{
   Label = new string[ NumberOfLoci ];
}

/**
 * Sets the name of this composite locus (usually from the 
 * allelefreqs.txt file
 *
 * newlabel - the name of this composite locus
 */
void CompositeLocus::SetLabel( int index, string newlabel )
{
   Label[index] = newlabel;
}

/**
 * Gets the name of this composite locus
 *
 * returns:
 * the name of this composite locus
 */
string CompositeLocus::GetLabel(int index)
{
   return( Label[index] );
}

/*
 * Given the ancestry of mother and father, and the haplotype pairs compatible with the genotype, returns one of
 * these haplotype pairs.
 * 
 * HapPairs - a vector of possible haplotype pairs (coded) compatible with genotype
 *
 * ancestry - a two-element vector of parental ancestry (e.g. 1,0 
 *   might represent european paternal and african maternal).
 *
 */
//TODO: change ancestry to an int *, change hap to hapPair
void CompositeLocus::SampleHapPair(int hap[2], std::vector<hapPair > &HapPairs, Vector_i ancestry){
  double Probs[HapPairs.size()];//1way array of hap.pair probs

  for(unsigned k = 0; k < HapPairs.size(); ++k){
    Probs[k] = HapPairProbs[ HapPairs[k].haps[0] * NumberOfStates * Populations * Populations + 
			     HapPairs[k].haps[1] * Populations * Populations +
			     ancestry(0) * Populations  + ancestry(1)];
  }

  int h = SampleFromDiscrete(Probs, HapPairs.size());

  //hap should really be a HapPair object
  // so we can do
  //return HapPairs[h];
  //or
  //hap = HapPairs[h];
  hap[0] = HapPairs[h].haps[0];
  hap[1] = HapPairs[h].haps[1];
}

/**
 * Called every time the haplotype frequencies change. Sets values in the 
 * four-dimensional matrix of probabilities of pairs of haplotypes.
 */

void CompositeLocus::SetHapPairProbs(){
  for(int h0 = 0; h0 < NumberOfStates; ++h0){
    for(int h1 = 0; h1 < NumberOfStates; ++h1){
      for(int k0 = 0; k0 < Populations; ++k0){
	for(int k1 = 0; k1 < Populations; ++k1)
	  HapPairProbs[h0 * NumberOfStates * Populations * Populations +
		       h1 * Populations * Populations +
		       k0 * Populations + k1] = AlleleProbs[h0][k0] * AlleleProbs[h1][k1];
      }
    }
  }
}

/*
  returns array of marginal alleleprobs for each allele of each locus, for ancestry state k,  within the composite locus
  P must be of correct dimensions ie a ragged array with dimensions #loci and #alleles
*/
void CompositeLocus::getLocusAlleleProbs(double **P, int k){

 for(int j = 0; j < NumberOfLoci; ++j)
    for(int jj = 0 ; jj < NumberOfAlleles[j]; ++jj)
      P[j][jj] = 0.0;

  int *hA =  new int[NumberOfLoci];

  for(int h = 0; h < NumberOfStates; ++h){//loop over all haplotypes
    decodeIntAsHapAlleles(h, hA);
    //compute marginal probs by summing over relevant hap probs
    for(int j = 0; j < NumberOfLoci; ++j)
      P[j][hA[j]-1] += AlleleProbs[h][k]; 
  }
  delete[] hA;
}

/**
 * Given a list of possible haplotype pairs, returns sums of probabilities of these haplotypes
 * given each possible ordered pair of locus ancestry states 
 * 
 * Haplotypes - a list of possible haplotypes pair labels compatible with the observed genotypes
 *
 * chibindicator is only to facilitate the Chib algorithm in Individual; instructs CompositeLocus to use HapPairProbsMAP
 * instead of HapPairProbs when allelefreqs are not fixed.
 *
 * RandomAlleleFreqs - indicates whether the allele frequencies are random 
 *
 * Probs
 * a k x k array with rows and cols indexing paternal and maternal ancestry. For 
 *   example, for African and European populations:
 *
 *       | AFR | EUR |
 *   ----|-----|-----|
 *   AFR | 0.5 | 0.2 |
 *   ----|-----|-----|
 *   EUR | 0.2 | 0.1 |
 *   ----|-----|-----|
 *
 *   n.b. the sum of all probabilities might not equal 1.0 - but
 *   probabilities are in correct proportions.
 */
// void CompositeLocus::GetGenotypeProbs(double **Probs, std::vector<hapPair > &HapPairs, bool chibindicator){
//   for(int k0 = 0; k0 < Populations; ++k0)for(int k1 = 0; k1 < Populations; ++k1){
//     Probs[k0][k1] = 0.0;
//     for(unsigned int h = 0; h < HapPairs.size() ; ++h)
//       if( chibindicator && RandomAlleleFreqs )
// 	Probs[k0][k1] += HapPairProbsMAP[HapPairs[h].haps[0]][HapPairs[h].haps[1]][k0][k1];
//       else
// 	Probs[k0][k1] += HapPairProbs[HapPairs[h].haps[0]][HapPairs[h].haps[1]][k0][k1];
//   }
// }

// doesn't really calculate posterior mode
// just sets to current value of hap freqs.  ok for Chib algorithm if strong prior
//this either misnamed or misdefined
void CompositeLocus::setHaplotypeProbsMAP()
{
  //HaplotypeProbs = HaplotypeProbsMAP;
  for(int h0 = 0; h0 < NumberOfStates * NumberOfStates * Populations * Populations; ++h0)
    HapPairProbs[h0] = HapPairProbsMAP[h0]; 
}

// arguments: integer, length of bit array
// returns: 1D array of bits representing integer
void CompositeLocus::intToBits(int n, const int length, bool *bits) 
{
  for( int i = length - 1; i >= 0; i-- ) {
    bits[i] = n % 2;
    // right shift n by 1 to get next bit 
    n >>= 1; 
  }
}

// updates 1D array of counting bases used to increment hap code
// should call once to initialize CompositeLocus object 
void CompositeLocus::setBaseForHapCode()
{
  delete[] base;
  base = new int[NumberOfLoci];
  base[NumberOfLoci - 1] = 1;
  for( int i = NumberOfLoci - 2; i >= 0; i-- ) {
    base[i] = base[i + 1] * NumberOfAlleles[i + 1];
  }
}

// updates 2D array of counting bases used to increment permMissing
void CompositeLocus::setBaseMissing(const int *missingLoci, const int numMissingLoci, int baseMissing[][2])
{
  baseMissing[numMissingLoci - 1][1] = 1;
  baseMissing[numMissingLoci - 1][0] = NumberOfAlleles[ missingLoci[numMissingLoci - 1] ];
  if( numMissingLoci > 1 ) {
    for( int i = numMissingLoci - 2; i >= 0; i-- ) {
      baseMissing[i][1] = baseMissing[i + 1][0] * NumberOfAlleles[ missingLoci[i + 1] ]; 
      baseMissing[i][0] = baseMissing[i][1] * NumberOfAlleles[ missingLoci[i] ]; 
    }
  }
}

// updates 2D array MissingAlleles with alleles specified by permMissing
void CompositeLocus::setMissingAlleles(const int baseMissing[][2], const int numMissingLoci, const int permMissing, int MissingAlleles[][2]) 
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
// arguments: NumberOfAlleles, base, and hapAlleles are 1D arrays of length NumberOfLoci
// returns: hap code as integer
int CompositeLocus::codeHapAllelesAsInt(const int *hapAlleles)
{
  int h = hapAlleles[NumberOfLoci - 1] - 1;
  for(int i = NumberOfLoci - 2; i >= 0; i-- ) {
    h += (hapAlleles[i] - 1) * base[i];
  }
  return( h );
}

// decode hap code to array of integers
// updates: 1D array hapAlleles
void CompositeLocus::decodeIntAsHapAlleles(const int h, int *hapAlleles)
{
  int remainder = h;
  // loop over loci to extract alleles by integer division
  for(int i = 0; i < NumberOfLoci; i++ ) {
    hapAlleles[i] = 1 + remainder / base[i];
    remainder = remainder % base[i];
  }
}


// arguments: HapAllelesPair is 2D array with NumberOfLoci rows and 2 cols
// updates: 1D array of 2 integers coding haplotypes  
void CompositeLocus::codeHapAllelesPairAsIntPair(const int HapAllelesPair[][2], int *hpair)
{
  hpair[0] = HapAllelesPair[NumberOfLoci - 1][0]- 1;
  hpair[1] = HapAllelesPair[NumberOfLoci - 1][1] - 1;
  for( int i = NumberOfLoci - 2; i >= 0; i-- ) {
    hpair[0] += (HapAllelesPair[i][0] - 1) * base[i];
    hpair[1] += (HapAllelesPair[i][1] - 1) * base[i];
  }
}

// arguments: genotype as 2D array, hetLoci as array of col nums of het loci, isHet as array of length equal to NumberOfLoci
// updates: haplotype pair array with permHet th permutation of alleles at heterozygous loci  
void CompositeLocus::permuteHetLoci(const bool *isHet, const int numHetLoci, const int permHet, 
				   unsigned short **Genotype, int HapAllelesPair[][2])
{
  //recode permHet as array of bits, with length equal to NumHetLoci
  bool permbits[numHetLoci];
  intToBits(permHet, numHetLoci, permbits);
  
  int hetLocusOffset = 0; // used to loop over heterozygous loci 
  //loop over loci to assign HapAllelesPair and swap alleles where element of permbits is 1;
  for(int locus = 0; locus < NumberOfLoci; locus++ ) {
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
}
  
void CompositeLocus::permuteMissingLoci(const bool *isMissing, const int numMissingLoci, const int permMissing, 
			const int HapAllelesPair[][2], const int baseMissing[][2], int HapAllelesPairNoMissing[][2]) 
{
  int MissingAlleles[numMissingLoci][2];
  setMissingAlleles(baseMissing, numMissingLoci, permMissing, MissingAlleles);
  int missingLocusOffset = 0; // used to loop over missing loci 
  //loop over loci to assign HapAllelesPairNoMissing
  for(int locus = 0; locus < NumberOfLoci; locus++ ) {
    if( ! isMissing[locus] ) { // assign alleles from HapAllelesPair
      HapAllelesPairNoMissing[locus][0] = HapAllelesPair[locus][0];
      HapAllelesPairNoMissing[locus][1] = HapAllelesPair[locus][1];
    } else { // assign alleles from MissingAlleles
      HapAllelesPairNoMissing[locus][0] = MissingAlleles[missingLocusOffset][0];
      HapAllelesPairNoMissing[locus][1] = MissingAlleles[missingLocusOffset][1];
    }
    missingLocusOffset += isMissing[locus]; 
  } 
}

// arguments: genotype is 2D array of alleles with 2 rows and NumberOfLoci cols
// updates: PossibleHapPairs, an stl vector of arrays of 2 integers
// call once for each individual at start of program 
void CompositeLocus::setPossibleHaplotypePairs(unsigned short **Genotype, vector<hapPair> &PossibleHapPairs)
{
  setBaseForHapCode();
  int numHetLoci = 0;
  int numMissingLoci = 0;
  int numPermsHet = 1;
  int numPermsMissing = 1;
  bool isHet[NumberOfLoci];
  bool isMissing[NumberOfLoci];
  //loop over loci to determine which are heterozygous or missing
  for( int i = 0; i < NumberOfLoci; i++ ) {
    isHet[i] = false;
    isMissing[i] = false;
    if( (Genotype[i][0] == 0) | (Genotype[i][1] == 0)  ) { // missing genotype
      isMissing[i] = true;
      numMissingLoci ++;
      numPermsMissing *= NumberOfAlleles[i] * NumberOfAlleles[i];
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
  int hetLoci[numHetLoci];
  int missingLoci[numMissingLoci];
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
  
  int HapAllelesPair[NumberOfLoci][2];
  int HapAllelesPairNoMissing[NumberOfLoci][2];
  hapPair hpair;
  // loop over all possible ordered haplotype pairs compatible with unphased genotype
  // loop over permsHet to permute alleles at heterozygous loci
  for (int permHet = 0; permHet < numPermsHet; permHet++ ) {
    permuteHetLoci(isHet, numHetLoci, permHet, Genotype, HapAllelesPair);
    if( numMissingLoci == 0 ) {
      codeHapAllelesPairAsIntPair(HapAllelesPair, hpair.haps);
      PossibleHapPairs.push_back(hpair);
    } else {
      int baseMissing[numMissingLoci][2];
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

//--------------------------------------
//      SCORE TEST FUNCTIONS


/**
 * Sets behaviour to not merge any haplotypes.
 * (Usually, rare haplotypes can be merged).
 * This function is only used for monitoring.
 */
void CompositeLocus::SetNoMergeHaplotypes()
{
   MergeHaplotypes = new int[ NumberOfStates ];
   for( int i = 0; i < NumberOfStates; i++ )
      MergeHaplotypes[i] = i;
   NumberOfMergedHaplotypes = NumberOfStates;
}

/**
 * Decides which haplotypes to merge for score test, based on
 * frequencies <=0.01.  alpha is proportion of admix of each
 * population (used for weighting).  This function is only used for
 * monitoring.
 */

void CompositeLocus::SetDefaultMergeHaplotypes( double *alpha)
{
   int count = 0, count2 = 0;
   double p, sumalpha = 0.0, sumAlleleProbs[Populations];
   int temp[ NumberOfStates ], Merged[ NumberOfStates ];
   fill(Merged, Merged + NumberOfStates, 0);
   MergeHaplotypes = new int[ NumberOfStates ];
   for( int j = 0; j < Populations; j++ ){
     sumalpha +=alpha[j];
     sumAlleleProbs[j] = 0.0;
   }

   for( int i = 0; i < NumberOfStates - 1; i++ ){
      p = 0;
      for( int j = 0; j < Populations; j++ ){
	p += alpha[j] * AlleleProbs[i][j];
	sumAlleleProbs[j] += AlleleProbs[i][j];
      }
      p /= sumalpha;
      if( p > 0.01 ){
         temp[ i ] = count;
         count++;
      }
      else
         temp[i] = NumberOfStates;
   }
   p = 0;
   for( int j = 0; j < Populations; j++ )
      p += alpha[j] * ( 1 - sumAlleleProbs[j] );
   p /= sumalpha;
   if( p > 0.01 ){
      temp[ NumberOfStates - 1 ] = count;
      count++;
   }
   else
      temp[ NumberOfStates - 1 ] = NumberOfStates;

   for( int i = 0; i < NumberOfStates; i++ ){
      if( temp[i] == NumberOfStates )
         MergeHaplotypes[i] = count;
      else{
         MergeHaplotypes[i] = temp[i];
         Merged[ count2 ] = i;
         count2++;
      }
   }

   if( count == NumberOfStates )
      NumberOfMergedHaplotypes = NumberOfStates;
   else
      NumberOfMergedHaplotypes = count + 1;

   HapLabels = new int[ NumberOfMergedHaplotypes * NumberOfLoci ];
   int *haplotype =  new int[NumberOfLoci];
   for( int i = 0; i < NumberOfMergedHaplotypes - 1; i++ ){
      decodeIntAsHapAlleles(Merged[i], haplotype);
      for(int j = 0; j < NumberOfLoci; ++j)HapLabels[i*NumberOfLoci +j] = haplotype[j]-1;
   }
   delete[] haplotype;
   for( int j = 0; j < NumberOfLoci; j++ )
     HapLabels[ (NumberOfMergedHaplotypes - 1)*NumberOfLoci + j ] = 99;
}

const int *CompositeLocus::GetHapLabels( int i )const
{
   return( HapLabels+i*NumberOfLoci );
}

/**
 * Returns number of haplotypes that have been merged.
 * This function is only used for monitoring.
 */
int CompositeLocus::GetNumberOfMergedHaplotypes()
{
   return( NumberOfMergedHaplotypes );
}

/**
 * Given haplotype, returns merged haplotype.
 * This function is only used for monitoring.
 */
int CompositeLocus::GetMergedHaplotype( int i )
{
   return( MergeHaplotypes[i] );
}

// apparently calculates contribution of allele freqs to marginal likelihood of model
// by subtracting log prior density from log posterior
// in current version, this function is not called anywhere
double GetMarginalLikelihood( Vector_d PriorAlleleFreqs, Vector_d AlleleCounts )
{
   double f = gsl_sf_lngamma( PriorAlleleFreqs.Sum() ) -
      gsl_sf_lngamma( PriorAlleleFreqs.Sum() + AlleleCounts.Sum() );
   for( int i = 0; i < PriorAlleleFreqs.GetNumberOfElements(); i++ )
      f += gsl_sf_lngamma( PriorAlleleFreqs(i) + AlleleCounts(i) )
         - gsl_sf_lngamma( PriorAlleleFreqs(i) );
   return(f);
}












