/** 
 *   ADMIXMAP
 *   CompositeLocus.cc 
 *   Class to represent a composite locus
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "CompositeLocus.h"
#include "rand.h"
#include <math.h>
#include <stdlib.h>
#include "functions.h"
#include <numeric>

using namespace std;

bool CompositeLocus::RandomAlleleFreqs;
int CompositeLocus::Populations;

/// Constructor 
CompositeLocus::CompositeLocus()
{
  NumberOfLoci = 0;
  NumberOfStates = 1;
  base = 0;
  
  Populations = 0;
  RandomAlleleFreqs = false;
  NumberOfMergedHaplotypes = 0;
  AlleleProbs = 0;
  SumAlleleProbs = 0;
#ifndef PARALLEL
  HapPairProbs = 0;
  HapPairProbsMAP = 0;
#endif

  MergeHaplotypes = 0;
  HapLabels = 0;
}

/// Destructor 
CompositeLocus::~CompositeLocus()
{
  delete[] base;
#ifndef PARALLEL
  if(HapPairProbsMAP != HapPairProbs)  delete[] HapPairProbsMAP;
  delete[] HapPairProbs;
#endif
  delete[] MergeHaplotypes;
  delete[] HapLabels;
  free_matrix( SumAlleleProbs, Populations);
}

// ******** Initialisation *************************
/**
 * sets number of simple loci in this composite locus, then sets each
 * locus to be diallelic (as for a SNP).
 *
 * NewNumberOfLoci - the number of simple loci to be represented by this
 * object.
 */
// void CompositeLocus::SetNumberOfLoci( int NewNumberOfLoci )
// {
//    NumberOfLoci = NewNumberOfLoci;
//    NumberOfAlleles = new int[NumberOfLoci];
//    NumberOfAlleles[0] = 2;
//    NumberOfStates = (int)pow( 2.0, NumberOfLoci );
//    base = new int[NumberOfLoci];
//    setBaseForHapCode();
// }

/**
 * sets number of alleles/haplotypes in this composite locus.
 *
 */

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
// void CompositeLocus::SetNumberOfAllelesOfLocus( int locus, int alleles )
// {
//    if( locus > NumberOfLoci - 1 ){
//       cout << "Input to SetNumberOfAllelesOfLocus > NumberOfLoci\n";
//       exit(1);
//    }

//    NumberOfStates /= NumberOfAlleles[ locus ];
//    NumberOfAlleles[ locus ] = alleles;
//    NumberOfStates *= alleles;
// }

/**
 * Extends the composite locus by adding one simple locus with given
 * number of alleles to the end of the composite locus.
 *
 * alleles - the number of alleles in the simple locus to be added
 */

void CompositeLocus::AddLocus( int alleles, string label = "")
{ 
  NumberOfAlleles.push_back(alleles);
  ++NumberOfLoci;
  NumberOfStates *= alleles;
  Label.push_back(label);
}

/**
 * Sets the label of a locus
 *
 * newlabel - the name of this composite locus
 */
// void CompositeLocus::SetLabel( int locus, string newlabel )
// {
//   if(locus < Label.size())
//     Label[locus] = newlabel;
// }

// ****** Accessors ***************************************
/**
 * Gets the number of simple loci in this composite locus.
 * 
 * Returns:
 * the number of simple loci
 */
int CompositeLocus::GetNumberOfLoci()const
{
   return( NumberOfLoci );
}

/**
 * Returns the number of states at this composite 
 * locus. If this composite locus comprises a single locus, the 
 * number of states will be equal to the number of alleles. If there 
 * are more than one locus in this composite, the number of states 
 * will be equal to the number of possible haplotypes.
 */
int CompositeLocus::GetNumberOfStates()const 
{
   return( NumberOfStates );
}


/**
 * Gets the number of alleles at a given simple locus.
 * Exits with an error if the simple locus does not exist.
 *
 */
int CompositeLocus::GetNumberOfAllelesOfLocus( int locus )const
{
   if( locus > NumberOfLoci - 1 ){
      cout << "Input to GetNumberOfAllelesOfLocus > NumberOfLoci\n";
      exit(0);}

   return( NumberOfAlleles[ locus ] );
}

/**
 * Gets the name of this composite locus
 */
const string CompositeLocus::GetLabel(int index)const
{
  if(index < NumberOfLoci)
    return( Label[index] );
  else return Label[NumberOfLoci-1];
}

/**
  returns array of marginal probs for each allele at each simple locus, for ancestry state k,  within the composite locus.
  P must be of correct dimensions: ie a ragged array with dimensions numloci and numalleles
*/
void CompositeLocus::getLocusAlleleProbs(double **P, int k)const{

 for(int j = 0; j < NumberOfLoci; ++j)
    for(int jj = 0 ; jj < NumberOfAlleles[j]; ++jj)
      P[j][jj] = 0.0;

  int *hA =  new int[NumberOfLoci];

  for(int h = 0; h < NumberOfStates; ++h){//loop over all haplotypes
    decodeIntAsHapAlleles(h, hA);
    //compute marginal probs by summing over relevant hap probs
    for(int j = 0; j < NumberOfLoci; ++j)
      P[j][hA[j]-1] += AlleleProbs[k*NumberOfStates + h]; 
  }
  delete[] hA;
}

void CompositeLocus::InitialiseHapPairProbs(const double* const AFreqs){
  AlleleProbs = AFreqs;//set AlleleProbs to point to allele freqs in AlleleFreqs
  AlleleProbsMAP = 0; //AFreqs;
  SumAlleleProbs = alloc2D_d(Populations, NumberOfStates);//allocates and fills with zeros
#ifndef PARALLEL
  //set size of array of haplotype pair probs
  HapPairProbs = new double[NumberOfStates * NumberOfStates * Populations * Populations];
  if(!RandomAlleleFreqs)AccumulateAlleleProbs();//if allelefreqs are fixed, SumAlleleProbs are initialised to AlleleProbs(==Allelefreqs)
  SetHapPairProbs();
#endif
}

void CompositeLocus::InitialiseHapPairProbsMAP(){
#ifndef PARALLEL
  HapPairProbsMAP = new double[NumberOfStates * NumberOfStates * Populations * Populations];
  //   //Initialise HapPairProbsMAP to values in HapPairProbs
  //   for(int h0 = 0; h0 < NumberOfStates * NumberOfStates * Populations * Populations; ++h0){
  //     HapPairProbsMAP[h0] = HapPairProbs[h0];
  SetHapPairProbsMAP();
#endif
}

/**
   Sets AlleleProbsMAP to point to FreqsMAP in AlleleFreqs . 
*/
void CompositeLocus::setAlleleProbsMAP(const double* const FreqsMAP) {
  AlleleProbsMAP = FreqsMAP;
}

/**
   Sets haplotype pair probabilities.
   Called every time the haplotype frequencies change. Sets elements in the  
   array of probabilities of ordered haplotype pairs for each ordered pair of ancestry states
*/
#ifndef PARALLEL
void CompositeLocus::SetHapPairProbs(){
  SetHapPairProbs(AlleleProbs, HapPairProbs);
}


/**
   SetsHapPairProbsMAP using AlleleProbsMAP
   Not done in parallel version as HapPairProbs and HapPairProbsMAP are not stored.
*/
void CompositeLocus::SetHapPairProbsMAP() {
#ifndef PARALLEL
  SetHapPairProbs( AlleleProbsMAP, HapPairProbsMAP);
#endif

}

/// private function
void CompositeLocus::SetHapPairProbs(const double* alleleProbs, double* hapPairProbs) {
  for(int h0 = 0; h0 < NumberOfStates; ++h0){
    for(int h1 = 0; h1 < NumberOfStates; ++h1){
      for(int k0 = 0; k0 < Populations; ++k0){
	for(int k1 = 0; k1 < Populations; ++k1)
	  hapPairProbs[h0 * NumberOfStates * Populations * Populations +
		       h1 * Populations * Populations +
		       k0 * Populations + k1] = alleleProbs[k0*NumberOfStates + h0] * alleleProbs[k1*NumberOfStates + h1];
      }
    }
  }
}
#endif

void CompositeLocus::SetHapPairProbsToPosteriorMeans(int iterations){
  if(RandomAlleleFreqs){//if fixed allele freqs there is nothing to do
#ifndef PARALLEL
    double *mu = new double[Populations * NumberOfStates];
#endif
    for(int k = 0; k < Populations; ++k){
      for(int h = 0; h < NumberOfStates; ++h)SumAlleleProbs[k][h] /= (double)iterations;
#ifdef PARALLEL
      //in parallel version, we don't store the HapPairProbs so instead of setting them to posterior means,
      //simply set AlleleProbs to their PM and they will bw used to compute genotype probs required for deviance at PM
      //this is ok as this is done at end, after main loop
      softmax(NumberOfStates, (double*)AlleleProbs+k*NumberOfStates, SumAlleleProbs[k]);
#else
      softmax(NumberOfStates, mu+k*NumberOfStates, SumAlleleProbs[k]);
#endif
      //if(RandomAlleleFreqs)for(int h = 0; h < NumberOfStates; ++h)SumAlleleProbs[k][h] *= (double)iterations;
    }

#ifndef PARALLEL
    SetHapPairProbs(mu, HapPairProbs);//sets HapPairProbs using posterior means of Haplotype Probs
    delete[] mu;
#endif
  }
}

void CompositeLocus::AccumulateAlleleProbs(){
  for (int k = 0; k < Populations; k++ ) {
    double* temp = new double[NumberOfStates];
    try{
      inv_softmax(NumberOfStates, AlleleProbs+k*NumberOfStates, temp);
    }
    catch(string str){
      string err = "Error accumulating alleleprobs: ";
      err.append(str);
      throw(err);
    }
    for(int h = 0; h < NumberOfStates; ++h)SumAlleleProbs[k][h] += temp[h];
    delete[] temp;
  }
}

/**
 * samples hap pair given locus ancstry
 * HapPairs - a vector of possible haplotype pairs (coded) compatible with genotype
 * ancestry - a two-element vector of parental ancestry (e.g. 1,0 
 *   might represent european paternal and african maternal).
 *
 */
#ifndef PARALLEL 
void CompositeLocus::SampleHapPair(hapPair* hap, const std::vector<hapPair > &HapPairs, const int ancestry[2])const{
  double* Probs = new double[HapPairs.size()];//1way array of hap.pair probs
  double* p = Probs;
  happairiter end = HapPairs.end();
  happairiter hiter = HapPairs.begin();
  for( ; hiter != end ; ++hiter, ++p) {
    *p = HapPairProbs[ hiter->haps[0] * NumberOfStates * Populations * Populations + 
		       hiter->haps[1] * Populations * Populations +
		       ancestry[0] * Populations  + ancestry[1]];
  }

  int h = Rand::SampleFromDiscrete(Probs, HapPairs.size());
  delete[] Probs;

  hap->haps[0] = HapPairs[h].haps[0];
  hap->haps[1] = HapPairs[h].haps[1];
}
#endif

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

/// updates 2D array of counting bases used to increment permMissing
void CompositeLocus::setBaseMissing(const vector<int> missingLoci, const int numMissingLoci, vector<int> baseMissing[2])
{
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
void CompositeLocus::setMissingAlleles(const vector<int> baseMissing[2], int numMissingLoci, int permMissing, 
				       vector<int> MissingAlleles[2]) 
{
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

void CompositeLocus::setMissingAlleles(int numMissingLoci, int permMissing, vector<int>& MissingAlleles, 
				       const vector<int>& MissingLoci) 
{
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
int CompositeLocus::codeHapAllelesAsInt(const int *hapAlleles)
{
  int h = hapAlleles[NumberOfLoci - 1] - 1;
  for(int i = NumberOfLoci - 2; i >= 0; i-- ) {
    h += (hapAlleles[i] - 1) * base[i];
  }
  return( h );
}

/// decode hap code to array of integers
/// updates: 1D array hapAlleles
void CompositeLocus::decodeIntAsHapAlleles(const int h, int *hapAlleles)const
{
  int remainder = h;
  // loop over loci to extract alleles by integer division
  for(int i = 0; i < NumberOfLoci; i++ ) {
    hapAlleles[i] = 1 + remainder / base[i];
    remainder = remainder % base[i];
  }
}

/// updates 1D array of 2 integers coding haplotypes  
/// arguments: HapAllelesPair is 2D array with NumberOfLoci rows and 2 cols
void CompositeLocus::codeHapAllelesPairAsIntPair(const vector<int> HapAllelesPair[2], int *hpair)
{
  hpair[0] = HapAllelesPair[0][NumberOfLoci - 1]- 1;
  hpair[1] = HapAllelesPair[1][NumberOfLoci - 1] - 1;
  for( int i = NumberOfLoci - 2; i >= 0; i-- ) {
    hpair[0] += (HapAllelesPair[0][i] - 1) * base[i];
    hpair[1] += (HapAllelesPair[1][i] - 1) * base[i];
  }
}

///returns a vector of length NumberOfLoci of numbers of copies of allele a at this locus, as encoded in haplotype pair h
const vector<int> CompositeLocus::getAlleleCounts(int a, const int* happair)const{
  vector<int> counts(NumberOfLoci);
  if(NumberOfStates == 2){
    counts[0] = (happair[0] == a-1) + (happair[1] == a-1);
  }
  else{
    int* hap1Alleles = new int[NumberOfLoci];
    int* hap2Alleles = new int[NumberOfLoci];
    decodeIntAsHapAlleles(happair[0], hap1Alleles);
    decodeIntAsHapAlleles(happair[1], hap2Alleles);
    
    fill(counts.begin(), counts.end(), 0);
    for(int l = 0; l < NumberOfLoci; ++l) counts[l] += (hap1Alleles[l]==a) + (hap2Alleles[l]==a);
    
    delete[] hap1Alleles;
    delete[] hap2Alleles;
  }
  return counts;    
}

///given a haplotype pair, returns a vector of haplotype counts
const vector<int> CompositeLocus::getHaplotypeCounts(const int* happair){
  vector<int> counts(NumberOfMergedHaplotypes);
  fill(counts.begin(), counts.end(), 0);

  counts[GetMergedHaplotype(happair[0]) ]++;
  counts[GetMergedHaplotype(happair[1]) ]++;

  return counts;
}

/// updates: haplotype pair array with permHet th permutation of alleles at heterozygous loci.  
/// arguments: genotype as 2D array, hetLoci as array of col nums of het loci, isHet as array of length equal to NumberOfLoci
void CompositeLocus::permuteHetLoci(const vector<bool> isHet, const int numHetLoci, const int permHet, 
				   const vector<vector<unsigned short> > Genotype, vector<int> HapAllelesPair[2])
{
  //recode permHet as array of bits, with length equal to NumHetLoci
  bool* permbits = new bool[numHetLoci];
  intToBits(permHet, numHetLoci, permbits);
  
  int hetLocusOffset = 0; // used to loop over heterozygous loci 
  //loop over loci to assign HapAllelesPair and swap alleles where element of permbits is 1;
  for(int locus = 0; locus < NumberOfLoci; locus++ ) {
    if( !isHet[locus] || !permbits[hetLocusOffset] ) { // assign alleles without swapping
      HapAllelesPair[0][locus] = Genotype[locus][0];
      HapAllelesPair[1][locus] = Genotype[locus][1];
    } else {
      if( permbits[hetLocusOffset] ) { // swap alleles
	HapAllelesPair[0][locus] = Genotype[locus][1];
	HapAllelesPair[1][locus] = Genotype[locus][0];
      }
    }
    hetLocusOffset += isHet[locus]; 
  } 
  delete[] permbits;
}
  
void CompositeLocus::permuteMissingLoci(const vector<bool> isMissing, const int numMissingLoci, const int permMissing, 
			const vector<int> HapAllelesPair[2], const vector<int> baseMissing[2], vector<int> HapAllelesPairNoMissing[2]) 
{
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

void CompositeLocus::permuteMissingLoci(const vector<bool>& isMissing, const int numMissingLoci, const int permMissing, 
					const vector<int>& HapAlleles,  const vector<int>& MissingLoci, 
					vector<int>& HapAllelesNoMissing) 
{
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

/// updates: PossibleHapPairs, an stl vector of arrays of 2 integers.
/// arguments: genotype is 2D array of alleles with 2 rows and NumberOfLoci cols. 
/// call once for each individual at start of program 
void CompositeLocus::setPossibleHaplotypePairs(const vector<vector<unsigned short> > Genotype, vector<hapPair> &PossibleHapPairs)
{
  setBaseForHapCode();
  if(Genotype[0].size()==1){//haploid genotypes
    setPossibleXHaplotypes(Genotype, PossibleHapPairs);
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
    if( (Genotype[i][0] == 0) || (Genotype[i][1] == 0)  ) { // missing genotype
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
    permuteHetLoci(isHet, numHetLoci, permHet, Genotype, HapAllelesPair);
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
void CompositeLocus::setPossibleXHaplotypes(const vector<vector<unsigned short> > Genotype, vector<hapPair> &PossibleHapPairs){
  //setBaseForHapCode();
  int numMissingLoci = 0;
  int numPermsMissing = 1;
  vector<bool> isMissing(NumberOfLoci);
  vector<int> HapAlleles = vector<int>(NumberOfLoci);
  vector<int> HapAllelesNoMissing = vector<int>(NumberOfLoci);

  for( int i = 0; i < NumberOfLoci; i++ ) {
    isMissing[i] = false;
    HapAlleles[i] = Genotype[i][0];
    if( (Genotype[i][0] == 0)   ) { // missing genotype
      isMissing[i] = true;
      numMissingLoci ++;
      numPermsMissing *= NumberOfAlleles[i];
    }
  }

  vector<int> missingLoci(numMissingLoci);
  int offsetMissingLoci = 0;
  for( int i = 0; i < NumberOfLoci; i++ ) {
    if( isMissing[i] ) {
      missingLoci[offsetMissingLoci] = i;
      offsetMissingLoci++;
    }
  }

  hapPair hpair;
  hpair.haps[1] = -1;//using -1 to denote *haplotype* rather than happair (must be nagative as alleles are counted from 0)

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
//--------------------------------------
//      SCORE TEST FUNCTIONS


// /**
//  * Sets behaviour to not merge any haplotypes.
//  * (Usually, rare haplotypes can be merged).
//  * This function is only used for monitoring.
//  */
// void CompositeLocus::SetNoMergeHaplotypes()
// {
//    MergeHaplotypes = new int[ NumberOfStates ];
//    for( int i = 0; i < NumberOfStates; i++ )
//       MergeHaplotypes[i] = i;
//    NumberOfMergedHaplotypes = NumberOfStates;
// }

/**
 * Decides which haplotypes to merge for score test, based on
 * frequencies <=0.01.  alpha is proportion of admix of each
 * population (used for weighting).  This function is only used for
 * monitoring.
 */
void CompositeLocus::SetDefaultMergeHaplotypes( const double* const alpha)
{
   int count = 0, count2 = 0;
   double p, sumalpha = 0.0;
   vector<int> temp( NumberOfStates ), Merged( NumberOfStates );
   fill(Merged.begin(), Merged.end(), 0);
   MergeHaplotypes = new int[ NumberOfStates ];
   for( int j = 0; j < Populations; j++ ){
     sumalpha +=alpha[j];
   }

   for( int i = 0; i < NumberOfStates - 1; i++ ){
      p = 0;
      for( int j = 0; j < Populations; j++ ){
	p += alpha[j] * AlleleProbs[j*NumberOfStates +i];
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
      p += alpha[j] * ( AlleleProbs[j*NumberOfStates +NumberOfStates-1] );
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
///returns index of ith merged haplotype, coded as integer
const int *CompositeLocus::GetHapLabels( int i )const
{
   return( HapLabels+i*NumberOfLoci );
}

/**
 * Returns number of haplotypes that have been merged.
 * This function is only used for monitoring.
 */
int CompositeLocus::GetNumberOfMergedHaplotypes()const
{
   return( NumberOfMergedHaplotypes );
}

/**
 * Given haplotype, returns merged haplotype.
 * This function is only used for monitoring.
 */
int CompositeLocus::GetMergedHaplotype( int i )const
{
   return( MergeHaplotypes[i] );
}













