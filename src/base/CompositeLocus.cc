/** 
 *   CompositeLocus.cc 
 *   Class to represent a composite locus
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "CompositeLocus.h"
#include "bclib/rand.h"
#include <cmath>
#include <cstdlib>
#include <numeric>
#include "bclib/misc.h"

using namespace std;

bool CompositeLocus::RandomAlleleFreqs;
int CompositeLocus::Populations;
int CompositeLocus::PopulationsSquared;
int CompositeLocus::PopulationsSquared_x_3;

CompositeLocus::CompositeLocus()
{
  NumberOfLoci = 0;
  NumberOfStates = 1;
  
  Populations = 0;
  PopulationsSquared = 0;
  RandomAlleleFreqs = false;
  NumberOfMergedHaplotypes = 0;
  AlleleProbs = 0;
  SumAlleleProbs = 0;
  HapPairProbs = 0;
  HapPairProbsMAP = 0;

  MergeHaplotypes = 0;
  HapLabels = 0;
}

CompositeLocus::~CompositeLocus()
{
  if(HapPairProbsMAP != HapPairProbs)  delete[] HapPairProbsMAP;
  delete[] HapPairProbs;
  delete[] MergeHaplotypes;
  delete[] HapLabels;
  if(RandomAlleleFreqs)
    delete[] SumAlleleProbs;
}

// ******** Initialisation *************************

/// sets number of alleles/haplotypes in this composite locus.
void CompositeLocus::SetNumberOfStates( int newNumberOfStates )
{
   NumberOfStates= newNumberOfStates;
}

void CompositeLocus::SetNumberOfPopulations(int pops){
  Populations = pops;
  /*
   * Caching values for efficiency in function
   * getFirstAndLastConditionalGenotypeProbs
   */
  PopulationsSquared = pops * pops;
  PopulationsSquared_x_3 = pops * pops * 3;
}
void CompositeLocus::SetRandomAlleleFreqs(bool b){
  RandomAlleleFreqs = b;
}
/** Extends the composite locus by adding one simple locus with given
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
  HaplotypeSetter.AddLocus(alleles);
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
/** Gets the number of simple loci in this composite locus.
 * 
 * Returns:
 * the number of simple loci
 */
int CompositeLocus::GetNumberOfLoci()const
{
   return( NumberOfLoci );
}

/** Returns the number of states at this composite 
 * locus.
 * 
 * If this composite locus comprises a single locus, the 
 * number of states will be equal to the number of alleles. If there 
 * are more than one locus in this composite, the number of states 
 * will be equal to the number of possible haplotypes.
 */
int CompositeLocus::GetNumberOfStates()const 
{
   return( NumberOfStates );
}


/** Gets the number of alleles at a given simple locus.
 * 
 * Exits with an error if the simple locus does not exist.
 */
int CompositeLocus::GetNumberOfAllelesOfLocus( int locus )const
{
   if( locus > NumberOfLoci - 1 ){
      cout << "Input to GetNumberOfAllelesOfLocus > NumberOfLoci\n";
      exit(0);}

   return( NumberOfAlleles[ locus ] );
}

///Gets the name of this composite locus
const string CompositeLocus::GetLabel(int index)const
{
  if(index < NumberOfLoci)
    return( Label[index] );
  else return Label[NumberOfLoci-1];
}

/** returns array of marginal probs for each allele at each simple locus
 * 
 * for ancestry state k,  within the composite locus.
 * P must be of correct dimensions: ie a ragged array with dimensions
 * numloci and numalleles
 */
void CompositeLocus::getLocusAlleleProbs(double **P, int k)const{

 for(int j = 0; j < NumberOfLoci; ++j)
    for(int jj = 0 ; jj < NumberOfAlleles[j]; ++jj)
      P[j][jj] = 0.0;

  int *hA =  new int[NumberOfLoci];

  for(int h = 0; h < NumberOfStates; ++h){//loop over all haplotypes
    HaplotypeSetter.decodeIntAsHapAlleles(h, hA);
    //compute marginal probs by summing over relevant hap probs
    for(int j = 0; j < NumberOfLoci; ++j)
      P[j][hA[j]-1] += AlleleProbs[k*NumberOfStates + h]; 
  }
  delete[] hA;
}

void CompositeLocus::InitialiseHapPairProbs(const double* const AFreqs, bool AllHaploid){
  AlleleProbs = AFreqs;//set AlleleProbs to point to allele freqs in AlleleFreqs
  AlleleProbsMAP = 0; //AFreqs;
  SumAlleleProbs = const_cast<double*>(AlleleProbs);
  if(RandomAlleleFreqs){
    SumAlleleProbs = new double[Populations*NumberOfStates];
    fill(SumAlleleProbs, SumAlleleProbs+Populations*NumberOfStates, 0.0);
  }

  if(!AllHaploid){//some diploid data, need HapPairProbs
    //set size of array of haplotype pair probs
    HapPairProbs = new double[NumberOfStates * NumberOfStates * Populations * Populations];
    SetHapPairProbs();
  }

}

void CompositeLocus::InitialiseHapPairProbsMAP(){
  HapPairProbsMAP = new double[NumberOfStates * NumberOfStates * Populations * Populations];
  //   //Initialise HapPairProbsMAP to values in HapPairProbs
  //   for(int h0 = 0; h0 < NumberOfStates * NumberOfStates * Populations * Populations; ++h0){
  //     HapPairProbsMAP[h0] = HapPairProbs[h0];
  SetHapPairProbsMAP();
}

///Sets AlleleProbsMAP to point to FreqsMAP in AlleleFreqs.
void CompositeLocus::setAlleleProbsMAP(const double* const FreqsMAP) {
  AlleleProbsMAP = FreqsMAP;
}

/** Sets haplotype pair probabilities.
 * 
 * Called every time the haplotype frequencies change. Sets elements
 * in the array of probabilities of ordered haplotype pairs for each
 * ordered pair of ancestry states
 */
void CompositeLocus::SetHapPairProbs(){
  SetHapPairProbs(AlleleProbs, HapPairProbs);
}


/**
   SetsHapPairProbsMAP using AlleleProbsMAP
*/
void CompositeLocus::SetHapPairProbsMAP() {
  SetHapPairProbs( AlleleProbsMAP, HapPairProbsMAP);

}

/// private function
void CompositeLocus::SetHapPairProbs(const double* alleleProbs, double* hapPairProbs) {
  for(int h0 = 0; h0 < NumberOfStates; ++h0){
    for(int h1 = 0; h1 < NumberOfStates; ++h1){
      for(int k0 = 0; k0 < Populations; ++k0){
        for(int k1 = 0; k1 < Populations; ++k1) {
          hapPairProbs[h0 * NumberOfStates * Populations * Populations +
		       h1 * Populations * Populations +
		       k0 * Populations + k1] = alleleProbs[k0*NumberOfStates + h0] * alleleProbs[k1*NumberOfStates + h1];
        }
      }
    }
  }
}


void CompositeLocus::SetHapPairProbsToPosteriorMeans(int iterations){
  if(RandomAlleleFreqs){//if fixed allele freqs there is nothing to do

    double *mu = new double[Populations * NumberOfStates];

    for(int k = 0; k < Populations; ++k){
      for(int h = 0; h < NumberOfStates; ++h)SumAlleleProbs[k*NumberOfStates +h] /= (double)iterations;
      bclib::softmax(NumberOfStates, mu+k*NumberOfStates, SumAlleleProbs+k*NumberOfStates);

      //if(RandomAlleleFreqs)for(int h = 0; h < NumberOfStates; ++h)SumAlleleProbs[k*NumberOfStates+h] *= (double)iterations;
    }

    SetHapPairProbs(mu, HapPairProbs);//sets HapPairProbs using posterior means of Haplotype Probs
    delete[] mu;
  }
}

void CompositeLocus::AccumulateAlleleProbs(){
  for (int k = 0; k < Populations; k++ ) {
    double* temp = new double[NumberOfStates];
    bool* b = new bool[NumberOfStates];
    for(int s = 0; s < NumberOfStates; ++s)b[s] = (bool)(AlleleProbs[k*NumberOfStates+s]>0.0);
    try{
      bclib::inv_softmax(NumberOfStates, AlleleProbs+k*NumberOfStates, temp, b);
    }
    catch(string str){
      string err = "Error accumulating alleleprobs: ";
      err.append(str);
      throw(err);
    }
    for(int h = 0; h < NumberOfStates; ++h)SumAlleleProbs[k*NumberOfStates+h] += temp[h];
    delete[] temp;
    delete[] b;
  }
}


/**
 * Returns probabilities of ordered hap pairs conditional on hidden states
 */
void CompositeLocus::getConditionalHapPairProbs(bclib::pvector<double>& Probs, const std::vector<hapPair > &PossibleHapPairs, const int ancestry[2])const{
  //Note: probs should have length equal to the total number of possible diploid states ie  NumberOfStates^2 .
  // (in haploid case, we get the probs directly from alleleprobs/allelefreqs
  
  // Following check turned off because of speed:
  // In a 5-testing iterations this function is being called
  // over 86 million times.
//  if((int)Probs.size() != NumberOfStates*NumberOfStates)throw string("Wrongly sized vector passed to CompositeLocus::getConditionalHapPairProbs");
  fill(Probs.begin(), Probs.end(), 0.0);//fill with zeros

  const happairiter end = PossibleHapPairs.end();
  happairiter hiter = PossibleHapPairs.begin();//hiter points to elements of PossibleHapPairs
  int PopSq_x_NumberOfStates = PopulationsSquared * NumberOfStates;
      
  for( ; hiter != end ; ++hiter) {
    if(hiter->haps[1] >= 0){//diploid (haps (as opposed to happairs) have 2nd element -1
      //retrieve required element from HapPairProbs
      Probs[hiter->haps[0]*NumberOfStates + hiter->haps[1]] = 
          HapPairProbs[ hiter->haps[0] * PopSq_x_NumberOfStates + 
                        hiter->haps[1] * PopulationsSquared +
                        ancestry[0] * Populations  + ancestry[1]];
    }
    else{//haploid
      throw string("CompositeLocus::getConditionalHapPairProbs: ERROR: attempting to compute posterior genotypes probs for haploid individual!");
//       const double prob = AlleleProbs[ancestry[0]*NumberOfStates + hiter->haps[0]];
//       Probs[hiter->haps[0]] = prob;
//       sum += prob;  
     }
  }

	/*
	 * There's a need to normalize this vector only when some happairs
	 * are `not possible'. When all happairs are possible, this vector
	 * is already normalized, with precision to about 1e-5 (from what
	 * could be seen at runtime).
	 */
  if ((int)PossibleHapPairs.size() < NumberOfStates * NumberOfStates) {
    Probs.normalize();
  }
}

/** Simplified version of getConditionalHapPairProbs.
 * 
 * Sets only two values of the Probs vector.
 * 
 * This function is one of the main performance bottlenecks
 * of the allelic association score test. It's being called
 * K^2 (e.g. 64) times for each individual at each locus.
 */
double CompositeLocus::getFirstConditionalHapPairProbs(const int ancestry[2]) const
{
  return HapPairProbs[ /* 0 * PopSq_x_NumberOfStates + */ 
		      /* 0 * PopulationsSquared + */
		      ancestry[0] * Populations  + ancestry[1]];
}

double CompositeLocus::getLastConditionalHapPairProbs( const int ancestry[2]) const
{
  return HapPairProbs[ PopulationsSquared_x_3 + 
		       ancestry[0] * Populations  + ancestry[1]];
}

/**
 * samples hap pair given hidden states
 * HapPairs - a vector of possible haplotype pairs (coded) compatible with genotype
 * ancestry - a two-element vector of parental ancestry (e.g. 1,0 
 *   might represent european paternal and african maternal).
 *
 */
void CompositeLocus::SampleHapPair(hapPair* hap, const std::vector<hapPair > &PossibleHapPairs, const int ancestry[2])const{
  double* Probs = new double[PossibleHapPairs.size()];//1way array of hap.pair probs
  // getConditionalHapPairProbs(Probs, PossibleHapPairs, ancestry);
  happairiter end = PossibleHapPairs.end();
  happairiter hiter = PossibleHapPairs.begin();
  for( unsigned j = 0; j < PossibleHapPairs.size() ; ++j) {
    Probs[j] = HapPairProbs[ PossibleHapPairs[j].haps[0] * NumberOfStates * Populations * Populations + 
			     PossibleHapPairs[j].haps[1] * Populations * Populations +
			     ancestry[0] * Populations  + ancestry[1]];
  }
  //no need to renormalize for SampleFromDiscrete
  int h = bclib::Rand::SampleFromDiscrete(Probs, PossibleHapPairs.size());
  delete[] Probs;

  hap->haps[0] = PossibleHapPairs[h].haps[0];
  hap->haps[1] = PossibleHapPairs[h].haps[1];
}


///returns a vector of length NumberOfLoci of numbers of copies of allele a at this locus, as encoded in haplotype pair h
const vector<int> CompositeLocus::getAlleleCounts(int a, const int* happair)const{
  vector<int> counts(NumberOfLoci);
  if(NumberOfStates == 2){
    counts[0] = (happair[0] == a-1) + (happair[1] == a-1);
  }
  else{
    fill(counts.begin(), counts.end(), 0);

    int* hapAlleles = new int[NumberOfLoci];
    HaplotypeSetter.decodeIntAsHapAlleles(happair[0], hapAlleles);
    for(int l = 0; l < NumberOfLoci; ++l) counts[l] += (hapAlleles[l]==a);

    if(happair[1] >-1){
      HaplotypeSetter.decodeIntAsHapAlleles(happair[1], hapAlleles);
      for(int l = 0; l < NumberOfLoci; ++l) counts[l] += (hapAlleles[l]==a);
    }
    
    delete[] hapAlleles;
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
      HaplotypeSetter.decodeIntAsHapAlleles(Merged[i], haplotype);
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

/** Returns number of haplotypes that have been merged.
 * This function is only used for monitoring.
 */
int CompositeLocus::GetNumberOfMergedHaplotypes()const
{
   return( NumberOfMergedHaplotypes );
}

/** Given haplotype, returns merged haplotype.
 * This function is only used for monitoring.
 */
int CompositeLocus::GetMergedHaplotype( int i )const
{
   return( MergeHaplotypes[i] );
}














