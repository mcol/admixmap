/**
 *   HAPMIXMAP 
 *   HapMixIndividual.cc 
 *   Class to represent an individual in a hapmixmodel
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "HapMixIndividual.h"
#include "InputHapMixData.h"
#include "HapMixOptions.h"
#include <algorithm>
//#include "../base/DebugMacros.h"

const FreqArray* HapMixIndividual::HaploidGenotypeProbs;
const FreqArray* HapMixIndividual::DiploidGenotypeProbs;
GenotypeProbIterator HapMixIndividual::GPI;
HapMixGenome* HapMixIndividual::pG;

/** Possible haplotype pairs */
struct phps_t {
  int hpp_size; //< Size of possible happairs vector
  int h0; //< haplotype (haps[0])
  int h1; //< haplotype (haps[1])
};

/** Compare two sets of hpp_size, h0, h1.
 * 
 * Returns true if left (a) is less than right (b)
 */
struct phps_t_cmp {
  bool operator () (const phps_t& a, const phps_t& b) const {
    if      (a.hpp_size < b.hpp_size) { return true;  }
    else if (a.hpp_size > b.hpp_size) { return false; }
    else if       (a.h0 < b.h0)       { return true;  }
    else if       (a.h0 > b.h0)       { return false; }
    else if       (a.h1 < b.h1)       { return true;  }
    else                              { return false; }
  }
};

/** Unordered genotype probabilities */
struct ungp_t {
  int phpp0; //< probability of unordered genotype 0 (1, 1)
  int phpp1; //< probability of unordered genotype 1 (1, 2), (2, 1)
  int phpp2; //< probability of unordered genotype 2 (2, 2)
};

map<int, int> HapMixIndividual::ord2unord;

//Note: assuming SetStaticMembers is called first
HapMixIndividual::HapMixIndividual(int number, const HapMixOptions* const options, 
				   InputHapMixData* const Data, const double* GlobalTheta)
  :Individual(number){

  GenotypesMissing = new bool*[numChromosomes];
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    GenotypesMissing[j] = new bool[ Loci->GetSizeOfChromosome(j) ];
  }  

  //retrieve genotypes
  Data->GetGenotype(number, *Loci, &genotypes, GenotypesMissing);
  DeleteGenotypes();

  //isHaploid = (bool)(genotypes[0][0].size()==1);//note: assumes at least one autosome before X-chr
  isHaploid = Data->GetHapMixGenotype( number, *Loci, &hgenotypes, GenotypesMissing);
  Individual::Initialise(options, Data);

  Theta = const_cast<double*>(GlobalTheta);

  unsigned numCompositeLoci = Loci->GetNumberOfCompositeLoci();
  // loop over composite loci to set possible haplotype pairs compatible with genotype 
  for(unsigned j = 0; j < numCompositeLoci; ++j) {
    ploidy p = ( !isHaploid && (!Loci->isXLocus(j) || SexIsFemale)) ? diploid : haploid;
    SetPossibleHaplotypePairs(hgenotypes.begin()+j, PossibleHapPairs[j], p); 
    //HapMixGenotypeIterator G(hgenotypes.begin() + j, p);
    //(*Loci)(j)->HaplotypeSetter.setPossibleHaplotypePairs(&G, PossibleHapPairs[j]);
    
    // initialise sampledHapPairs with the first of the possible happairs. 
    // if only one possible happair or if annealing (which uses hamiltonian sampler), sampling of hap pair will be skipped.
    sampledHapPairs.push_back(PossibleHapPairs[j][0]);
  }
  //Now the PossibleHapPairs have ben determined and missing genotype indicators have been set, 
  //the genotypes are deleted as they are no longer needed 
  if( options->getHWTestIndicator())SetMissingGenotypes();
  
  /* Allocate space for unordered genotype probs
     They have form of vector of vectors of vectors of doubles.
     condition is : if genotypes have been masked (because then we are using UnorderedProbs to output posterior predictive Genotype Probs)
                    OR this is a test individual and the score test is on (only test individuals are included in the score test)
                    OR there are no test individuals and the score test is on (in this case, all individuals are included in the score test).
  */
  if(( options->OutputCGProbs() && Data->IsTestIndividual(number)) ||
     ( options->getTestForAllelicAssociation() && 
       (Data->IsTestIndividual(number) || !Data->getNumberOfTestIndividuals())
       )
     ){
    vector<double> v1 = vector<double>(1);
    vector<vector<double> > v3 = vector<vector<double> >(3, v1);
    UnorderedProbs = vector<vector<vector<double> > >(numCompositeLoci, v3);

    /*
     * Initialize values for observed genotypes. Information
     * is read from PossibleHapPairs and decoded into probabilities.
     * 
     * In other words, we're inferring unordered genotype probabilities
     * from PossibleHapPairs.
     */
#define SET_UNORD_PROBS(LOCUS_NO, V0, V1, V2) \
    UnorderedProbs[LOCUS_NO][0][0] = V0; \
    UnorderedProbs[LOCUS_NO][1][0] = V1; \
    UnorderedProbs[LOCUS_NO][2][0] = V2;
    
    /** Macro to define a mapping from number of haplotype pairs
     * and haplotype values to unordered genotype probabilities.
     */
#define SET_PHPS(HPP_SIZE, H0, H1, P0, P1, P2) \
    phps.hpp_size = HPP_SIZE; \
    phps.h0 = H0; \
    phps.h1 = H1; \
    ungp.phpp0 = P0; \
    ungp.phpp1 = P1; \
    ungp.phpp2 = P2; \
    phps_map[phps] = ungp

    map<phps_t, ungp_t, phps_t_cmp> phps_map;
    ungp_t ungp;
    phps_t phps;

    /** Mapping from number of haplotype pairs and haplotype values
     * to unordered genotype probabilities.
     */
    SET_PHPS(1, 0, -1, /* => */ 1, 0, 0); // haploid
    SET_PHPS(1, 1, -1, /* => */ 0, 0, 1); // haploid
    SET_PHPS(1, 0,  0, /* => */ 1, 0, 0);
    SET_PHPS(1, 1,  1, /* => */ 0, 0, 1);
    SET_PHPS(2, 0,  1, /* => */ 0, 1, 0);
    SET_PHPS(2, 1,  0, /* => */ 0, 1, 0);

    map<phps_t, ungp_t, phps_t_cmp>::const_iterator phpp_i;

    for (unsigned j = 0; j < numCompositeLoci; ++j) {
      if (!GenotypeIsMissing(j)) {
        // Genotype is present, set the unordered probability
        phps.hpp_size = PossibleHapPairs[j].size();
        phps.h0 = PossibleHapPairs[j][0].haps[0];
        phps.h1 = PossibleHapPairs[j][0].haps[1];
        phpp_i = phps_map.find(phps);
        if (phps_map.end() == phpp_i) {
          throw string("HapMixIndividual(): Something's wrong with PossibleHapPairs.");
        }
        SET_UNORD_PROBS(j,  phpp_i->second.phpp0,
                            phpp_i->second.phpp1,
                            phpp_i->second.phpp2);
      }
    }
    
    /*
     * Map from ordered to unordered genotype probabilities
     */
    if (ord2unord.size() == 0) {
      ord2unord[0] = 0;
      ord2unord[1] = 1;
      ord2unord[2] = 1;
      ord2unord[3] = 2;
    }
  }
}

///sets possible hap pairs for a single SNP
//TODO?? extend to compound loci
void HapMixIndividual::SetPossibleHaplotypePairs(const vector<unsigned short>::const_iterator GI, vector<hapPair> &PossibleHapPairs, ploidy p){
  //if(Genotype.size()!=1)throw string("Invalid call to Individual::SetPossibleHapPairs()");
  //  hapPair hpair;
  PossibleHapPairs.clear();
  if(p == haploid){
    //only one possibility - the observed haplotype (-1 denotes haplotype not happair)
    PossibleHapPairs.push_back(hapPair(*GI - 1 ,-1));
  }
  else{//diploid
    switch(*GI)
      {
      case(0):{//missing genotype -> 4 possibilities
        PossibleHapPairs.push_back(hapPair(0,0));
        PossibleHapPairs.push_back(hapPair(0,1));
        PossibleHapPairs.push_back(hapPair(1,0));
        PossibleHapPairs.push_back(hapPair(1,1));
        break;
      }
      case(1):{//1,1 -> only 1 possibility
        PossibleHapPairs.push_back(hapPair(0,0));
        break;
      }
      case(2):{//2,2 -> only 1 possibility
        PossibleHapPairs.push_back(hapPair(1,1));
        break;
      }
      case(3):{//1,2 -> 2 possibilities
        PossibleHapPairs.push_back(hapPair(0,1));
        PossibleHapPairs.push_back(hapPair(1,0));
        break;
      }
      default:{
        throw string("ERROR: Invalid genotype passed to HapMixIndividual::SetPossibleHapPairs");
      }
      }
  }
}
void HapMixIndividual::SetMissingGenotypes(){
  //allocates and sets an array of bools indicating whether genotypes at each locus are missing
  //used in HW score test; NB call before genotypes are deleted
  if(hgenotypes.size()==0)throw string("determining missing genotypes after genotypes have been deleted");
  missingGenotypes = new bool[Loci->GetTotalNumberOfLoci()];
  unsigned index = 0;
  for(unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); ++j)
    //for(int k = 0; k < Loci->getNumberOfLoci(j); ++k){
      missingGenotypes[index++] = (hgenotypes[j] == 0);
  //}
}

HapMixIndividual::~HapMixIndividual(){

}

// void HapMixIndividual::SetStaticMembers( Genome* const pLoci, const Options* const options){
//   Individual::SetStaticMembers(pLoci, options);
// }

void HapMixIndividual::SetGenotypeProbs(HapMixGenome* const G, const FreqArray& haploidGenotypeProbs, const FreqArray& diploidGenotypeProbs){
  pG = G;
  HaploidGenotypeProbs = &haploidGenotypeProbs;
  DiploidGenotypeProbs = &diploidGenotypeProbs;
}

///Indicates whether genotype is missing at all simple loci within a composite locus
bool HapMixIndividual::GenotypeIsMissing(unsigned int locus)const {
  return hgenotypes[locus]==0;
}
///Indicates whether genotype is missing at a simple locus
//used by HW score test
bool HapMixIndividual::simpleGenotypeIsMissing(unsigned locus)const{
  return hgenotypes[locus]==0;
}

/***
    Updates inputs to HMM for chromosome j
*/
void HapMixIndividual::UpdateHMMInputs(unsigned int j, const Options& , 
				       const double* const , const vector<double> ) {

  HapMixChromosome* C = pG->getHapMixChromosome(j);
  const bool diploid = !isHaploid && (j!=X_posn || SexIsFemale);
  const unsigned locus0 = C->GetLocus(0);//index of first locus on chromosome

  if(diploid){
    //pass diploid genotype probs for this locus and this individual's observed genotype
    //genotypes array doubles as indicator for missing genotype (missing = 0)
    //if(genotypes[locus0]){//genotype not missing
    GPI.assign(DiploidGenotypeProbs, &hgenotypes, 3, locus0);
    //}else//missing genotype, it doesn't matter what probs we pass as they will not be used
      //C->SetGenotypeProbs(0, GenotypesMissing[j]);
  }
  else{//haploid
    //pass allele freqs and this individual's observed genotype
    GPI.assign( HaploidGenotypeProbs, &hgenotypes, 2, locus0);
  }

  C->HMM->SetGenotypeProbs( GPI, GenotypesMissing[j]);
  logLikelihood.HMMisOK = false;
}

/// calculate genotype probs as an average over conditional probs of hidden states.
void HapMixIndividual::calculateUnorderedGenotypeProbs(const Options& options){
  unsigned locus = 0; 
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    //update HMM if required
     if( !logLikelihood.HMMisOK ) {
       UpdateHMMInputs(j, options, Theta, _rho);
     } 
     for(unsigned jj = 0; jj < Loci->GetSizeOfChromosome(j); ++jj, ++locus){
      // if genotype is missing, calculate
      if (GenotypeIsMissing(locus)) {
	calculateUnorderedGenotypeProbs(locus);
      }
      /*
       * If it's non-missing, its UnorderedGenotypeProbs were
       * inferred in the constructor.
       */
    }
  }
  logLikelihood.HMMisOK = true;
}

/// Same as calculateUnorderedProbs(void), but for j^th locus only
void HapMixIndividual::calculateUnorderedGenotypeProbs(unsigned j){
  // The following checks are turned off because of speed concerns
//  if (isHaploidIndividual()) {
//    string s = "HapMixIndividual::calculateUnorderedGenotypeProbs(int j) not implemented for haploid individuals";
//    throw(s);
//  }
//  if (!GenotypeIsMissing(j)) {
//    throw string("HapMixIndividual::calculateUnorderedGenotypeProbs(unsigned): cannot call this function for a non-missing genotype.");
//  }
//  if (not (Loci->GetNumberOfStates(j) == 2)) {
//    throw string("Trying to calculate UnorderedProbs but Loci->GetNumberOfStates(j) != 2");
//    return;
//  }
  int anc[2];
  
  bclib::pvector<double> orderedStateProbs = getStateProbs( Loci->GetChrNumOfLocus(j),
							     Loci->getRelativeLocusNumber(j));
  
  // set UnorderedProbs[j][*][0] to 0;
  for (vector<vector<double> >::iterator gi = UnorderedProbs[j].begin();
       gi != UnorderedProbs[j].end();
       ++gi)
    {
      (*gi)[0] = 0;
    }
  
  int ospIdx;
  
  /* Possible optimization: if the probability of the state
   * (orderedStateProbs[ospIdx]) is close to zero, it might have
   * a very little effect on the results, this state could
   * be skipped. A threshold of 1e-7 is too high.
   */
  
  //  orderedStateProbs.snapToZero();
  
  for (anc[0] = 0; anc[0] < NumHiddenStates; ++anc[0]) {
    /*
     * Calculating first index of ordered state probabilities.
     * It will be incremented in the inner loop after each iteration.
     */
    ospIdx = anc[0] * NumHiddenStates;
    for (anc[1] = 0; anc[1] < NumHiddenStates; ++anc[1], ++ospIdx) {
      
      if (orderedStateProbs[ospIdx] == 0) continue;
      
      /*
       * Calling a simplified version of getConditionalHapPairProbs
       * which sets only 0th and 3rd element of orderedGenotypeProbs
       * vector. 
       * 
       * PossibleHapPairs[j] are 
       * assumed to be all four: 0,0; 1,0; 0,1; 1,1.
       */
      /*
       * multiply result by conditional probs of anc and accumulate
       * result in array genotype probs (size 3 x number of loci)
       */
      UnorderedProbs[j][0][0] +=
	(*Loci)(j)->getFirstConditionalHapPairProbs(anc) * orderedStateProbs[ospIdx];
      UnorderedProbs[j][2][0] +=
	(*Loci)(j)->getLastConditionalHapPairProbs(anc) * orderedStateProbs[ospIdx];
    }
  }
  /*
   * Calculate the probability of heterozygous state as complement
   * to one.
   * (up0,   up1,         up2)
   *  calc.  complement   calc.
   */
  UnorderedProbs[j][1][0] = 1.0 - UnorderedProbs[j][0][0] - UnorderedProbs[j][2][0];
}

/** Get probabilities of hidden states from HMM */

const bclib::pvector<double>& HapMixIndividual::getStateProbs(int chromosome, int locus)const{
  return pG->getHapMixChromosome(chromosome)->HMM->GetHiddenStateProbs(!isHaploid && (chromosome!=(int)X_posn || SexIsFemale), locus);
}

/**
 * Return unordered probs as a vector of vectors of doubles.
 * 
 */
const vector<vector<double> >& HapMixIndividual::getUnorderedProbs(unsigned int j)const{
  return UnorderedProbs[j];
}

void HapMixIndividual::SampleJumpIndicators(int* SumArrivalCounts){
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    HapMixChromosome* C = pG->getHapMixChromosome(j);
    // don't need to sample jump indicators if globalrho and no conjugate update of admixture this iteration
    //sample number of arrivals, update SumNumArrivals and SumLocusAncestry
    //if( !Loci->isXChromosome(j) )
    C->HMM->SampleJumpIndicators(LocusAncestry[j], gametes[j], SumArrivalCounts);
    //    else 
    /// C->SampleJumpIndicators(LocusAncestry[j], gametes[j], SumArrivalCounts_X);
  } //end chromosome loop
}
/**
   Accumulate counts of arrivals of each state. ConcordanceCounts is a L * 2K  array, where the first K elements
   in each row are counts of discordant loci and the remaining K are counts of concordant loci
*/
void HapMixIndividual::AccumulateConcordanceCounts(int* ConcordanceCounts)const{
  unsigned locus = 0;
  const unsigned  K = NumHiddenStates;
  // const unsigned KSq = NumHiddenStates * NumHiddenStates;
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    const Chromosome* C = Loci->getChromosome(j);
    ++locus;//skip first locus on each chromosome
    for(unsigned locus = 1; locus < C->GetSize(); ++locus){
      
      //first gamete
      if( LocusAncestry[j][locus-1] != LocusAncestry[j][locus]) //discordant loci
        ++ConcordanceCounts[ locus*2*K + LocusAncestry[j][locus] ];
      else//concordant loci
        ++ConcordanceCounts[ locus*2*K + K + LocusAncestry[j][locus] ];
      
      //second gamete
      if(!j==X_posn || SexIsFemale){
        if( LocusAncestry[j][C->GetSize() + locus-1] != LocusAncestry[j][C->GetSize() + locus])
          ++ConcordanceCounts[locus*2*K + LocusAncestry[j][C->GetSize()+locus]];
        else
          ++ConcordanceCounts[locus*2*K + K + LocusAncestry[j][C->GetSize()+locus]];
 	
      }
      ++locus;
    }
  } //end chromosome loop
}
