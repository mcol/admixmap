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
#include <algorithm>
//#include "../common/DebugMacros.h"

const FreqArray* HapMixIndividual::HaploidGenotypeProbs;
const FreqArray* HapMixIndividual::DiploidGenotypeProbs;
GenotypeProbIterator HapMixIndividual::GPI;

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
  int phpp0; //< probability of unordered genotype 0 (0, 0)
  int phpp1; //< probability of unordered genotype 1 (1, 0), (0, 1)
  int phpp2; //< probability of unordered genotype 2 (1, 1)
};

map<int, int> HapMixIndividual::ord2unord;

// No default constructor
//HapMixIndividual::HapMixIndividual(){
//
//}

//Note: assuming SetStaticMembers is called first
HapMixIndividual::HapMixIndividual(int number, const Options* const options, const InputData* const Data, const double* GlobalTheta){

  GenotypesMissing = new bool*[numChromosomes];
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    GenotypesMissing[j] = new bool[ Loci->GetSizeOfChromosome(j) ];
  }  

  //retrieve genotypes
  Data->GetGenotype(number, options->getgenotypesSexColumn(), *Loci, &genotypes, GenotypesMissing);
  DeleteGenotypes();

  //isHaploid = (bool)(genotypes[0][0].size()==1);//note: assumes at least one autosome before X-chr
  isHaploid = Data->GetHapMixGenotype( number, options->getgenotypesSexColumn(), *Loci, &hgenotypes, GenotypesMissing);
  Individual::Initialise(number, options, Data);

  Theta = const_cast<double*>(GlobalTheta);

  unsigned numCompositeLoci = Loci->GetNumberOfCompositeLoci();
  // loop over composite loci to set possible haplotype pairs compatible with genotype 
  for(unsigned j = 0; j < numCompositeLoci; ++j) {
    ploidy p = ( !isHaploid && (!Loci->isXLocus(j) || SexIsFemale)) ? diploid : haploid;
    //#ifdef PARALLEL
    //cannot use function in CompositeLocus because workers have no CompositeLocus objects
    SetPossibleHaplotypePairs(hgenotypes.begin()+j, PossibleHapPairs[j], p); 
    //#else
    //HapMixGenotypeIterator G(hgenotypes.begin() + j, p);
    //(*Loci)(j)->HaplotypeSetter.setPossibleHaplotypePairs(&G, PossibleHapPairs[j]);
    //#endif
    
    // initialise sampledHapPairs with the first of the possible happairs. 
    // if only one possible happair or if annealing (which uses hamiltonian sampler), sampling of hap pair will be skipped.
    sampledHapPairs.push_back(PossibleHapPairs[j][0]);
  }
  //Now the PossibleHapPairs have ben determined and missing genotype indicators have been set, 
  //the genotypes are deleted as they are no longer needed 
  if( options->getHWTestIndicator())SetMissingGenotypes();
  
  // Allocate space for unordered genotype probs
  // They have form of vector of vectors of vectors of doubles.
  // FIXME: allocation condition
  // UnorderedProbs should be allocated also for individuals who are
  // case/control. How to check if an individual is case/control?
  if(true or (options->getHapMixModelIndicator() && options->getTestForAllelicAssociation())){
    vector<double> v1 = vector<double>(1);
    vector<vector<double> > v3 = vector<vector<double> >(3, v1);
    UnorderedProbs = vector<vector<vector<double> > >(numCompositeLoci, v3);

    /*
     * Initialize values for observed genotypes. Information
     * is read from PossibleHapPairs and decoded into probabilities.
     * 
     * In other words, we're infering unordered genotype probabilities
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

void HapMixIndividual::SetStaticMembers( Genome* const pLoci, const Options* const options, 
                                         const FreqArray& haploidGenotypeProbs, const FreqArray& diploidGenotypeProbs){
  HaploidGenotypeProbs = &haploidGenotypeProbs;
  DiploidGenotypeProbs = &diploidGenotypeProbs;

  Individual::SetStaticMembers(pLoci, options);
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
    also sets Diploid flag in Chromosome (last arg of SetStateArrivalProbs)
*/
void HapMixIndividual::UpdateHMMInputs(unsigned int j, const Options* const options, 
				 const double* const , const vector<double> ) {

  Chromosome* C = Loci->getChromosome(j);
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

  C->SetGenotypeProbs( GPI, GenotypesMissing[j]);
  C->SetStateArrivalProbs(options->isRandomMatingModel(), diploid);
  logLikelihood.HMMisOK = false;//because forward probs in HMM have been changed
}

/** 
    function to calculate genotype probs as an average
    over conditional probs of hidden states.
    
    call this function from IndividualCollection.cc just after SampleLocusAncestry
    call new function in Chromosome to get state probs
    unchanged state probs from HMM
 */

void HapMixIndividual::calculateUnorderedGenotypeProbs(void){
  unsigned int numberCompositeLoci = Loci->GetNumberOfCompositeLoci();
  for(unsigned int j = 0; j < numberCompositeLoci; ++j) {
    // if genotype is missing, calculate
    if (GenotypeIsMissing(j)) {
      calculateUnorderedGenotypeProbs(j);
    }
    /*
     * If it's non-missing, its UnorderedGenotypeProbs were
     * inferred in HapMixIndividual constructor.
     */
  }
  return;
}

/**
   Same as Individual::calculateUnorderedProbs(void),
   but for j^th locus only
 */
void HapMixIndividual::calculateUnorderedGenotypeProbs(unsigned j){
  // The following checks are turned off because of speed concerns
//  if (isHaploidIndividual()) {
//    string s = "Individual::calculateUnorderedGenotypeProbs(int j) not implemented for haploid individuals";
//    throw(s);
//  }
//  if (!GenotypeIsMissing(j)) {
//    throw string("HapMixIndividual::calculateUnorderedGenotypeProbs(unsigned): Thou shalt not call this function for a non-missing genotype.");
//  }
//  if (not (Loci->GetNumberOfStates(j) == 2)) {
//    throw string("Trying to calculate UnorderedProbs but Loci->GetNumberOfStates(j) != 2");
//    return;
//  }
  int anc[2];
  
  pvector<double> orderedStateProbs = getStateProbs(
      not this->isHaploidIndividual(),
      Loci->getChromosomeNumber(j),
      Loci->getRelativeLocusNumber(j));

  // set UnorderedProbs[j][*][0] to 0;
  vector<vector<double> >::iterator gi;
  for (gi = UnorderedProbs[j].begin(); gi != UnorderedProbs[j].end(); ++gi) {
    (*gi)[0] = 0;
  }

  int ospIdx;
  orderedGenotypeProbs.assign(4, 0);
  
  /* Possible optimization: if the probability of the state
   * (orderedStateProbs[ospIdx]) is close to zero, it might have
   * a very little effect on the results, so this state could
   * be skipped. Unfortunately, threshold of 1e-7 is
   * still too high.
   */
  
//  orderedStateProbs.snapToZero();

  for (anc[0] = 0; anc[0] < NumHiddenStates; ++anc[0]) {
    for (anc[1] = 0; anc[1] < NumHiddenStates; ++anc[1]) {
      ospIdx = anc[0] * NumHiddenStates + anc[1];
      
      if (orderedStateProbs[ospIdx] == 0) continue;
      
      /*
       * Calling a simplified version of getConditionalHapPairProbs
       * which sets only 0th and 3rd element of orderedGenotypeProbs
       * vector.
       */
      (*Loci)(j)->getFirstAndLastConditionalHapPairProbs(
          orderedGenotypeProbs,
          PossibleHapPairs[j],
          anc);

      /*
       * multiply result by conditional probs of anc and accumulate
       * result in array genotype probs (size 3 x number of loci)
       * `ogpi' stands for ordered genotype probabilities index
       * 
       * Increment ogpi by 3 to skip the middle two elements:
       * (p0,  p1,   p2,   p3)
       *  get  skip  skip  get
       */
      for (int ogpi = 0; ogpi < 4; ogpi += 3) {
        UnorderedProbs[j][ord2unord[ogpi]][0] +=
            orderedGenotypeProbs[ogpi] * orderedStateProbs[ospIdx];
      }
      /*
       * Calculate the probability of heterozygous state as complement
       * to one.
       * (up0,   up1,         up2)
       *  calc.  complement   calc.
       */
			UnorderedProbs[j][1][0] = 1.0
					- UnorderedProbs[j][0][0]
					- UnorderedProbs[j][2][0];
    }
  }
}

/** Get probabilities of hidden states from HMM */

const pvector<double>& HapMixIndividual::getStateProbs(const bool isDiploid,const int chromosome,const int locus)const{
  return Loci->getChromosome(chromosome)->getHiddenStateProbs(isDiploid, locus);
}

/**
 * Return unordered probs as a vector of vectors of doubles.
 * Function AncestryAssocTest::Update() wants them this way.
 */
vector<vector<double> >& HapMixIndividual::getUnorderedProbs(
  const unsigned int j)
{
  return UnorderedProbs[j];
}

void HapMixIndividual::SampleJumpIndicators(int* SumArrivalCounts){
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    Chromosome* C = Loci->getChromosome(j);
    // don't need to sample jump indicators if globalrho and no conjugate update of admixture this iteration
    //sample number of arrivals, update SumNumArrivals and SumLocusAncestry
    //if( !Loci->isXChromosome(j) )
      C->SampleJumpIndicators(LocusAncestry[j], gametes[j], SumArrivalCounts);
    //    else 
    /// C->SampleJumpIndicators(LocusAncestry[j], gametes[j], SumArrivalCounts_X);
  } //end chromosome loop
}
