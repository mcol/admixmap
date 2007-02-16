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

const FreqArray* HapMixIndividual::HaploidGenotypeProbs;
const FreqArray* HapMixIndividual::DiploidGenotypeProbs;
GenotypeProbIterator HapMixIndividual::GPI;


HapMixIndividual::HapMixIndividual(){

}

//Note: assuming SetStaticMembers is called first
HapMixIndividual::HapMixIndividual(int number, const Options* const options, const InputData* const Data){

  Individual::Initialise(number, options, Data);

//    GenotypesMissing = new bool*[numChromosomes];
//    for( unsigned int j = 0; j < numChromosomes; j++ ){
//      GenotypesMissing[j] = new bool[Loci->GetSizeOfChromosome(j)];
//    }  
   isHaploid = Data->GetHapMixGenotype( myNumber, options->getgenotypesSexColumn(), *Loci, &hgenotypes, GenotypesMissing);

//   Individual::Initialise(options, Data);

//   // loop over composite loci to set possible haplotype pairs compatible with genotype 
//   int numCompositeLoci = Loci->GetNumberOfCompositeLoci();
//   for(unsigned j = 0; j < (unsigned)numCompositeLoci; ++j) {
//     SetPossibleHaplotypePairs(j, genotypes.begin()+j, PossibleHapPairs[j]); 
    
//     // initialise sampledHapPairs with the first of the possible happairs. 
//     // if only one possible happair or if annealing (which uses hamiltonian sampler), sampling of hap pair will be skipped.
//     sampledHapPairs.push_back(PossibleHapPairs[j][0]);
//   }

}

HapMixIndividual::~HapMixIndividual(){

}

void HapMixIndividual::SetStaticMembers(Genome* const pLoci, const Options* const options,
                            const FreqArray& haploidGenotypeProbs, const FreqArray& diploidGenotypeProbs){
  HaploidGenotypeProbs = &haploidGenotypeProbs;
  DiploidGenotypeProbs = &diploidGenotypeProbs;

  Individual::SetStaticMembers(pLoci, options);

}

void HapMixIndividual::SetPossibleHaplotypePairs(unsigned locus, vector<unsigned short>::const_iterator g, vector<hapPair> &PossibleHapPairs ){
#ifdef PARALLEL
    //NOTE: X data not yet supported in parallel version


  //TODO

#else
    //  Note: assuming  only SNPs. Otherwise would require incrementing iterator by NumLoci in comp locus
  if(isHaploid || (locus==X_posn && !SexIsFemale) )//haploid genotype
      (*Loci)(locus)->setPossibleHaplotypes(g, PossibleHapPairs);

  else{//diploid genotype
    ;//TODO
    }
#endif
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
