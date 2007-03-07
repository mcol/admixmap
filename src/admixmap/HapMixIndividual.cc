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

const IFreqArray* HapMixIndividual::HaploidGenotypeProbs;
const IFreqArray* HapMixIndividual::DiploidGenotypeProbs;
GenotypeProbIterator HapMixIndividual::GPI;

// No default constructor
//HapMixIndividual::HapMixIndividual(){
//
//}

//Note: assuming SetStaticMembers is called first
HapMixIndividual::HapMixIndividual(
    int number,
    const IOptions* const options,
    const IInputData* const Data)
//try
{

  GenotypesMissing = new bool*[numChromosomes];
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    GenotypesMissing[j] = new bool[ Loci->GetSizeOfChromosome(j) ];
  }  
//  missingGenotypes = 0;//allocated later, if needed
  //retrieve genotypes
  Data->GetGenotype(number, options->getgenotypesSexColumn(), *Loci, &genotypes, GenotypesMissing);
  DeleteGenotypes();
//   for( unsigned int j = 0; j < numChromosomes; j++ ){
//     delete[] GenotypesMissing[j];
//   }  
//   delete[] GenotypesMissing;
//   GenotypesMissing = new bool*[numChromosomes];
//   for( unsigned int j = 0; j < numChromosomes; j++ ){
//     GenotypesMissing[j] = new bool[ Loci->GetSizeOfChromosome(j) ];
//   }  

  //isHaploid = (bool)(genotypes[0][0].size()==1);//note: assumes at least one autosome before X-chr

  isHaploid = Data->GetHapMixGenotype( number, options->getgenotypesSexColumn(), *Loci, &hgenotypes, GenotypesMissing);
  Individual::Initialise(number, options, Data);


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
  //DeleteGenotypes();

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

void HapMixIndividual::SetStaticMembers(
    IGenome* const pLoci,
    const IOptions* const options,
    const IFreqArray& haploidGenotypeProbs,
    const IFreqArray& diploidGenotypeProbs)
{
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

  IChromosome* C = Loci->getChromosome(j);
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
