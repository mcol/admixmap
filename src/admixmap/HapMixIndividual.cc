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

const array_of_allelefreqs* HapMixIndividual::HaploidGenotypeProbs;
const double* HapMixIndividual::DiploidGenotypeProbs;


HapMixIndividual::HapMixIndividual(){

}

HapMixIndividual::HapMixIndividual(int number, const Options* const options, const InputData* const Data){

  GenotypesMissing = new bool*[numChromosomes];
  for( unsigned int j = 0; j < numChromosomes; j++ ){

  }  

  Individual::Initialise(options, Data);
  myNumber = number;

  unsigned long numhaploid = 0, numdiploid = 0, numhaploidX = 0, numdiploidX = 0;
  //unsigned numCompositeLoci = Loci->GetNumberOfCompositeLoci();
  //  for(unsigned j  = 0; j < numCompositeLoci; ++j){
  unsigned locus = 0;
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    GenotypesMissing[j] = new bool[ Loci->GetSizeOfChromosome(j) ];
    for(unsigned jj = 0; jj < Loci->GetSizeOfChromosome(j); ++jj){

      std::vector<unsigned short> g = Data->GetGenotype(locus, myNumber, options->getgenotypesSexColumn());
      
      if(g.size()==1){//haploid
        genotypes.push_back(g[0]);
        if(g[0]==0)GenotypesMissing[j][jj] = true;
        if(j== X_posn)++numhaploidX;
        else ++numhaploid;
        
      }
      else if(g.size()==2){//diploid
        switch(g[0]+g[1]){
        case 0:{//0,0
          genotypes.push_back(0);
          GenotypesMissing[j][jj] = true;
          break;
        }
        case 2:{//1,1
          genotypes.push_back(1);
          break;
        }
        case 3:{//1,2
          genotypes.push_back(3);
          break;
        }
        case 4:{//2,2
          genotypes.push_back(2);
          break;
        }
        default:{
          throw string("Invalid allele coding");
          break;
        }
          
          
        }
        if(j== X_posn)++numdiploidX;
        else ++numdiploid;
      }
      
    }
  }
  Data->CheckGenotypes(numhaploid, numdiploid, numhaploidX, numdiploidX, myNumber);
  isHaploid = (bool)(numdiploid+numdiploidX == 0); 
}

HapMixIndividual::~HapMixIndividual(){

}

void HapMixIndividual::SetStaticMembers(Genome* const pLoci, const Options* const options, const array_of_allelefreqs* const haploidGenotypeProbs, const double* const diploidGenotypeProbs){
  HaploidGenotypeProbs = haploidGenotypeProbs;
  DiploidGenotypeProbs = diploidGenotypeProbs;

  Individual::SetStaticMembers(pLoci, options);

}

///Indicates whether genotype is missing at all simple loci within a composite locus
bool HapMixIndividual::GenotypeIsMissing(unsigned int locus)const {
  return genotypes[locus]==0;
}
///Indicates whether genotype is missing at a simple locus
//used by HW score test
bool HapMixIndividual::simpleGenotypeIsMissing(unsigned locus)const{
  return genotypes[locus]==0;
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
    if(genotypes[locus0])//genotype not missing
      C->SetGenotypeProbs(DiploidGenotypeProbs + locus0*Populations*3 + genotypes[locus0]-1/*-1 to offset first genotype back to zero*/, GenotypesMissing[j]);
    else//missing genotype, it doesn't matter what probs we pass as they will not be used
      C->SetGenotypeProbs(0, GenotypesMissing[j]);
  }
  else{//haploid
    //pass allele freqs for this locus and this individual's observed allele at this locus
    C->SetGenotypeProbs((*HaploidGenotypeProbs)[locus0] + genotypes[locus0]-1 , GenotypesMissing[j]);
  }

  C->SetStateArrivalProbs(options->isRandomMatingModel(), diploid);
  logLikelihood.HMMisOK = false;//because forward probs in HMM have been changed
}
