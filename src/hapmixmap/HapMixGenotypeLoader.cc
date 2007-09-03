/* 
 *   HAPMIXMAP
 *   HapMixGenotypeLoader.cc 
 *   class to load and assign genotypes
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "HapMixGenotypeLoader.h"
#include "bclib/LogWriter.h"
#include "Genome.h"
#include "bclib/DataReader.h"
#include "bclib/StringConvertor.h"

HapMixGenotypeLoader::HapMixGenotypeLoader(){
  NumTestIndividuals = 0;
}

void HapMixGenotypeLoader::ReadTestGenotypes(const char* filename, bclib::LogWriter& Log){
  bclib::DataReader::ReadData(filename, testGeneticData_, Log, true);
  if(testGeneticData_.size()){
    NumTestIndividuals = testGeneticData_.size() - 1;
    FindTypedLoci();
  }
}

///fills a hapmix individual's genotype vector
bool HapMixGenotypeLoader::GetHapMixGenotype(int i, const Genome &Loci, 
					     vector<unsigned short>* genotypes, bool** Missing){
  const bool isTestIndiv = IsTestIndividual(i);

  unsigned long numhaploid = 0, numdiploid = 0, numhaploidX = 0, numdiploidX = 0;
  //unsigned numCompositeLoci = Loci.GetNumberOfCompositeLoci();
  //  for(unsigned j  = 0; j < numCompositeLoci; ++j){
  unsigned locus = 0;
  unsigned cclocus = 0;
  
  for( unsigned int j = 0; j < Loci.GetNumberOfChromosomes(); j++ ){
    bool isXChr = Loci.isXChromosome(j);
    for(unsigned jj = 0; jj < Loci.GetSizeOfChromosome(j); ++jj){
      
      std::vector<unsigned short> g = isTestIndiv ?
        GetTestGenotype(locus, &cclocus,i)//function will increment cclocus if this locus is typed
        : GenotypeLoader::GetGenotype(locus, i);

      //for backward-compatibility, allow diploid X-chr genotypes for males
      //    if(isXChr && !isFemale && g.size()==2)
      //  g.pop_back();
      
      
      if(g.size()==1){//haploid
        if(g[0] > 2)throwGenotypeError(i, locus,  Loci(locus)->GetLabel(0), g[0], 0, 2);//only SNPs allowed
        
        genotypes->push_back(g[0]);
        if(g[0]==0)Missing[j][jj] = true;
        else{//exclude missing genotypes from counts
          if(isXChr)++numhaploidX;
          else ++numhaploid;
        }
        
      }
      else if(g.size()==2){//diploid
        //check for bad genotype coding
        if( g[0]>2 || g[1]>2  || (g[0]==0 && g[1]!=0) || (g[0]!=0 && g[1]==0) )
          throwGenotypeError(i, locus, Loci(locus)->GetLabel(0), g[0], g[1], 2);
        
        unsigned short gg = g[0]+g[1];
        if(gg>0){//not missing
          if(isXChr)++numdiploidX;
          else ++numdiploid;
        }
        switch(gg){
          case 0:{//0,0
            genotypes->push_back(0);
            Missing[j][jj] = true;
            break;
          }
          case 2:{//1,1
            genotypes->push_back(1);
            break;
          }
          case 3:{//1,2
            genotypes->push_back(3);
            break;
          }
          case 4:{//2,2
            genotypes->push_back(2);
            break;
          }
          default:{
            throwGenotypeError(i, locus, Loci(locus)->GetLabel(0), g[0], g[1], 2);
            break;
          }
            
        }
      } //end if diploid
      else {//bad formatting eg 1,1,1
	cerr << "Unrecognized genotype format - "
	     << geneticData_[i][locus]
	     << " - at locus " << Loci(locus)->GetLabel(0) << " for individual " << geneticData_[i][0] << std::endl;
	exit(1);
      }
     ++locus;
    }
  }

  const string ID = isTestIndiv ? testGeneticData_[i - NumIndividuals][0] : geneticData_[i][0];

  CheckGenotypes((numhaploid + numdiploid + numhaploidX + numdiploidX), 
		 numhaploid, numdiploid, numhaploidX, numdiploidX, i, ID);
  bool isHaploid = (bool)(numdiploid+numdiploidX == 0); 

  return isHaploid;
}

unsigned HapMixGenotypeLoader::getNumberOfTestIndividuals()const{
  return NumTestIndividuals;
}

unsigned HapMixGenotypeLoader::getNumberOfIndividuals()const{
  return (NumIndividuals + NumTestIndividuals);
}

unsigned HapMixGenotypeLoader::getNumTypedLoci()const{
  return testGeneticData_[0].size()-1;//-1 to exclude header
}

bool HapMixGenotypeLoader::IsTestIndividual(unsigned i)const{
  return (bool)(i > getNumberOfIndividuals() - getNumberOfTestIndividuals() );
}

bool HapMixGenotypeLoader::isTypedLocus(unsigned locus)const{
  if(isTypedSNP.size())
    return isTypedSNP[locus];
  else
    return false;
}
void HapMixGenotypeLoader::clear(){
  GenotypeLoader::clear();
  for(unsigned i = 0; i < testGeneticData_.size(); ++i)
    testGeneticData_[i].clear();
  testGeneticData_.clear();
}

///search the loci in genotypesfile for loci in a testgenotypesfile
void HapMixGenotypeLoader::FindTypedLoci(){
  for(Vector_s::const_iterator j = geneticData_[0].begin()+1; j != geneticData_[0].end(); ++j){
    isTypedSNP.push_back(bclib::StringConvertor::isListedString(*j, testGeneticData_[0]));
  }
}

///gets a hapmix case-control genotype from the testgenotypes file
vector<unsigned short> 
HapMixGenotypeLoader::GetTestGenotype(unsigned locus, unsigned* obslocus,
				      int individual)const{
  vector<unsigned short> g;
  //if(individual <= NumIndividuals || !NumTestIndividuals) {//not case or control
    //throw string("Trying to read case-control genotype but there are none!");
  //}
  int col = 1 + SexColumn + *obslocus;
  if (IsPedFile)col = 1 + SexColumn + 2* (*obslocus);

  if(isTypedSNP[locus]){//is a typed locus
    g = GenotypeLoader::GetGenotype(testGeneticData_[individual-NumIndividuals][col]);
    ++(*obslocus);
  }
  else 
    //g has one zero. 
    //Doesn't matter whether haploid or diploid as a missing genotype 
    //is coded the same way and missing genotypes are ignored in most places anyway.
    g.push_back(0);

  return g;
}

void HapMixGenotypeLoader::GetGenotype(int i, const Genome &Loci,  
				       vector<genotype>* genotypes, bool** Missing)const{
  //these next lines should be removed (but cannot be yet as HapMixIndividual still calls this function)
  if(i > NumIndividuals && NumTestIndividuals) {
    GetTestGenotype(i, Loci, genotypes, Missing);
  }
  else
    GenotypeLoader::GetGenotype(i, Loci, genotypes, Missing);
}

///obsolete function for retrieving a hapmix case/control genotype, called by admixmap's GetGenotype function
//should be removed
void HapMixGenotypeLoader::GetTestGenotype(int i, const Genome &Loci, vector<genotype>* genotypes, bool** Missing)const{
  unsigned int simplelocus = 0;//simple locus counter
  unsigned complocus = 0;
  unsigned long numhaploid = 0;
  unsigned long numdiploid = 0;
  unsigned long numhaploidX = 0;
  unsigned long numdiploidX = 0;
  unsigned long cclocus = 0;//case-control locus counter
  //  unsigned numXloci = 0;

  //  const std::vector<std::string>& CCLoci = testGeneticData_[0];//labels of loci in case-control genotypesfile
  for(unsigned c = 0; c < Loci.GetNumberOfChromosomes(); ++c){
    bool isXchrm = Loci.isXChromosome(c);
    for(unsigned int j = 0; j < Loci.GetSizeOfChromosome(c); ++j){
      genotype G;
      // loop over composite loci to store genotype strings as pairs of integers in stl vector genotype
      const int numLoci = Loci.getNumberOfLoci(complocus);
      unsigned int count = 0;
      //  if(isXchrm)numXloci += numloci;
      for (int locus = 0; locus < numLoci; locus++) {
	const int numalleles = 2;
	vector<unsigned short> g;
	if(isTypedSNP[simplelocus]){
	  int col = 1 + SexColumn + cclocus;
	  if (IsPedFile)col = 1 + SexColumn + 2*cclocus;
	  
	  g = GenotypeLoader::GetGenotype(testGeneticData_[i-NumIndividuals][col]);
	  ++cclocus;
	  }
	else g.assign(2,0);//set genotypes at untyped loci 0
	if(g.size()==2)
	  if( (g[0] > numalleles) || (g[1] > numalleles))
	    throwGenotypeError(i, simplelocus, Loci(complocus)->GetLabel(0), 
			       g[0], g[1], numalleles );
	  else if (g.size()==1)
	    if( (g[0] > numalleles))
	      throwGenotypeError(i, simplelocus, Loci(complocus)->GetLabel(0), 
				 g[0], 0, numalleles );

	if(isXchrm){
	  if(g.size()==1)++numhaploidX;
	  else {//diploid X genotype
	    if(!isFemale(i)){//males cannot have diploid X genotypes
	      //cerr << "Genotype error in Individual " << i << ". Only females can have diploid X-chromosome genotypes.";
	      //exit(1);
	      //NOTE: allowing this for backward compatibility, for now
	      //instead remove second element
	      g.pop_back();
	    }
	    ++numdiploidX;
	  }
	}
	else{
	  if(g.size()==1)++numhaploid;
	  else ++numdiploid;
	}
	simplelocus++;
	G.push_back(g);
	count += g[0];
      }
      
      Missing[c][j] = (count == 0);
      
      genotypes->push_back(G);
      ++complocus;
    }
  }
  //check genotypes are valid
  if(numhaploidX + numdiploidX > 0){//some X genotypes present
    if(numdiploidX>0 && !isFemale(i)){//male with diploid X data
      cerr << "Genotype error in Individual " << i << ". Only females can have diploid X-chromosome genotypes.";
      exit(1);
    }
    if(numhaploidX>0 && isFemale(i)){
      if(numdiploidX >0){//female with haploid and diploid X data
	cerr << "Genotype error in Individual " << i << ". Females should have diploid X-chromosome genotypes.";
	exit(1);
      }
      if(numdiploid>0){//phased X data but unphased autosomal genotypes
	cerr << "Genotype error in Individual " << i << ". Female with diploid autosomes ands haploid X-Chromosome."; 
	exit(1);
      }
    }
  }
  else{//only autosomes
    if( numhaploid>0 && numdiploid > 0 ){//mixed haploid/diploid data
      cerr << "Genotype error in Individual " << i << ". Both haploid and diploid genotypes and no X chromosome.";
      exit(1);
    }
  }

}
