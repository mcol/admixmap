/*
 *   GenotypeLoader.cc 
 *   class to load and assign genotypes
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "GenotypeLoader.h"
#include "bcppcl/LogWriter.h"
#include "Genome.h"
#include "bcppcl/DataReader.h"
#include "bcppcl/StringConvertor.h"

GenotypeLoader::GenotypeLoader(){
  NumIndividuals = 0;
  numDiploid = 0;
}

void GenotypeLoader::Read(const char* filename, LogWriter& Log){
  DataReader::ReadData(filename, geneticData_, Log);
  IsPedFile = determineIfPedFile();
  NumIndividuals = geneticData_.size() - 1;
}

unsigned GenotypeLoader::getNumberOfIndividuals()const{
  return NumIndividuals;
}

unsigned GenotypeLoader::NumLoci()const{
  return geneticData_[0].size();
}

/**
   extracts genotype for given individual at given locus
   \param locus the locus index
   \param individual individual index (row number)
   \param SexColumn index of sex column, or 0 if none
   \return genotype as vector of unsigned short integers
*/
vector<unsigned short> GenotypeLoader::GetGenotype(unsigned locus, int individual, int SexColumn)const{
  int col = 1 + SexColumn + locus;
  if (IsPedFile)col = 1 + SexColumn + 2*locus;
  return GetGenotype(geneticData_[individual][col]);
}


///convert a genotype string to a vector of unsigned ints
vector<unsigned short> GenotypeLoader::GetGenotype(const string genostring)const{
  vector<unsigned short> g;
 
  //strip quotes from string
  const std::string str = StringConvertor::dequote(genostring);
  if(str.length()==0){
    //if empty string, interpret as missing genotype for backward compatibility
    g.push_back(0);
    g.push_back(0);
    return g;
  }
  //look for , or / 
  string::size_type i = str.find_first_of(",/");
  //extract first allele as portion of string up to first ',' or '/'
  //NOTE: if string consists only of ',' or '/' or anything non-numeric, genotype is taken as missing
  g.push_back(atoi(str.substr(0,i).c_str()));

  if( i != string::npos){// , or / found
      //extract second allele as portion of string after first ',' or '/'
      // NOTE: if nothing after, allele is taken as 0
      g.push_back(atoi(str.substr(i+1,str.length()-i).c_str()));
  }

  return g;  
}

///gets genotypes in admixmap model (hapmix genotypes are coded differently)
void GenotypeLoader::GetGenotype(int i, int SexColumn, const Genome &Loci,  vector<genotype>* genotypes, bool** Missing)
const
{

  unsigned int simplelocus = 0;//simple locus counter
  unsigned complocus = 0;
  unsigned long numhaploid = 0;
  unsigned long numdiploid = 0;
  unsigned long numhaploidX = 0;
  unsigned long numdiploidX = 0;
  unsigned long numObserved = 0;
  //  unsigned numXloci = 0;

  for(unsigned c = 0; c < Loci.GetNumberOfChromosomes(); ++c){
    bool isXchrm = Loci.isXChromosome(c);

    for(unsigned int j = 0; j < Loci.GetSizeOfChromosome(c); ++j){
      genotype G;
      // loop over composite loci to store genotype strings as pairs of integers in stl vector genotype
      const int numLoci = Loci.getNumberOfLoci(complocus);
      unsigned int count = 0;
      //  if(isXchrm)numXloci += numloci;
      for (int locus = 0; locus < numLoci; locus++) {
	const int numalleles = Loci(complocus)->GetNumberOfAllelesOfLocus(locus);

	vector<unsigned short> g = GetGenotype(simplelocus, i, SexColumn);
	if(g.size()==2)
	  if( (g[0] > numalleles) || (g[1] > numalleles))
	    throwGenotypeError(i, simplelocus, Loci(complocus)->GetLabel(0), 
			       g[0], g[1], numalleles );
	  else if (g.size()==1)
	    if( (g[0] > numalleles))
	      throwGenotypeError(i, simplelocus, Loci(complocus)->GetLabel(0), 
				 g[0], 0, numalleles );

	if(isXchrm){
	  if(g.size()==1){
	    ++numhaploidX;
	  }
	  else {//diploid X genotype
	    if(!isFemale(i)){//males cannot have diploid X genotypes
	      //cerr << "Genotype error in Individual " << i << ". Males cannot have diploid X-chromosome genotypes.";
	      //exit(1);
	      //NOTE: allowing this for backward compatibility, for now
	      //instead remove second element
	      g.pop_back();
	    }
	    ++numdiploidX;
	  }
	}
	else{//autosome
	  if(g.size()==1)++numhaploid;
	  else ++numdiploid;
	}
	simplelocus++;
	G.push_back(g);
	count += g[0];
	numObserved += g[0];
      }//end locus loop
      
      Missing[c][j] = (count == 0);
      
      genotypes->push_back(G);
      ++complocus;
    }
  }
  CheckGenotypes(numObserved, numhaploid, numdiploid, numhaploidX, numdiploidX, i, geneticData_[i][0]);
}

void GenotypeLoader::CheckGenotypes(unsigned long numObserved, unsigned long numhaploid, unsigned long numdiploid, 
				    unsigned long numhaploidX, unsigned long numdiploidX, unsigned i, const string& ID)const{

///check for no observed genotypes
    if(numObserved ==0){
        cerr << endl << "WARNING: Individual " << ID << " has no observed genotypes." << endl;
        //exit(1);
    }

/// check for male with diploid X data
  if(numhaploidX + numdiploidX > 0){//some X genotypes present
//     if(numdiploidX>0 && !isFemale(i)){
//       cerr << endl << "Genotype error in Individual " << ID << ": Males cannot have diploid X-chromosome genotypes." << endl;
//       exit(1);
//     }

    if(numhaploidX>0 && isFemale(i)){
///check for female with haploid and diploid X data
      if(numdiploidX >0){
          cerr << endl << "Genotype error in Individual " << ID << ": Females should have diploid X-chromosome genotypes." << endl;
	exit(1);
      }
///check for phased X data but unphased autosomal genotypes
      if(numdiploid>0){
          cerr << endl << "Genotype error in Individual " << ID << ": Female with diploid autosomes and haploid X-Chromosome." << endl; 
	exit(1);
      }
    }
  }
  else{//only autosomes
//check for mixed haploid/diploid data
    if( numhaploid>0 && numdiploid > 0 ){
        cerr << endl << "Genotype error in Individual " << ID << ": Both haploid and diploid genotypes and no X chromosome." << endl;
      exit(1);
    }
  }

}

///write an error message to stderr when a genotype has an invalid allele number
// Doesn't actually throw anything and has to use cerr as logfile is not available yet
void GenotypeLoader::throwGenotypeError(int ind, int locus, std::string label, int g0, int g1, int numalleles)const{

  cerr << "Error in genotypes file:\n"
       << "Individual " << ind << " at locus " << locus <<" (" << label << ")"
       << " has genotype " << g0 << ", " << g1 << " \n"
       << "Number of allelic states at locus = " << numalleles << "\n";
  if(ind == NumIndividuals)
    exit(1);
}

bool GenotypeLoader::isFemale(unsigned i)const{
  //if (options->getgenotypesSexColumn() == 1) {
  int sex = StringConvertor::toInt(geneticData_[i][1]);
  if (sex > 2) {
    cout << "Error: sex must be coded as 0 - missing, 1 - male or 2 - female.\n";
    exit(0);
  }        
  //}
  return (bool)(sex==2);
}

/** Determine if genotype table is in pedfile format. 
    Does this by testing if number of strings in row 1 equals
    twice the number of strings in the header row minus one. 
*/ 
bool GenotypeLoader::determineIfPedFile()const {
  if (geneticData_.size() <= 0) {
    throw string("GenotypeLoader::determineIfPedFile(): geneticData_ has size zero.");
  }
  const bool isPedFile = (bool)(2*geneticData_[0].size() - 1 == geneticData_[1].size());

  return (isPedFile);
}
bool GenotypeLoader::isPedFile()const{
  return IsPedFile;
}


const std::vector<std::string>& GenotypeLoader::getHeader()const{
  return geneticData_[0];
}

void GenotypeLoader::clear(){
  for(unsigned i = 0; i < geneticData_.size(); ++i)
    geneticData_[i].clear();
  geneticData_.clear();

}

