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
#include "bclib/LogWriter.h"
#include "Genome.h"
#include "bclib/DataReader.h"
#include "bclib/StringConvertor.h"
#include "bclib/LogWriter.h"

//delimiters for diploid genotypes
#define GENOTYPE_DELIMS ",/:"

GenotypeLoader::GenotypeLoader(){
  NumIndividuals = 0;
  numDiploid = 0;
  SexColumn = 0;
}

void GenotypeLoader::Read(const char* filename, unsigned NumLociInLocusfile, bclib::LogWriter& Log){
  bclib::DataReader::ReadData(filename, geneticData_, Log);
  IsPedFile = determineIfPedFile();
  DetermineSexColumn( NumLociInLocusfile, Log);
  NumIndividuals = geneticData_.size() - 1;
}

unsigned GenotypeLoader::getNumberOfIndividuals()const{
  return NumIndividuals;
}

unsigned GenotypeLoader::NumLoci()const{
  return geneticData_[0].size();
}


///checks number of loci in genotypes file is the same as in locusfile, 
///and determines if there is a sex column
void GenotypeLoader::DetermineSexColumn(unsigned NumLociInLocusfile, bclib::LogWriter& Log){
   // Determine if "Sex" column present in genotypes file.
  if (NumLociInLocusfile == this->NumLoci() - 1) {
    SexColumn = 0;//no sex col
  } else if (NumLociInLocusfile == this->NumLoci() - 2) {
    SexColumn  = 1;//sex col
  } else {//too many cols
    Log << "Error: " << NumLociInLocusfile << " loci in locus file but " 
	 <<  this->NumLoci() - 1 << " loci in genotypes file.\n";
    exit(1);
  }
}

int GenotypeLoader::getSexColumn()const{
  return SexColumn;
}
/**
   extracts genotype for given individual at given locus
   \param locus the locus index
   \param individual individual index (row number)
   \param SexColumn index of sex column, or 0 if none
   \return genotype as vector of unsigned short integers
*/
vector<unsigned short> GenotypeLoader::GetGenotype(unsigned locus, int individual)const{
  int col = 1 + SexColumn + locus;
  if (IsPedFile)col = 1 + SexColumn + 2*locus;
  return GetGenotype(geneticData_[individual][col]);
}


///convert a genotype string to a vector of unsigned ints
vector<unsigned short> GenotypeLoader::GetGenotype(const string genostring)const{
  vector<unsigned short> g;
 
  //strip quotes from string
  const std::string str = bclib::StringConvertor::dequote(genostring);
  if(str.length()==0){
    //if empty string, interpret as missing genotype for backward compatibility
    g.push_back(0);
    g.push_back(0);
    return g;
  }
  //look for , or / or :
  string::size_type i = str.find_first_of(GENOTYPE_DELIMS);
  //extract first allele as portion of string up to first ',' or '/' or ':'
  //NOTE: if string consists only of ',' or '/' or ':' or anything non-numeric, genotype is taken as missing
  g.push_back(atoi(str.substr(0,i).c_str()));

  if( i != string::npos){// , or / found
      //extract second allele as portion of string after first ',' or '/'
      // NOTE: if nothing after, allele is taken as 0
      g.push_back(atoi(str.substr(i+1,str.length()-i).c_str()));
  }

  return g;  
}

///gets genotypes in admixmap model (hapmix genotypes are coded differently)
void GenotypeLoader::GetGenotype(int i, const Genome &Loci,  vector<genotype>* genotypes, bool** Missing)
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

	vector<unsigned short> g = GetGenotype(simplelocus, i);
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
  if (SexColumn == 1){
    int sex = bclib::StringConvertor::toInt(geneticData_[i][1]);
    if (sex > 2) {
      cout << "Error: sex must be coded as 0 - missing, 1 - male or 2 - female.\n";
      exit(0);
    }        
    
    return (bool)(sex==2);
  }
  else return false;
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

bool GenotypeLoader::CheckForUnobservedAlleles(const bclib::DataMatrix& LocusData, bclib::LogWriter& Log){
  //sanity checks
  if(!(geneticData_.size()))
    throw runtime_error("CheckForUnobservedAlleles called when genotype data are not available\n");
  if(LocusData.nRows() != geneticData_[0].size()-1)
    throw string("Error in CheckForUnobservedAlleles: number of loci in locusfile and genotypesfile do not match");

  //offset cols by 2 if sexcol, 1 (for ID) otherwise
  const unsigned offset = (SexColumn>0) ? 2 :1;
  vector<string> MMLoci;//vector of names of monomorphic loci
  for(unsigned locus = 0; locus < LocusData.nRows(); ++locus){
    //create vector of counts of length number of alleles
    vector<unsigned> AlleleCounts((unsigned)(LocusData.get(locus, 0)), 0);

    //loop over individuals, counting alleles
    for(unsigned i = 1; i < geneticData_.size(); ++i){
      const string g = bclib::StringConvertor::dequote(geneticData_[i][locus+offset]);
      //find separator in diploid genotypes
      string::size_type sep = g.find_first_of(GENOTYPE_DELIMS);
      const unsigned allele1 = atoi(g.substr(0, sep).c_str());
      const unsigned allele2 = atoi(g.substr(sep+1).c_str());
      //check allele codes do not exceed number of alleles of this locus
      if(allele1 > AlleleCounts.size() 
          || (sep != string::npos && allele2 > AlleleCounts.size()))
          throwGenotypeError(i, locus+1, geneticData_[0][locus+1], 
                             allele1, allele2, AlleleCounts.size());
      //increment count of first allele
      if(allele1)++AlleleCounts[allele1-1];
      //increment count of second allele
      if(sep != string::npos)      
	if(allele2)++AlleleCounts[allele2-1 ];
    }
    for(vector<unsigned>::const_iterator a = AlleleCounts.begin(); a != AlleleCounts.end(); ++a){
      if(*a ==0 ){//found an unobserved allele
	MMLoci.push_back(geneticData_[0][locus+1]);
	break;
      }
    }

  }

  bool allClear = true;
  if(MMLoci.size()){
    allClear = false;
    Log << "WARNING: found " << (int)(MMLoci.size());
    if(MMLoci.size()==1) Log << " locus "; else Log << " loci ";
    Log << "with unobserved alleles:\n";
    if(MMLoci.size() < 10){
      for(vector<string>::const_iterator i = MMLoci.begin(); i < MMLoci.end(); ++i)
	Log << *i << "\n";
    }
    else
      Log << "The first is " << MMLoci[0] << " and the last is " << MMLoci[MMLoci.size()-1] << "\n";
  }
  return allClear;
}
