/** 
 *   ADMIXMAP
 *   Genome.cc (formerly GeneticArray.cc) 
 *   Class to hold and access (pointers to) Composite Locus objects and information about the genome.
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include <stdlib.h>
#include <sstream>
#include "Genome.h"
#include "Chromosome.h"
#include "DataMatrix.h"

using namespace std;


Genome::Genome()
{
  NumberOfCompositeLoci = 0;
  NumberOfChromosomes = 0;
  TheArray = 0;
  LengthOfGenome = 0;
  LengthOfXchrm = 0;
  TotalLoci = 0;
  Distances = 0;
  SizesOfChromosomes = 0;
}

//used by Chromosome
Genome::Genome( int size )
{
  if ( size < 1 ){
    size = 1;
  }
  NumberOfCompositeLoci = size;
  TheArray = new CompositeLocus*[ NumberOfCompositeLoci ];
  Distances = new double[ NumberOfCompositeLoci ];
  LengthOfGenome = 0;
  LengthOfXchrm = 0;
  TotalLoci = 0;
  SizesOfChromosomes = 0;
}


Genome::~Genome()
{
  if(NumberOfCompositeLoci > 0) 
    //The CompositeLocus objects are owned by a pure Genome object (Loci)
    //Chromosomes only have an array of pointers so they are not allowed to delete any CompositeLoci
    //otherwise the destruction of the lead Genome object would trigger an illegal second attempt at their destruction
    for(unsigned int i = 0; i < NumberOfCompositeLoci; i++){
      delete TheArray[i];
    }
  delete[] TheArray;
  delete[] SizesOfChromosomes;
  delete[] Distances; 
  for(unsigned i = 0; i < NumberOfChromosomes; i++){
    delete C[i];
  }
  delete[] C;
}

//sets the labels of all composite loci in TheArray using SetLabel function in CompositeLocus
//all loci within each composite locus are assigned the same label from labels vector
void Genome::SetLabels(const vector<string> &labels, const vector<double> &distances)
{
    int index = -1; // counts through composite loci
    int locus = 0;//counts through loci on a comp locus

    for (size_t count = 0; count < labels.size(); ++count) {
        const string& label = labels[count];
        if (distances[count] > 0.0) {//new comp locus
            locus = 0;
	    ++index;
        } 
	TheArray[index]->SetLabel(locus++,label);
    }
}

//gets contents of locusfile and genotypesfile and creates CompositeLocus array
void Genome::loadAlleleStatesAndDistances(const InputData* const data_){

  DataMatrix locifileData =  data_->getLocusMatrix();//locus file converted to doubles
  
  //determine number of composite loci
  NumberOfCompositeLoci = data_->getNumberOfCompositeLoci();
  
  //set up CompositeLocus objects
  InitialiseCompositeLoci();
  
  // Set number of alleles at each locus
  unsigned row = 0;//counts lines in locusfile
  size_t next_line = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    ++next_line;
    
    //get chromosome labels from col 4 of locusfile, if there is one    
    const Vector_s& m = data_->getLocusData()[next_line];
    if (m.size() == 4) ChrmLabels.push_back(m[3]);

    //set numbers of alleles and distances for each locus
    TheArray[i]->AddLocus( (int)locifileData.get( row, 0 ) );//sets number of alleles of first locus
    SetDistance( i, locifileData.get( row, 1 ) );//sets distance between locus i and i-1

    //loop through lines in locusfile for current complocus
    while( row < locifileData.nRows() - 1 && locifileData.get( row + 1, 1 ) == 0 ){
      ++next_line;
      TheArray[i]->AddLocus( (int)locifileData.get( row + 1, 0 ) );//adds locus with number of alleles given as argument
      row++;
    }
    
    TheArray[i]->SetNumberOfLabels();
    row++;
    if(TheArray[i]->GetNumberOfLoci()>8) cerr << "WARNING: Composite locus with >8 loci\n";
  }

  SetLabels(data_->getLocusLabels(), locifileData.getCol(1));
}

//Creates an array of pointers to Chromosome objects, sets their labels
//also determines length of genome, NumberOfCompositeLoci, TotalLoci, NumberOfChromosomes, SizesOfChromosomes, 
//LengthOfXChrm and chrmandlocus
void Genome::GetChromosomes( int populations)
{
  int *cstart = new int[NumberOfCompositeLoci];
  int *cfinish = new int[NumberOfCompositeLoci];
  //since we don't know the number of chromsomes yet, these arrays have the maximum number, num. Comp. Loci, as size
  //TODO: should use vectors and expand as necessary

  int cnum = -1; //cnum = number of chromosomes -1
  int lnum = 0;
  LocusTable.resize(NumberOfCompositeLoci);
  X_data = false;
  
  TotalLoci = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++){
    
    LocusTable[i].resize(2);
    
    if (GetDistance(i) >= 100){//new chromosome
      cnum++;
      lnum = 0; 
      cstart[cnum] = i; //locus number of first locus on new chromosome
    } else if(cnum==-1){
      cerr << "first locus should have distance of >=100, but doesn't" << endl;
    }else{
      lnum++;
    }
    LocusTable[i][0] = cnum;//chromosome on which locus i is located
    LocusTable[i][1] = lnum;//number on chromosome cnum of locus i
    cfinish[cnum] = i;//locus number of last locus on currrent chromosome
  }
  
  NumberOfChromosomes = cnum +1;
  SizesOfChromosomes = new unsigned int[NumberOfChromosomes];//array to store lengths of the chromosomes
  
  C = new Chromosome*[cnum+1]; 
  //C is an array of chromosome pointers

  for(unsigned i = 0; i < NumberOfChromosomes; i++){//loop over chromsomes
    int size = cfinish[i] - cstart[i] + 1;//number of loci on chromosome i

    //set chromosome label
    string label;
    //default (if none supplied), numbered in sequence from 1
    if( ChrmLabels.size() == 0 ){
      stringstream labelstr;
      labelstr << "\"" << i+1 << "\"";
      label = labelstr.str();
    }
    else
      label = ChrmLabels[cstart[i]];
    //determine if X Chromosome
    bool isX = false;
    string s1("\"X\""), s2("X");
    isX = ( (label == s1) || (label == s2));

    C[i] = new Chromosome(size, cstart[i], populations, isX);
    //C[i] is a pointer to Chromosome

    C[i]->SetLabel(label);

    for(int j = 0; j < size; j++){//loop over loci on chromosome

      C[i]->SetDistance(j,GetDistance(cstart[i]+j));

      if( j != 0 ){
	if( !isX ){
	  LengthOfGenome += GetDistance(cstart[i]+j);
	  //              cout << i << " " << j << " " << GetDistance(cstart[i]+j) << endl;
	  //NB length of genome does not include X chromosome
	}
	//case of X chromosome
	else{
	  LengthOfXchrm += GetDistance(cstart[i]+j);
	  X_data = true;
	  XChromosomeIndex = cstart[i];
	}
      }
    }
    SizesOfChromosomes[i] = C[i]->GetSize(); 
  }
  
  delete[] cstart;
  delete[] cfinish;
}

const Chromosome* const* Genome::getChromosomes()const{
  return C;
}
Chromosome* Genome::getChromosome(unsigned j){
  return C[j];
}

//accesses a composite locus
CompositeLocus*& Genome::operator() ( int ElementNumber ) const
{
  if ( ElementNumber >= (int)NumberOfCompositeLoci ){
    cout << "WARNING: Genome::operator() Element Number";
    cout << ElementNumber << " > NumberOfCompositeLoci: " << NumberOfCompositeLoci << " accessed." << endl;
  }
  if ( ElementNumber < 0){
    cout << "WARNING: Genome::operator() ElementNumber";
    cout << ElementNumber << " < 0 accessed." << endl;
  }

  return TheArray[ElementNumber];
}

//creates array of CompositeLocus objects
void Genome::InitialiseCompositeLoci()
{
  TheArray = new CompositeLocus*[ NumberOfCompositeLoci ];
  Distances = new double[ NumberOfCompositeLoci ];
  for(unsigned int i = 0; i < NumberOfCompositeLoci;i++){
    //each element in the array is a pointer to a CompositeLocus object
    TheArray[i] = new CompositeLocus();
  }
}

void Genome::SetDistance( int locus, double distance )
{
  Distances[ locus ] = distance;
}

void Genome::SetSizes(LogWriter &Log){
  TotalLoci = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    TotalLoci += TheArray[i]->GetNumberOfLoci();
  }
  Log.setDisplayMode(Quiet);
  Log << "\n" << TotalLoci << " simple loci\n"
      << NumberOfCompositeLoci << " compound loci; "
      << NumberOfChromosomes << " chromosome"; if(NumberOfChromosomes > 1) Log << "s";
  Log << "\n";

  Log << "Effective length of autosomes under study: " << LengthOfGenome << " Morgans.\n";

  if( isX_data() ){
    Log << "Effective length of X chromosome under study: " << LengthOfXchrm << " Morgans.\n";
   }
}

//Accessors
unsigned int Genome::GetNumberOfCompositeLoci()const
{
  return NumberOfCompositeLoci;
}

int Genome::getNumberOfLoci(int j)const{
  return TheArray[j]->GetNumberOfLoci();
}
unsigned int Genome::GetNumberOfChromosomes()const{
  return NumberOfChromosomes;
}

//returns total number of simple loci
unsigned int Genome::GetTotalNumberOfLoci()const{
  return TotalLoci;
}
//returns int array of chromosome sizes
const unsigned int *Genome::GetSizesOfChromosomes()const{
  return SizesOfChromosomes;
}
unsigned Genome::GetSizeOfChromosome(unsigned j)const{
  return SizesOfChromosomes[j];
}
const double *Genome::GetDistances()const
{
  return( Distances );
}

double Genome::GetDistance( int locus )const
{
  return( Distances[ locus ] );
}

//returns number of states of a comp locus
int Genome::GetNumberOfStates(int locus)const{
  return TheArray[locus]->GetNumberOfStates();
}
//returns total number of states accross all comp loci
int Genome::GetNumberOfStates()const
{
  int ret = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    ret += TheArray[i]->GetNumberOfStates();
  }
  return ret;
}

const vector<int> Genome::GetChrmAndLocus( int j )const{
  return LocusTable[j];
}

const vector<vector<int > > Genome::GetChrmAndLocus()const{
  return LocusTable;
}
bool Genome::isX_data()const
{
   return X_data;
}

double Genome::GetLengthOfGenome()const
{
   return LengthOfGenome;
}
double Genome::GetLengthOfXchrm()const
{
   return LengthOfXchrm;
}
unsigned Genome::isXChromosome(unsigned j){
  return C[j]->isXChromosome();
}

unsigned Genome::getFirstXLocus()const{
  if(X_data)return XChromosomeIndex;
  else return NumberOfCompositeLoci;
}

void Genome::InitialiseLociCorr(const vector<double> rho){
  for( unsigned int j = 0; j < NumberOfChromosomes; j++ ) {
    C[j]->InitialiseLociCorr(rho);
  }
}
void Genome::SetLociCorr(const vector<double> rho){
  for( unsigned int j = 0; j < NumberOfChromosomes; j++ ) {
    C[j]->SetLociCorr(rho);
  }
}
void Genome::InitialiseLociCorr(double rho){
  for( unsigned int j = 0; j < NumberOfChromosomes; j++ ) {
    C[j]->InitialiseLociCorr(rho);
  }
}
void Genome::SetLociCorr(double rho){
  for( unsigned int j = 0; j < NumberOfChromosomes; j++ ) {
    C[j]->SetLociCorr(rho);
  }
}




