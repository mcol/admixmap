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
#include "DataMatrix.h"

using namespace std;


Genome::Genome()
{
  NumberOfCompositeLoci = 0;
  NumberOfChromosomes = 0;
  LocusArray = 0;
  LengthOfGenome = 0;
  LengthOfXchrm = 0;
  TotalLoci = 0;
  Distances = 0;
  SizesOfChromosomes = 0;
}

Genome::~Genome()
{
  delete[] LocusArray;
  delete[] SizesOfChromosomes;
  delete[] Distances; 
  for(unsigned i = 0; i < NumberOfChromosomes; i++){
    delete C[i];
  }
  delete[] C;
}

//gets contents of locusfile and genotypesfile, creates CompositeLocus array and creates Chromosome array
void Genome::Initialise(const InputData* const data_, int populations, LogWriter &Log, int rank){
  // in parallel version: rank0 needs distances for updating sumintensities
  //                      rank1 needs CompLocus objects
  //                      other ranks need chromosomes and distances

  const DataMatrix& locifileData =  data_->getLocusMatrix();//locus file converted to doubles
  const Vector_s& locusLabels = data_->getLocusLabels();
  bool isWorker = (rank <0 || rank>1);//a worker updates individuals, so needs chromosomes
  
  //determine number of composite loci
  NumberOfCompositeLoci = data_->getNumberOfCompositeLoci();
  //NumberOfChromosomes = data_->getNumberOfChromosomes();
  
  //create array of CompositeLocus objects
  if(rank ==1 || rank==-1)LocusArray = new CompositeLocus[ NumberOfCompositeLoci ];
  if(rank!=1)Distances = new double[ NumberOfCompositeLoci ];
  
  // Set number of alleles at each locus
  unsigned row = 0;//counts lines in locusfile
  vector<unsigned int> cstart;
 
  int cnum = -1; //cnum = number of chromosomes -1
  int lnum = 0;
  LocusTable.resize(NumberOfCompositeLoci);
  X_data = false;
  
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    LocusTable[i].resize(2);
    
    //retrieve first row of this comp locus from locusfile
    const Vector_s& m = data_->getLocusData()[row+1];//+1 because LocusData has a header, LocusMatrix doesn't
    //get chromosome labels from col 4 of locusfile, if there is one   
    if (m.size() == 4) ChrmLabels.push_back(m[3]);

    if(rank!=0)SetDistance( i, locifileData.get( row, 1 ) );//sets distance between locus i and i-1

    if(locifileData.isMissing(row, 1) ){//new chromosome, triggered by missing value
      cnum++;
      lnum = 0; 
      cstart.push_back(i);//locus number of first locus on new chromosome
    } else{
      lnum++;//one more locus on chromosome
    }
    LocusTable[i][0] = cnum;//chromosome on which locus i is located
    LocusTable[i][1] = lnum;//number on chromosome cnum of locus i
    
    //set number of alleles of first locus in comp locus
    if(rank ==1 || rank==-1)LocusArray[i].AddLocus( (int)locifileData.get( row, 0), locusLabels[row] );
    //loop through lines in locusfile for current complocus
    while( row < locifileData.nRows() - 1 && !locifileData.isMissing( row + 1, 1 ) && locifileData.get( row + 1, 1 ) == 0 ){
      if(rank ==1 || rank==-1)LocusArray[i].AddLocus( (int)locifileData.get( row+1, 0 ), locusLabels[row+1] );
      //adds locus with number of alleles given as argument
      row++;
    }
    row++;
    if(rank ==1 || rank==-1){
      if(LocusArray[i].GetNumberOfLoci()>8) Log << "WARNING: Composite locus with >8 loci\n";
      TotalLoci += LocusArray[i].GetNumberOfLoci();
    }
  }//end comp loci loop
  NumberOfChromosomes = cnum +1;
  
  SizesOfChromosomes = new unsigned int[NumberOfChromosomes];//array to store lengths of the chromosomes
  cstart.push_back(NumberOfCompositeLoci);
  for(unsigned c = 0; c < NumberOfChromosomes; ++c) SizesOfChromosomes[c] = cstart[c+1] - cstart[c];
  if(isWorker){
    InitialiseChromosomes(cstart, populations);
  }
  
  PrintSizes(Log);//prints length of genome, num loci, num chromosomes
}

//Creates an array of pointers to Chromosome objects, sets their labels
//also determines length of genome, NumberOfCompositeLoci, TotalLoci, NumberOfChromosomes, SizesOfChromosomes, 
//LengthOfXChrm 
void Genome::InitialiseChromosomes(const vector<unsigned> cstart, int populations){
  C = new Chromosome*[NumberOfChromosomes]; 
  //C is an array of chromosome pointers

  for(unsigned i = 0; i < NumberOfChromosomes; i++){//loop over chromsomes
    int size = cstart[i+1] - cstart[i];//number of loci on chromosome i

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
  }
}

const Chromosome* const* Genome::getChromosomes()const{
  return C;
}
Chromosome* Genome::getChromosome(unsigned j){
  return C[j];
}

//accesses a composite locus
CompositeLocus* Genome::operator() ( int ElementNumber ) const
{
  if ( ElementNumber >= (int)NumberOfCompositeLoci ){
    cout << "WARNING: Genome::operator() Element Number";
    cout << ElementNumber << " > NumberOfCompositeLoci: " << NumberOfCompositeLoci << " accessed." << endl;
  }
  if ( ElementNumber < 0){
    cout << "WARNING: Genome::operator() ElementNumber";
    cout << ElementNumber << " < 0 accessed." << endl;
  }

  return &(LocusArray[ElementNumber]);
}

void Genome::SetDistance( int locus, double distance )
{
  Distances[ locus ] = distance;
}

void Genome::PrintSizes(LogWriter &Log)const{
  Log.setDisplayMode(Quiet);
  Log << "\n" << TotalLoci << " simple loci\n"
      << NumberOfCompositeLoci << " compound loci; "
      << NumberOfChromosomes << " chromosome"; if(NumberOfChromosomes > 1) Log << "s";
  Log << "\n";

  Log << "Effective length of autosomes under study: " << LengthOfGenome << " Morgans.\n";

  if( isX_data() ){
    Log << "Effective length of X chromosome under study: " << LengthOfXchrm << " Morgans.\n";
   }
  Log << "\n";
}

//Accessors
unsigned int Genome::GetNumberOfCompositeLoci()const
{
  return NumberOfCompositeLoci;
}

int Genome::getNumberOfLoci(int j)const{
  return LocusArray[j].GetNumberOfLoci();
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
  return LocusArray[locus].GetNumberOfStates();
}
//returns total number of states accross all comp loci
int Genome::GetNumberOfStates()const
{
  int ret = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    ret += LocusArray[i].GetNumberOfStates();
  }
  return ret;
}

const vector<int> Genome::GetChrmAndLocus( int j )const{
  return LocusTable[j];
}

const vector<vector<int > > Genome::GetChrmAndLocus()const{
  return LocusTable;
}
void Genome::GetChrmAndLocus(unsigned locus, unsigned* c, unsigned* l){
  *c = LocusTable[locus][0];
  *l = LocusTable[locus][1];
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

void Genome::InitialiseLocusCorrelation(const vector<double> rho){
  for( unsigned int j = 0; j < NumberOfChromosomes; j++ ) {
    C[j]->InitialiseLocusCorrelation(rho);
  }
}
// void Genome::InitialiseLocusCorrelation(double rho){
//   for( unsigned int j = 0; j < NumberOfChromosomes; j++ ) {
//     C[j]->InitialiseLocusCorrelation(rho);
//   }
// }

//set global locus correlation across all chromosomes
void Genome::SetLocusCorrelation(const vector<double> rho){
  for( unsigned int j = 0; j < NumberOfChromosomes; j++ ) {
    //in case of global rho model (rho has length 1), sets f globally across loci
    //in hapmixmodel, sets locus-specific f
    C[j]->SetLocusCorrelation(rho, (rho.size()==1), false);
  }
}

void Genome::SetLocusCorrelation(double rho){
  for( unsigned int j = 0; j < NumberOfChromosomes; j++ ) {
    C[j]->SetLocusCorrelation(rho);
  }
}

void Genome::PrintLocusTable(const char* filename)const{
  ofstream outfile(filename);
  outfile << "LocusName\tNumHaps\tMapPosition\tChromosome" << endl;
  unsigned locus = 0;
  for(unsigned c = 0; c < NumberOfChromosomes; ++c){
    double mapPosition = 0.0;
    for(unsigned j  = 0; j < SizesOfChromosomes[c]; ++j){
      mapPosition += Distances[locus];
      outfile << LocusArray[locus].GetLabel(0) << "\t" << LocusArray[locus].GetNumberOfStates() << "\t" 
	      << mapPosition*100.0 << "\t" << c+1 << endl;
      ++locus;
    }
  }
  outfile.close();
}


