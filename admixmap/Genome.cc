/** 
 *   ADMIXMAP
 *   Genome.cc (formerly GeneticArray.cc) 
 *   Class to hold and access (pointers to) Composite Locus objects and information about the genome.
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
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
}


Genome::~Genome()
{
//   if(NumberOfCompositeLoci > 0)
//     for(unsigned int i = 0; i < NumberOfCompositeLoci; i++){
//       delete TheArray[i];
//     }
  delete[] TheArray;
  delete[] SizesOfChromosomes;
  delete[] Distances; 
}

//sets the labels of all composite loci in TheArray using SetLabel function in CompositeLocus
//all loci within each composite locus are assigned the same label from labels vector
void Genome::SetLabels(const vector<string> &labels, vector<double> temp)
{
    int index = -1; // counts through number of composite loci
    int locus = 0;

    for (size_t count = 2; count < labels.size(); ++count) {
        const string& label = labels[count];
        if (temp.size() == 1 || temp[count - 1]) {
            index++;
            locus = 0;
            TheArray[index]->SetLabel(0,label);
        } else if( count != 0 ) {
            locus++;
            TheArray[index]->SetLabel(locus,label);
        }
    }
}

//gets contents of locusfile and genotypesfile and creates CompositeLocus array
void Genome::loadAlleleStatesAndDistances(AdmixOptions *options,InputData *data_){

  DataMatrix locifileData =  data_->getLocusMatrix();
  
  //determine number of composite loci
  NumberOfCompositeLoci = data_->getNumberOfCompositeLoci();
  
  //set up CompositeLocus objects
  InitialiseCompositeLoci();
  
  // Set number of alleles at each locus
  unsigned index =0;
  size_t next_line = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    ++next_line;
    
    const Vector_s& m = data_->getLocusData()[next_line];
    
    //set numbers of alleles and distances for each locus
    if (m.size() == 4)
      ChrmLabels.push_back(m[3]);
    TheArray[i]->SetNumberOfAllelesOfLocus( 0, (int)locifileData.get( i, 0 ) );
    SetDistance( i, locifileData.get( index, 1 ) );
    while( index < locifileData.nRows() - 1 && locifileData.get( index + 1, 1 ) == 0 ){
      ++next_line;
      
      TheArray[i]->AddLocus( (int)locifileData.get( index + 1, 0 ) );
      index++;
    }
    
    TheArray[i]->SetNumberOfLabels();
    index++;
    //Log->logmsg(false,(*Loci)(i)->GetNumberOfLoci());
    //Log->logmsg(false," ");
  }

  // checks of input data files should be in class InputData
  if( options->getTextIndicator() ){

    Vector_s labels = data_->getGeneticData()[0];//header of genotypes file

    vector<double> vtemp = locifileData.getCol(1);
    vtemp.insert(vtemp.begin(), 0.0);// Forces SetLabels method to ignore first row of loci.txt 
    // Add a sex column if it is not included
    if( ! options->getgenotypesSexColumn() ){
      labels.insert(labels.begin(), "\"sexcol\"");
    }
    SetLabels(labels, vtemp);
  }

}

//Creates an array of pointers to Chromosome objects, sets their labels
//also determines length of genome, NumberOfCompositeLoci, TotalLoci, NumberOfChromosomes, SizesOfChromosomes, 
//LengthOfXChrm and chrmandlocus
Chromosome** Genome::GetChromosomes( int populations)
{
  int *cstart = new int[NumberOfCompositeLoci];
  int *cfinish = new int[NumberOfCompositeLoci];
  //since we don't know the number of chromsomes yet, these arrays have the maximum number, num. Comp. Loci, as size

  int cnum = -1; //number of chromosomes -1
  int lnum = 0;
  _chrmandlocus.resize(NumberOfCompositeLoci);
  X_data = false;
  
  TotalLoci = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++){
    
    _chrmandlocus[i].resize(2);
    
    if (GetDistance(i) >= 100){//new chromosome
      cnum++;
      lnum = 0; 
      cstart[cnum] = i; //locus number of first locus on new chromosome
    } else if(cnum==-1){
      cerr << "first locus should have distance of >=100, but doesn't" << endl;
    }else{
      lnum++;
    }
    _chrmandlocus[i][0] = cnum;
    _chrmandlocus[i][1] = lnum;
    cfinish[cnum] = i;//locus number of last locus on currrent chromosome
  }
  
  double LengthOfChrm[cnum+1];
  NumberOfChromosomes = cnum +1;
  SizesOfChromosomes = new unsigned int[NumberOfChromosomes];//array to store lengths of the chromosomes
  
  Chromosome **C = new Chromosome*[cnum+1]; 
  //C is an array of chromosome pointers

  for(int i = 0; i <= cnum; i++){//loop over chromsomes
    int size = cfinish[i] - cstart[i] + 1;//number of loci on chromosome i
    C[i] = new Chromosome(size, cstart[i], populations);
    //C[i] is a pointer to Chromosome

    //set chromosome label
    stringstream label;
    string result;
    //default (if none supplied)
    if( ChrmLabels.size() == 0 ){
      label << "\"" << i << "\"";
      result = label.str();
      C[i]->SetLabel(result);
    }
    else
      C[i]->SetLabel(ChrmLabels[cstart[i]]);

    for(int j = 0; j < size; j++){//loop over loci on chromosome
      //assign pointers to composite locus objects
      //NB: C[i] is a chromosome pointer, *(C[i]) is the chromosome it points to, 
      //(*(C[i]))(j) is an element of the array (of pointers) TheArray ie a locus on the chromosome
      //similarly, (*this)(j) is a locus on the genome
      //LHS and RHS are separate pointers but point to the SAME CompositeLocus objects
      //ie pointers are duplicated but CompositeLocus objects are not.

      (*(C[i]))(j) = (*this)(cstart[i]+j);
      C[i]->SetDistance(j,GetDistance(cstart[i]+j));

      if( j != 0 ){
	string s1("\"X\""), s2("X");
	if( C[i]->GetLabel(0) != s1 && C[i]->GetLabel(0) != s2){
	  LengthOfGenome += GetDistance(cstart[i]+j);
	  LengthOfChrm[i] += GetDistance(cstart[i]+j);
	  //              cout << i << " " << j << " " << GetDistance(cstart[i]+j) << endl;
	  //NB length of genome and LengthOfChrm do not include X chromosome
	}
	//case of X chromosome
	else{
	  LengthOfXchrm += GetDistance(cstart[i]+j);
	  X_data = true;
	  C[i]->ResetStuffForX();
	}
      }
    }
    SizesOfChromosomes[i] = C[i]->GetSize(); 
  }
  
  delete cstart;
  delete cfinish;
  return C;
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
    TheArray[i]->SetNumberOfLoci(1);//not necessary
  }
}

void Genome::SetDistance( int locus, double distance )
{
  Distances[ locus ] = distance;
}

void Genome::SetSizes(LogWriter *Log){
  TotalLoci = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    TotalLoci += TheArray[i]->GetNumberOfLoci();
  }
  Log->logmsg(false, TotalLoci);
  Log->logmsg(false, " simple loci\n");
  Log->logmsg(false, NumberOfCompositeLoci);
  Log->logmsg(false," compound loci; ");
  Log->logmsg(false, NumberOfChromosomes);
  Log->logmsg(false," chromosomes\n");

  Log->logmsg(true,"Effective length of autosomes under study: ");
  Log->logmsg(true, LengthOfGenome);
  Log->logmsg(true," Morgans.\n");

  if( isX_data() ){
    Log->logmsg(true,"Effective length of X chromosome under study: ");
    Log->logmsg(true, LengthOfXchrm);
    Log->logmsg(true," Morgans.\n");
   }
  Log->logmsg(true, "\n");
}

//Accessors
unsigned int Genome::GetNumberOfCompositeLoci()
{
  return NumberOfCompositeLoci;
}

int Genome::getNumberOfLoci(int j){
  return TheArray[j]->GetNumberOfLoci();
}
unsigned int Genome::GetNumberOfChromosomes(){
  return NumberOfChromosomes;
}

//returns total number of simple loci
unsigned int Genome::GetTotalNumberOfLoci(){
  return TotalLoci;
}
//returns int array of chromosome sizes
unsigned int *Genome::GetSizesOfChromosomes(){
  return SizesOfChromosomes;
}

double *Genome::GetDistances()
{
  return( Distances );
}

double Genome::GetDistance( int locus )
{
  return( Distances[ locus ] );
}

//returns number of states of a comp locus
int Genome::GetNumberOfStates(int locus){
  return TheArray[locus]->GetNumberOfStates();
}
//returns total number of states accross all comp loci
int Genome::GetNumberOfStates()
{
  int ret = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    ret += TheArray[i]->GetNumberOfStates();
  }
  return ret;
}

vector<int> Genome::GetChrmAndLocus( int j ){
  return _chrmandlocus[j];
}

vector<vector<int > > Genome::GetChrmAndLocus(){
  return _chrmandlocus;
}
bool Genome::isX_data()
{
   return X_data;
}

double Genome::GetLengthOfGenome()
{
   return LengthOfGenome;
}
double Genome::GetLengthOfXchrm()
{
   return LengthOfXchrm;
}








