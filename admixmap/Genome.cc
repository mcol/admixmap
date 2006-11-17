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
#include "utils/DataMatrix.h"
#include "utils/StringConvertor.h"
#include "Comms.h"
using namespace std;


Genome::Genome()
{
  NumberOfCompositeLoci = 0;
  NumberOfChromosomes = 0;
  LocusArray = 0;
  LengthOfGenome = 0;
  LengthOfXchrm = 0;
  XChromosomeIndex = 0; 
  TotalLoci = 0;
  Distances = 0;
  SizesOfChromosomes = 0;
  C = 0;
}

Genome::~Genome()
{
  delete[] LocusArray;
  delete[] SizesOfChromosomes;
  delete[] Distances; 
  if(C){
      for(unsigned i = 0; i < NumberOfChromosomes; i++){
	  delete C[i];
      }
      delete[] C;
  }
}

///gets contents of locusfile and creates CompositeLocus array and Chromosome array
void Genome::Initialise(const InputData* const data_, int populations, LogWriter &Log){
  // in parallel version: master(rank0) needs distances for updating sumintensities
  //                      freqsampler (rank1) needs CompLocus objects
  //                      workers (rank >1) need chromosomes and distances

  const DataMatrix& locifileData =  data_->getLocusMatrix();//locus file converted to doubles
  const Vector_s& locusLabels = data_->getLocusLabels();
  const bool isMaster = Comms::isMaster();
  const bool isWorker = Comms::isWorker();//a worker updates individuals, so needs chromosomes
  const bool isFreqUpdater = Comms::isFreqSampler();;
  
  //determine number of composite loci
  NumberOfCompositeLoci = data_->getNumberOfCompositeLoci();
  TotalLoci = data_->getNumberOfSimpleLoci();
  
  //create array of CompositeLocus objects
  if(isFreqUpdater)LocusArray = new CompositeLocus[ NumberOfCompositeLoci ];
  if(isMaster || isWorker)Distances = new double[ NumberOfCompositeLoci ];
  
  // Set number of alleles at each locus
  unsigned row = 0;//counts lines in locusfile
  vector<unsigned int> cstart;
 
  int cnum = -1; //cnum = number of chromosomes -1
  int lnum = 0;
  LocusTable.resize(NumberOfCompositeLoci);
  X_data = false;

  //determine if distances are given in Morgans or centimorgans
  GeneticDistanceUnit unit = data_->getUnitOfDistance();
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    LocusTable[i].resize(2);
    
    //retrieve first row of this comp locus from locusfile
    const Vector_s& m = data_->getLocusData()[row+1];//+1 because LocusData has a header, LocusMatrix doesn't
    //get chromosome labels from col 4 of locusfile, if there is one   
    if (m.size() == 4) ChrmLabels.push_back(StringConvertor::dequote(m[3]));

    if(isMaster || isWorker){
      Distances[ i ] = locifileData.get( row, 1 );
      if(unit == centimorgans)Distances[i] /= 100.0;//convert to Morgans
      //      SetDistance( i, locifileData.get( row, 1 ) );//sets distance between locus i and i-1
    }

    if(locifileData.isMissing(row, 1) || locifileData.get(row, 1)>=100.0){//new chromosome, triggered by missing value or value of >=100 for distance
      cnum++;
      lnum = 0; 
      cstart.push_back(i);//locus number of first locus on new chromosome
    } else{
      lnum++;//one more locus on chromosome
    }
    LocusTable[i][0] = cnum;//chromosome on which locus i is located
    LocusTable[i][1] = lnum;//number on chromosome cnum of locus i
    
    //set number of alleles of first locus in comp locus
    if(isFreqUpdater)LocusArray[i].AddLocus( (int)locifileData.get( row, 0), locusLabels[row] );
    //loop through lines in locusfile for current complocus
    while( row < locifileData.nRows() - 1 && !locifileData.isMissing( row + 1, 1 ) && locifileData.get( row + 1, 1 ) == 0 ){
      if(isFreqUpdater)LocusArray[i].AddLocus( (int)locifileData.get( row+1, 0 ), locusLabels[row+1] );
      //adds locus with number of alleles given as argument
      row++;
    }
    row++;
    if(isFreqUpdater){
#ifdef PARALLEL 
//at present, parallel version can only handle diallelic loci
      if(LocusArray[i].GetNumberOfStates()>2){
	stringstream err;
	err << "sorry, I can only handle diallelic loci. Composite locus " << i+1 << " has " 
	    << LocusArray[i].GetNumberOfLoci() << " simple loci and " << LocusArray[i].GetNumberOfStates() << " alleles/haplotypes";
	throw err.str();
      }
#endif
      if(LocusArray[i].GetNumberOfLoci()>8) Log << "WARNING: Composite locus with >8 loci\n";
      //TotalLoci += LocusArray[i].GetNumberOfLoci();
    }
  }//end comp loci loop

  NumberOfChromosomes = cnum +1;
  SizesOfChromosomes = new unsigned int[NumberOfChromosomes];//array to store lengths of the chromosomes
  cstart.push_back(NumberOfCompositeLoci);//add extra element for next line to work
  for(unsigned c = 0; c < NumberOfChromosomes; ++c) SizesOfChromosomes[c] = cstart[c+1] - cstart[c];
  if(isWorker){//create Chromosome objects
    InitialiseChromosomes(cstart, populations);
  }
  

  if(isMaster || isWorker ){
    PrintSizes(Log, unit);//prints length of genome, num loci, num chromosomes
  }
}

///Creates an array of pointers to Chromosome objects and sets their labels.
///Also determines length of genome, NumberOfCompositeLoci, TotalLoci, NumberOfChromosomes, SizesOfChromosomes, 
///LengthOfXChrm 
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
    string s1("X"), s2("x");
    isX = ( (label == s1) || (label == s2) );
    if(isX){
      X_data = true;
      XChromosomeIndex = cstart[i];//index of first locus on X chromosome
    }

    C[i] = new Chromosome(i, size, cstart[i], populations, isX);
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
	}
      }
    }
  }
}
///accesses the entire chromosome array
const Chromosome* const* Genome::getChromosomes()const{
  return C;
}
///accesses a chromosome
Chromosome* Genome::getChromosome(unsigned j){
  return C[j];
}

///accesses a composite locus
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

/// Writes numbers of loci and chromosomes and length of genome to Log and screen.
/// unit is the unit of measurement of the distances in the locusfile (Morgans/centiMorgans) 
void Genome::PrintSizes(LogWriter &Log, GeneticDistanceUnit u)const{
#ifdef PARALLEL
  ///1st worker tells master length of autosomes and xchrm
  ///(this is determined during creation of chromosomes, which master doesn't do)
  ///only master is allowed to write to logfile
  const int rank = MPI::COMM_WORLD.Get_rank();
  if(rank == 2) {
    MPI::COMM_WORLD.Send(&LengthOfGenome, 1, MPI::DOUBLE, 0, 0);
    MPI::COMM_WORLD.Send(&LengthOfXchrm, 1, MPI::DOUBLE, 0, 1);
  }
  if(rank==0){
    MPI::Status status;
    MPI::COMM_WORLD.Recv((double*)&LengthOfGenome, 1, MPI::DOUBLE, 2, 0, status);
    MPI::COMM_WORLD.Recv((double*)&LengthOfXchrm, 1, MPI::DOUBLE, 2, 1, status);
  }
#endif
  
  Log.setDisplayMode(Quiet);
  Log << "\n" << TotalLoci << " simple loci\n"
      << NumberOfCompositeLoci << " compound loci; "
      << NumberOfChromosomes << " chromosome"; if(NumberOfChromosomes > 1) Log << "s";
  Log << "\n";

  string unitstring;
  switch(u){
      case centimorgans:{
	  unitstring = " centimorgans";
	  break;
      }
      case Morgans:{
	  unitstring = " Morgans";
	  break;
      }
      case megabases:{
	  unitstring = " megabases";
	  break;
      }
      default:{
	  Log << "[unsupported unit]\n";
	  exit(1);
      }
  }

  Log << "Effective length of autosomes under study: ";
  if(u == centimorgans)Log << LengthOfGenome*100.0 ;
  else Log << LengthOfGenome;
  Log << unitstring << ".\n";

  if( isX_data() ){
    Log << "Effective length of X chromosome under study: ";
    if(u == centimorgans)Log << LengthOfXchrm*100.0;
    else Log << LengthOfXchrm;
    Log << unitstring << ".\n";


   }
  Log << "\n";
}

//Accessors

/// returns the number of composite loci
unsigned int Genome::GetNumberOfCompositeLoci()const
{
  return NumberOfCompositeLoci;
}
///returns the number of loci in a given composite locus
int Genome::getNumberOfLoci(int 
#ifdef PARALLEL
    )const{
  return 1;
#else
  j)const{
  return LocusArray[j].GetNumberOfLoci();
#endif
}
///returns the number of chromosomes
unsigned int Genome::GetNumberOfChromosomes()const{
  return NumberOfChromosomes;
}

///returns total number of simple loci
unsigned int Genome::GetTotalNumberOfLoci()const{
  return TotalLoci;
}
///returns int array of chromosome sizes
const unsigned int *Genome::GetSizesOfChromosomes()const{
  return SizesOfChromosomes;
}
///returns the number of loci on a given chromosome
unsigned Genome::GetSizeOfChromosome(unsigned j)const{
  return SizesOfChromosomes[j];
}
/// returns the vector of distances between loci
const double *Genome::GetDistances()const
{
  return( Distances );
}
///returns distance between a given locus and the previous one
double Genome::GetDistance( int locus )const
{
  return( Distances[ locus ] );
}

///returns number of states of a comp locus
int Genome::GetNumberOfStates(int
#ifdef PARALLEL
    )const{
    return 2;
#else
    locus)const{
  return LocusArray[locus].GetNumberOfStates();
#endif
}
///returns total number of states accross all comp loci
int Genome::GetNumberOfStates()const
{
#ifdef PARALLEL
  return NumberOfCompositeLoci*2;
#else
  int ret = 0;
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++ ){
    ret += LocusArray[i].GetNumberOfStates();
  }
  return ret;
#endif
}

unsigned Genome::GetChrNumOfLocus(unsigned locus){
  return LocusTable[locus][0];
}

const vector<int> Genome::GetChrmAndLocus( int j )const{
  return LocusTable[j];
}

const vector<vector<int > > Genome::GetChrmAndLocus()const{
  return LocusTable;
}
/// For a given (composite) locus, returns the number of the chromosome it is on and the position on that chromosome.
void Genome::GetChrmAndLocus(unsigned locus, unsigned* c, unsigned* l){
  *c = LocusTable[locus][0];
  *l = LocusTable[locus][1];
}
///indicates whether there is an X chromosome
bool Genome::isX_data()const
{
   return X_data;
}
///returns length of Genome in Morgans
double Genome::GetLengthOfGenome()const
{
   return LengthOfGenome;
}
//returns length of X chromosome in Morgans
double Genome::GetLengthOfXchrm()const
{
   return LengthOfXchrm;
}
///indicates if a chromosome is an X chromosome
unsigned Genome::isXChromosome(unsigned j)const{
  return C[j]->isXChromosome();
}
bool Genome::isXLocus(unsigned j)const{
  return C[LocusTable[j][0]]->isXChromosome();
}
/// returns index of X chromosome
unsigned Genome::getFirstXLocus()const{
  if(X_data)return XChromosomeIndex;
  else return NumberOfCompositeLoci;
}

///set global locus correlation across all chromosomes, case of vector-valued rho
void Genome::SetLocusCorrelation(const vector<double> rho){
  if(rho.size()==1) 
    for( unsigned int j = 0; j < NumberOfChromosomes; j++ ) {
      //in case of global rho model (rho has length 1), sets f globally across loci
      C[j]->SetLocusCorrelation(rho, true, false);
  }
  else{      //in hapmixmodel, sets locus-specific f
    if(rho.size()<NumberOfCompositeLoci-NumberOfChromosomes)throw string("Bad arguments passed to Chromosome::SetLocusCorr");
    vector<double>::const_iterator rho_iter = rho.begin();
    for( unsigned int j = 0 ; j < NumberOfChromosomes; j++ ) {
      C[j]->SetLocusCorrelation(rho_iter);
      rho_iter += C[j]->GetSize()-1;
    }
  }
}

///set global locus correlation across all chromosomes, case of global rho
void Genome::SetLocusCorrelation(double rho){
  for( unsigned int j = 0; j < NumberOfChromosomes; j++ ) {
    C[j]->SetLocusCorrelation(rho);
  }
}

///Prints table of cpmposite loci for R script to read
void Genome::PrintLocusTable(const char* filename, const vector<double>& Dist)const{
  //could use Distances array member but in parallel version the processor calling this function will not have this array
  //so we use the raw distances from the locusfile instead
  ofstream outfile(filename);
  outfile << "LocusName\tNumHaps\tMapPosition\tChromosome" << endl;
  unsigned locus = 0;//counter for composite locus
  unsigned simple_locus = 0;//need to count simple loci to step through vector of distances
  for(unsigned c = 0; c < NumberOfChromosomes; ++c){
    double mapPosition = 0.0;
    //first locus on chromosome
    outfile << LocusArray[locus].GetLabel(0) << "\t" << LocusArray[locus].GetNumberOfStates() << "\t" 
	    << 0.0 << "\t" << c+1 << endl;
    simple_locus += LocusArray[locus].GetNumberOfLoci();
    ++locus;
    for(unsigned j  = 1; j < SizesOfChromosomes[c]; ++j){//step through rest of loci on chromosome
      mapPosition += Dist[simple_locus];//increment map position by distance of first locus in complocus
      outfile << LocusArray[locus].GetLabel(0) << "\t" << LocusArray[locus].GetNumberOfStates() << "\t" 
	      << mapPosition*100.0 << "\t" << c+1 << endl;
      simple_locus += LocusArray[locus].GetNumberOfLoci();
      ++locus;
    }
  }
  outfile.close();
}


