/** 
 *   ADMIXMAP
 *   Chromosome.cc 
 *   Class to perform chromosome-wise updates and implement HMM class.
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "Chromosome.h"
#include "Individual.h"
#include <iostream>
#include "functions.h"

#define PR(x) cerr << #x << " = " << x << endl;

using namespace std;

Chromosome::Chromosome(){
    Distances = 0;
    NumberOfCompositeLoci = 0;
    populations = 0;
    isX = false;
    Diploid = true;
    f = 0;
    CodedStates = 0;
}

Chromosome::Chromosome(int size, int start, int inpopulations, bool isx = false) 
		      //size = number of comp loci on chromosome
{
  _startLocus = start;
  populations = inpopulations;
  isX = isx;
  if ( size < 1 ){
    size = 1;
  }
  NumberOfCompositeLoci = size;
  Distances = new double[ NumberOfCompositeLoci ];

  SampleStates.SetDimensions( size, populations );

  CodedStates = new int[size];
  f = new double[2*size];
}

Chromosome::~Chromosome()
{
  delete[] CodedStates;
  delete[] f;
  delete[] Distances; 
}

// ******** Chromosome information **********************************
bool Chromosome::isXChromosome()const
{
  return isX;
}

void Chromosome::SetLabel( string label )
{//sets label for this chromosome, usually "1", "2", etc or "X"
   _Label = label;
}

const string Chromosome::GetLabel( )const
{
   return _Label;
}

void Chromosome::SetDistance( int locus, double distance )
{
  Distances[ locus ] = distance;
}
const double *Chromosome::GetDistances()const
{
  return( Distances );
}

double Chromosome::GetDistance( int locus )const
{//returns distance between a locus and the previous locus
  return( Distances[ locus ] );
}

// Returns the number of the num'th compositelocus on this chromosome
// eg if chromosome 2 consists of compositeloci 5, 6, 7 and 8,
// GetLocus(i) returns 5 + i
int Chromosome::GetLocus(int num)const{
  return _startLocus + num;
}
//??: remove either of these functions
//returns number of composite loci in the chromosome
unsigned int Chromosome::GetSize()const{
  return NumberOfCompositeLoci;
}
unsigned int Chromosome::GetNumberOfCompositeLoci()const
{
  return NumberOfCompositeLoci;
}

// ****************** Setting of locus correlation, f *************************

//Initialises locus ancestry correlations f.  Necessary since individual-level parameters updated before global rho (in Latent)
void Chromosome::InitialiseLocusCorrelation(const double rho_){
  double rho = rho_;
  if(isX)rho *= 0.5;
  f[0] = f[1] = 0.0;
  for(unsigned int j = 1; j < NumberOfCompositeLoci; j++ )
    f[2*j] = f[2*j + 1] = ( -GetDistance( j ) * rho > -700) ? exp( -GetDistance( j ) * rho ) : 0.0;
}
void Chromosome::InitialiseLocusCorrelation(const vector<double> rho_){
  if(rho_.size() == 1)InitialiseLocusCorrelation(rho_[0]);
  else{
    f[0] = f[1] = 0.0;
    for(unsigned int j = 1; j < NumberOfCompositeLoci; j++ ){
      double rho = rho_[j + _startLocus];
      if(isX)rho *= 0.5;
      f[2*j] = f[2*j + 1] = ( -GetDistance( j ) * rho > -700) ? exp( -GetDistance( j ) * rho ) : 0.0;
    }
  }
}

//sets f for global rho, call after global rho is updated
void Chromosome::SetLocusCorrelation(const double rho_){
  double rho = rho_;
  if(isX)rho *= 0.5;
  for(unsigned j = 1; j < NumberOfCompositeLoci; j++ ){
    f[2*j] = f[2*j + 1] = exp( -GetDistance( j ) * rho );
  }
}
void Chromosome::SetLocusCorrelation(const vector<double> rho_, bool global, bool RandomMating=false){
  //rho should be either
  //(1) vector of length 1 containing a global (global rho model) or 
  //    individual-level (non globalrho and assortative mating) sumintensities
  //(2) vector of length 2 containing gamete-specific sumintensities for a single individual
  //(3) vector of length #loci, containing locus-specific sumintensities
  if(global){//global or individual-specific
    SetLocusCorrelation(rho_[0]);
  }
  else{
    if(RandomMating){//gamete-specific
      if(rho_.size()!=2)throw string("Bad arguments passed to Chromosome::SetLocusCorr");
      for( unsigned int jj = 1; jj < NumberOfCompositeLoci; jj++ ){
	f[2*jj] = exp( -GetDistance( jj ) * rho_[0] );
	f[2*jj + 1] = exp( -GetDistance( jj ) * rho_[1] );
	//TODO: ?? case of X chromosome: use 0.5*rho
      }
    }
    else{//locus-specific
      if(rho_.size()<NumberOfCompositeLoci)throw string("Bad arguments passed to Chromosome::SetLocusCorr");
      for(unsigned int j = 1; j < NumberOfCompositeLoci; j++ ){
	double rho = rho_[j+_startLocus];
	if(isX)rho *= 0.5;
	f[2*j] = f[2*j + 1] = exp( -GetDistance( j ) * rho );
      }
    }
  }
}

// ********** Interface to HMM ****************************************
void Chromosome::SetGenotypeProbs(double* const GenotypeProbs, bool* const GenotypesMissing) {
  SampleStates.SetGenotypeProbs(GenotypeProbs, GenotypesMissing);
}
void Chromosome::SetStateArrivalProbs(const double* const Admixture, bool RandomMating, bool diploid) {

  Diploid = diploid;//required for sampling of locus ancestry
  //construct StateArrivalProbs
  SampleStates.SetStateArrivalProbs(f, Admixture, RandomMating, Diploid);
}

void Chromosome::SampleLocusAncestry(int *OrderedStates){
  SampleStates.Sample(OrderedStates, Diploid);
}

std::vector<std::vector<double> > Chromosome::getAncestryProbs(const bool isDiploid, int j){
  //sets conditional probabilities of ancestry at locus j
  //One row per population, Cols 0,1,2 are probs that 0,1,2 of the 2 gametes have ancestry from that population
  //i.e. (i,2) = p_{ii}
  //     (i,1) = \sum_j{p_{ij}} +   \sum_j{p_{ji}} - 2.0*p_{ii}
  //     (i,0) = 1.0 - (i,1) - (i,2)
  //where p's are probs in StateProbs from HMM

  return SampleStates.Get3WayStateProbs(isDiploid, j);
}

//accessor for HMM Likelihood
double Chromosome::getLogLikelihood(const bool isDiploid)
{
  return SampleStates.getLogLikelihood(isDiploid);
}

//samples jump indicators xi for this chromosome, 
//updates SumLocusAncestry
void Chromosome::SampleJumpIndicators(const int* const LocusAncestry, const unsigned int gametes, 
				      int *SumLocusAncestry, vector<unsigned> &SumN, bool SampleArrivals)const {
  SampleStates.SampleJumpIndicators(LocusAncestry, gametes, SumLocusAncestry, SumN, SampleArrivals, _startLocus);
}


