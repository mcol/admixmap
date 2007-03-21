/** 
 *   ADMIXMAP
 *   Chromosome.cc 
 *   Class to perform chromosome-wise updates and implement HMM class.
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
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
#include "utils/misc.h"

#define PR(x) cerr << #x << " = " << x << endl;

using namespace std;

Chromosome::Chromosome(){
    Distances = 0;
    NumberOfCompositeLoci = 0;
    NumHiddenStates = 0;
    isX = false;
    //Diploid = true;
    f = 0;
    CodedStates = 0;
}

Chromosome::Chromosome(int n, int size, int start, int inNumHiddenStates, bool isx = false) 
		      //size = number of comp loci on chromosome
{
  Number = n;
  _startLocus = start;
  NumHiddenStates = inNumHiddenStates;
  isX = isx;
  if ( size < 1 ){
    size = 1;
  }
  NumberOfCompositeLoci = size;
  Distances = new double[ NumberOfCompositeLoci ];

  CodedStates = new int[size];
  f = new double[2*size];
  f[0] = f[1] = 0.0;

  SampleStates.SetDimensions( size, NumHiddenStates, f);
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
///sets label for this chromosome, usually "1", "2", etc or "X"
void Chromosome::SetLabel( std::string label )
{
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
///returns distance between a locus and the previous locus
double Chromosome::GetDistance( int locus )const
{
  return( Distances[ locus ] );
}

/// Returns the number of the num'th compositelocus on this chromosome
/// eg if chromosome 2 consists of compositeloci 5, 6, 7 and 8,
/// GetLocus(i) returns 5 + i
int Chromosome::GetLocus(int num)const{
  return _startLocus + num;
}
//??: remove either of these functions
///returns number of composite loci in the chromosome
unsigned int Chromosome::GetSize()const{
  return NumberOfCompositeLoci;
}
unsigned int Chromosome::GetNumberOfCompositeLoci()const
{
  return NumberOfCompositeLoci;
}

// ****************** Setting of locus correlation, f *************************

///Sets locus ancestry correlations f for locus-specific lambda (= rho*distance).
void Chromosome::SetLocusCorrelation(const vector<double>::const_iterator lambda_iter){
  for(unsigned int j = 1; j < NumberOfCompositeLoci; j++ ){
    double lambda = *(lambda_iter + j -1);//rho_[j +_startLocus];
    if(isX)lambda *= 0.5;
    f[2*j] = f[2*j + 1] = myexp( - lambda );
  }
}

///sets f for global rho, call after global rho is updated
void Chromosome::SetLocusCorrelation(const double rho_){
  double rho = rho_;
  if(isX)rho *= 0.5;
  for(unsigned j = 1; j < NumberOfCompositeLoci; j++ ){
    f[2*j] = f[2*j + 1] = myexp( -GetDistance( j ) * rho );
  }
}
/**
   Sets locus correlation, f.
   rho should be either
   (1) vector of length 1 containing a global (global rho model) or 
   individual-level (non globalrho and assortative mating) sumintensities
   (2) vector of length 2 containing gamete-specific sumintensities for a single individual
   (3) vector of length numloci, containing locus-specific sumintensities
*/
void Chromosome::SetLocusCorrelation(const std::vector<double> rho_, bool global, bool RandomMating=false){
  if(global){//global or individual-specific
    SetLocusCorrelation(rho_[0]);
  }
  else{
    if(RandomMating){//gamete-specific
      if(rho_.size()!=2)throw string("Bad arguments passed to Chromosome::SetLocusCorr");
      for( unsigned int j = 1; j < NumberOfCompositeLoci; j++ ){
	if(isX){
	  f[2*j] = myexp( -GetDistance( j ) * rho_[0]*0.5 );//sumintensities on Xchrm is set to half autosomal value
	  f[2*j + 1] = myexp( -GetDistance( j ) * rho_[1]*0.5 );
	}
	else{
	  f[2*j] = myexp( -GetDistance( j ) * rho_[0] );
	  f[2*j + 1] = myexp( -GetDistance( j ) * rho_[1] );
	}

      }
    }
    //else{//locus-specific
    //SetLocusCorrelation(rho_);
    //}
  }
}

// ********** Interface to HMM ****************************************
void Chromosome::SetGenotypeProbs(const GenotypeProbIterator& GenotypeProbs, const bool* const GenotypesMissing) {
  SampleStates.SetGenotypeProbs(GenotypeProbs, GenotypesMissing);
}
void Chromosome::SetHMMTheta(const MixturePropsWrapper& Admixture, const MixturePropsWrapper& ThetaSq,
			     const MixturePropsWrapper& ThetaSqInv){
    SampleStates.SetTheta(Admixture, ThetaSq, ThetaSqInv);
}

///sets state arrival probs in HMM, only required for diploid case
void Chromosome::SetStateArrivalProbs(bool RandomMating, bool isdiploid) {
  SampleStates.SetStateArrivalProbs(RandomMating, isdiploid);
}
///samples locus ancestry (hidden states in HMM)
void Chromosome::SampleLocusAncestry(int *OrderedStates, bool diploid){
  SampleStates.Sample(OrderedStates, diploid);
}

/**
   sets conditional probabilities of ancestry at locus j.
   One row per population, Cols 0,1,2 are probs that 0,1,2 of the 2 gametes have ancestry from that population
   i.e. (i,2) = p_{ii}

        (i,1) = \sum_j{p_{ij}} +   \sum_j{p_{ji}} - 2.0*p_{ii}

	(i,0) = 1.0 - (i,1) - (i,2)
	where p's are probs in StateProbs from HMM
*/
std::vector<std::vector<double> > Chromosome::getAncestryProbs(const bool isDiploid, int j){
  //return SampleStates.Get3WayStateProbs(isDiploid, j);
  const std::vector<double> probs = SampleStates.GetHiddenStateProbs(isDiploid, j);
  std::vector<std::vector<double> >AncestryProbs(3);

  if(isDiploid){
    for( int k1 = 0; k1 < NumHiddenStates; k1++ ){
      AncestryProbs[2].push_back(probs[ ( NumHiddenStates + 1 ) * k1 ]);//prob of 2 copies = diagonal element 
      AncestryProbs[1].push_back( 0.0 );
      for( int k2 = 0 ; k2 < NumHiddenStates; k2++ )
	AncestryProbs[1][k1] += probs[k1*NumHiddenStates +k2] + probs[k2*NumHiddenStates +k1];
      AncestryProbs[1][k1] -= 2.0*AncestryProbs[2][k1];
      AncestryProbs[0].push_back( 1.0 - AncestryProbs[1][k1] - AncestryProbs[2][k1] );
    }
  }
  else{//haploid case
    for( int k = 0; k < NumHiddenStates; k++ ){
      AncestryProbs[2].push_back( 0.0 );//cannot have two copies if haploid
      AncestryProbs[1].push_back(probs[k]);//probability of one copy
      AncestryProbs[0].push_back(1.0-probs[k]);//probability of no copies
    }
  }
  return AncestryProbs;
}

///returns HMM Likelihood
double Chromosome::getLogLikelihood(const bool isDiploid)
{
    return SampleStates.getLogLikelihood(isDiploid);
}

///samples jump indicators xi for this chromosome and updates SumLocusAncestry
void Chromosome::SampleJumpIndicators(const int* const LocusAncestry, const unsigned int gametes, 
				      int *SumLocusAncestry, std::vector<unsigned> &SumN, bool SampleArrivals)const {
  SampleStates.SampleJumpIndicators(LocusAncestry, gametes, SumLocusAncestry, SumN, SampleArrivals, _startLocus);
}
void Chromosome::SampleJumpIndicators(const int* const HiddenStates, const unsigned int gametes, 
				      int *SumHiddenStates)const {
  SampleStates.SampleJumpIndicators(HiddenStates, gametes, SumHiddenStates);
}

const pvector<double>& Chromosome::getHiddenStateProbs(const bool isDiploid, int t)
{
  return SampleStates.GetHiddenStateProbs(isDiploid, t);
}
