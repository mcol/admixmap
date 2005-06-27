/** 
 *   ADMIXMAP
 *   Chromosome.cc 
 *   Class to perform chrosome-wise updates and implement HMM class.
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
#include "Chromosome.h"
#include "Individual.h"
#include <iostream>
#include "functions.h"

#define PR(x) cerr << #x << " = " << x << endl;

using namespace std;

Chromosome::Chromosome(int size, int start, int inpopulations) : Genome(size)
{
  _startLocus = start;
  populations = inpopulations;
  D = populations * populations;

  SampleStates.SetDimensions( size, populations, true );

  Lambda = new double[size* populations * populations];

  CodedStates = new int[size];
  for(int j = 0; j < 2; ++j) f[j] = new double[size];
}

void Chromosome::ResetStuffForX()
{
  D = populations;
  SampleStates.SetDimensions( NumberOfCompositeLoci, populations, false );
}

void Chromosome::SetLabel( string label )
{
   _Label = label;
}

string Chromosome::GetLabel( int )
{
   return _Label;
}

Chromosome::~Chromosome()
{
  //TODO: delete these properly
  delete CodedStates;
  delete[] Lambda;
  delete[] f[0];
  delete[] f[1];
}


// Returns the number of the num'th compositelocus on this chromosome
// eg if chromosome 2 consists of compositeloci 5,6,7 and 8,
// GetLocus(i) returns 5 + i
int
Chromosome::GetLocus(int num){
  return _startLocus + num;
}

//returns number of composite loci in the chromosome
unsigned int Chromosome::GetSize(){
  return NumberOfCompositeLoci;
}

//Initialises f for global rho. Necessary since individual-level parameters updated before global rho (in Latent)
void Chromosome::InitialiseLociCorr(const double rho){
  for(unsigned int j = 1; j < NumberOfCompositeLoci; j++ )
    f[0][j] = f[1][j] = ( -GetDistance( j ) * rho > -700) ? exp( -GetDistance( j ) * rho ) : 0.0;
}

//sets f for global rho, called after Latent is updated
void Chromosome::SetLociCorr(const double rho){
    for(unsigned int jj = 1; jj < NumberOfCompositeLoci; jj++ ){
      f[0][jj] = f[1][jj] = exp( -GetDistance( jj ) * rho );
    }
}

void Chromosome::UpdateParameters(Individual* ind, double *Admixture, AdmixOptions* options,  
				  std::vector< double > _rho, bool chibindicator, bool diploid, bool randomAlleleFreqs ){
  //chibindicator is only to facilitate the Chib algorithm in Individual; instructs CompositeLocus to use HapPairProbsMAP
  //instead of HapPairProbs when allelefreqs are not fixed.
  //randomAlleleFreqs might not be necessary if CompositeLocus knows if AlleleFreqs are fixed or random

  // f0 and f1 are arrays of scalars of the form exp - rho*x, where x is distance between loci
  // required to calculate transition matrices 
  if( options->getRhoIndicator() ){

    for( unsigned int jj = 1; jj < NumberOfCompositeLoci; jj++ ){
      f[0][jj] = exp( -GetDistance( jj ) * _rho[0] );
      if( options->isRandomMatingModel() ){
	  f[1][jj] = exp( -GetDistance( jj ) * _rho[1] );
      }
      else
	f[1][jj] = f[0][jj];
    }
  }
  //global rho case already dealt with

  int locus = _startLocus;
  bool test = (options->getTestForAffectedsOnly() || options->getTestForLinkageWithAncestry());
  
  //construct Lambda
  for(unsigned int j = 0; j < NumberOfCompositeLoci; j++ ){
    if( !(ind->IsMissing(locus)) ){
      TheArray[j]->GetGenotypeProbs(Lambda+j*populations*populations, ind->getPossibleHapPairs(locus), chibindicator, randomAlleleFreqs);
      }
    else{
      for( int k = 0; k < populations*populations;k++) Lambda[j*populations*populations + k] = 1.0;
    }
    locus++;
  }

  if(diploid){
    //construct StateArrivalProbs
    SampleStates.SetStateArrivalProbs(f, Admixture, options->isRandomMatingModel());
 
    //Update Forward/Backward Probs in HMM
    SampleStates.UpdateForwardProbsDiploid(f, Lambda);
    if(test){
      SampleStates.UpdateBackwardProbsDiploid(f, Lambda);
    }
 
  }

  else{//haploid
    SampleStates.UpdateProbsHaploid(f, Admixture, Lambda, test);
  }

}

void Chromosome::SampleLocusAncestry(int *OrderedStates, double *Admixture, bool isdiploid){

  SampleStates.Sample(OrderedStates, Admixture, f, isdiploid);
}

void Chromosome::getAncestryProbs( int j, double AncestryProbs[][3] ){
  //
  //sets conditional probabilities of ancestry at locus j
  //One row per population, Cols 0,1,2 are probs that 0,1,2 of the 2 gametes have ancestry from that population
  //i.e. (i,2) = p_{ii}
  //     (i,1) = \sum_j{p_{ij}} +   \sum_j{p_{ji}} - 2.0*p_{ii}
  //     (i,0) = 1.0 - (i,1) - (i,2)
  //where p's are probs in StateProbs from HMM

  SampleStates.Get3WayStateProbs(j, AncestryProbs);
}

//accessor for HMM Likelihood
double Chromosome::getLogLikelihood()
{
  return SampleStates.getLikelihood();
}

//samples jump indicators xi for this chromosome, 
//updates sumxi and Sumrho0 and SumLocusAncestry
void Chromosome::SampleJumpIndicators(int *LocusAncestry, const unsigned int gametes, 
				      int *sumxi, double *Sumrho0, Matrix_i *SumLocusAncestry, Matrix_i *SumLocusAncestry_X, bool isX, 
				      unsigned int SumN[], unsigned int SumN_X[], bool RhoIndicator){

  SampleStates.SampleJumpIndicators(LocusAncestry, f, gametes, Distances, _startLocus, 
				    sumxi, Sumrho0, SumLocusAncestry, SumLocusAncestry_X, isX, 
				    SumN, SumN_X, RhoIndicator);
}


