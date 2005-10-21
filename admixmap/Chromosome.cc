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

double Chromosome::coolness = 1;//initialising to 1 samples from posterior by default
double *Chromosome::L_mod = 0;

Chromosome::Chromosome(int size, int start, int inpopulations) : Genome(size)
{
  isChromosome = true;
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

const string Chromosome::GetLabel( int )const
{
   return _Label;
}

Chromosome::~Chromosome()
{
  delete[] CodedStates;
  delete[] Lambda;
  delete[] f[0];
  delete[] f[1];
}

void Chromosome::setCoolness(double l, double *Lmod){
  coolness = l;
  L_mod = Lmod;
}

// Returns the number of the num'th compositelocus on this chromosome
// eg if chromosome 2 consists of compositeloci 5,6,7 and 8,
// GetLocus(i) returns 5 + i
int Chromosome::GetLocus(int num)const{
  return _startLocus + num;
}

//returns number of composite loci in the chromosome
unsigned int Chromosome::GetSize()const{
  return NumberOfCompositeLoci;
}

//Initialises f for global rho. Necessary since individual-level parameters updated before global rho (in Latent)
void Chromosome::InitialiseLociCorr(const double rho){
  for(unsigned int j = 1; j < NumberOfCompositeLoci; j++ )
    f[0][j] = f[1][j] = ( -GetDistance( j ) * rho > -700) ? exp( -GetDistance( j ) * rho ) : 0.0;
}

//sets f for global rho, call after global rho is updated
void Chromosome::SetLociCorr(const double rho){
    for(unsigned int jj = 1; jj < NumberOfCompositeLoci; jj++ ){
      f[0][jj] = f[1][jj] = exp( -GetDistance( jj ) * rho );
    }
}

void Chromosome::SetGenotypeProbs(const Individual* const ind, bool chibindicator){
  //chibindicator is only to facilitate the Chib algorithm in Individual; instructs CompositeLocus to use HapPairProbsMAP
  //instead of HapPairProbs when allelefreqs are not fixed.
  int locus = _startLocus;
  for(unsigned int j = 0; j < NumberOfCompositeLoci; j++ ){
    if( !(ind->IsMissing(locus)) ){
      TheArray[j]->GetGenotypeProbs(Lambda+j*populations*populations, ind->getPossibleHapPairs(locus), chibindicator);
      }
    else{
      for( int k = 0; k < populations*populations;k++) Lambda[j*populations*populations + k] = 1.0;
    }
    locus++;
  }
}

void Chromosome::SetGenotypeProbs(double *Probs){
  //Probs is a Populations x Populations array containing genotype probs obtained from CompositeLocus
  copy(Probs, Probs + populations * populations, Lambda);
}

void Chromosome::UpdateHMMForwardProbs(const double* const Admixture, const AdmixOptions* const options, 
				       const std::vector< double > _rho, bool diploid, bool annealindicator = false){
  //set annealindicator to true once per individual per iteration to accumulate unannealed loglikelihood stored in top level

  //_rho contains Individual sumintensities parameters, ignored if globalrho model

  //SetGenotypeProbs(ind, chibindicator);

  // f0 and f1 are arrays of scalars of the form exp - rho*x, where x is distance between loci
  // required to calculate transition matrices 
  if( !options->isGlobalRho() ){//non global rho case

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

  Diploid = diploid;//required for sampling of locus ancestry
  
  if(Diploid){
    //construct StateArrivalProbs
    SampleStates.SetStateArrivalProbs(f, Admixture, options->isRandomMatingModel());

    if(annealindicator) {
      SampleStates.UpdateForwardProbsDiploid(f, Lambda, 1.0);
      (*L_mod) += SampleStates.getLogLikelihood();
    }

    //if(coolness < 1.0)
    //Update Forward/Backward Probs in HMM
    SampleStates.UpdateForwardProbsDiploid(f, Lambda, coolness);

  }
  else{//haploid
    SampleStates.UpdateForwardProbsHaploid(f, Admixture, Lambda);
  }

}

void Chromosome::UpdateHMMBackwardProbs(const double* const hapAdmixture){
  //call only after a call to UpdateHMMForwardProbs
  //this is ok as whenever we need backward probs we also need forward probs but not vice versa
  if(Diploid)  SampleStates.UpdateBackwardProbsDiploid(f, Lambda);
  else SampleStates.UpdateBackwardProbsHaploid(f, hapAdmixture, Lambda);
}

void Chromosome::SampleLocusAncestry(int *OrderedStates, const double* const Admixture)const{

  SampleStates.Sample(OrderedStates, Admixture, f, Diploid);
}

void Chromosome::getAncestryProbs( int j, double AncestryProbs[][3] )const{
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
double Chromosome::getLogLikelihood()const
{
  return SampleStates.getLogLikelihood();
}

//samples jump indicators xi for this chromosome, 
//updates SumLocusAncestry
void Chromosome::SampleJumpIndicators(const int* const LocusAncestry, const unsigned int gametes, 
				      int *SumLocusAncestry, int *SumLocusAncestry_X, bool isX, 
				      unsigned int SumN[], unsigned int SumN_X[], bool isGlobalRho)const{

  SampleStates.SampleJumpIndicators(LocusAncestry, f, gametes,_startLocus, 
				    SumLocusAncestry, SumLocusAncestry_X, isX, 
				    SumN, SumN_X, isGlobalRho);
}


