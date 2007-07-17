/*
 *   HAPMIXMAP
 *   HapMixHMM.hh 
 *   Extension of ADMIXMAP's HMM class with locus-specific Hidden-State probs
 *   Copyright (c) 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "HapMixHMM.hh"
#include "bclib/rand.h"
#include <cmath> //for log

// HapMixHMM::HapMixHMM(){
//   SetNullValues();
// }

HapMixHMM::HapMixHMM( int inTransitions, int NumHiddenStates, const double* const f):
  HiddenMarkovModel(inTransitions, NumHiddenStates, f){
  SetNullValues();
  StateArrivalProbs[0] = 0;
  StateArrivalProbs[1] = 0;
  //  SetDimensions(inTransitions, NumHiddenStates, f);
}

//TODO: allocate beta here if required
//that is, if the allelic assoc score test is on or the PPGenotypes Probs are being written
void HapMixHMM::SetDimensions( int inTransitions, int NumHiddenStates, const double* const fin, bool diploid){
  //inTransitions = #transitions +1 = #Loci 
  p = new double[Transitions];
  f = fin;
  StateArrivalProbs[0] = new double[Transitions * K ];

  if(diploid){
    alpha = new double[Transitions*K*K];
    StateArrivalProbs[1] = new double[Transitions * K];
    ThetaThetaPrime = new double[Transitions * K * K];
    //  ThetaThetaInv = new double[Transitions * K * K];
    
    SetArraysForRecursionProbs(K);
  }
  else{
    alpha = new double[Transitions*K];
  }
}

HapMixHMM::~HapMixHMM(){
}

///set the pointers to genotype probs and missingGenotypes
void HapMixHMM::SetGenotypeProbs(const GenotypeProbIterator& lambdain, const bool* const missing){
  LambdaGPI = lambdain;
  missingGenotypes = missing;
  alphaIsBad = true;//new input so reset
  betaIsBad = true;
}

/**
   Set arrival rates and set pointer to ?? probs. . Requires f to have been set.
   theta is a Transitions * K array.
*/
void HapMixHMM::SetStateArrivalProbs(const double* const Theta, const int , const bool isdiploid){
  theta = Theta;
  alphaIsBad = true;//new input so reset
  betaIsBad = true;

  for(int t = 1; t < Transitions; t++ ){ 
       
    for(int j = 0; j < K; ++j){

      //state arrival probs, gamete 1
      StateArrivalProbs[0][t*K + j]    = (1.0 - f[2*t]) * theta[t*K + j];

      if(isdiploid){
	//state arrival probs, gamete 2
	StateArrivalProbs[1][t*K + j] = (1.0 - f[2*t + 1]) * theta[t*K + j ];
      }

    }//end loop over states
    p[t] = f[2*t] * f[2*t + 1];
  }//end loop over transitions

  //TODO: include this loop in the above. 
  if(isdiploid){
    for(int t = 0; t < Transitions; ++t){
      for(int j0 = 0; j0 < K; ++j0) {
	for(int j1 = 0; j1 < K; ++j1) {
	  ThetaThetaPrime[t*DStates + j0*K + j1] = theta[t*K + j0]*theta[t*K + j1];
	}
      }//end loop over diploid states

    }//end transition loop
  }//end if diploid

}

void HapMixHMM::SampleJumpIndicators(const int* const HiddenStates,  unsigned int gametes,
				     int *SumHiddenStates)const {
  //this does not require forward or backward probs, just state arrival probs
  if(LambdaGPI.isNull() || !theta || !f)
    throw ("Error: Call to HMM::SampleJumpIndicators when StateArrivalProbs are not set!");
  bool xi = false;//jump indicator
  double ProbJump = 0.0; // prob jump indicator is 1
  // first locus not included in loop below
  for( unsigned int g = 0; g < gametes; g++ ){
    SumHiddenStates[ HiddenStates[g*Transitions] ]++;
  }
  for( int t = 1; t < Transitions; t++ ) {
    for( unsigned int g = 0; g < gametes; g++ ){
      xi = true;
      if( HiddenStates[g*Transitions + t-1] == HiddenStates[g*Transitions + t] ){
	ProbJump = StateArrivalProbs[g][t*K + HiddenStates[t + g*Transitions]]; 
	xi = (bool)(ProbJump / (ProbJump + f[2*t+g]) > Rand::myrand());
      }
      if( xi ){ // increment sum if jump indicator is 1
	SumHiddenStates[ t*K + HiddenStates[t+g*Transitions] ]++;
      }//end if xi true
    }//end gamete loop
  } // ends loop over intervals
}

void HapMixHMM::UpdateForwardProbsDiploid(){
  if(LambdaGPI.isNull() || !theta || !f)throw ("Error: Call to HiddenMarkovModel when inputs are not set!");
 
  sumfactor = 0.0; // accumulates log-likelihood
  //sumfactor2 = 0.0;
  double Sum = 0.0;
  double scaleFactor = 0.0;

  for(int j = 0; j < DStates; ++j) {
    alpha[j] =  ThetaThetaPrime[j] * LambdaGPI(0, j);
  } 
  
  for( int t = 1; t < Transitions; ++t ){
    if(!missingGenotypes[t-1]) {
      Sum = 0.0;
      //scale previous alpha to sum to 1
      for( int j = 0; j <  DStates; ++j ) {
	Sum += alpha[(t-1)*DStates +j];
      }
      scaleFactor = 1.0 / Sum;
      for( int j = 0; j <  DStates; ++j ) {
	alpha[(t-1)*DStates +j] *= scaleFactor;
      }
      sumfactor += log(Sum);
      //if(t==1)sumfactor2 += log(Sum);
    }
    
    RecursionProbs(p[t], f + 2*t, StateArrivalProbs[0] + t*K, StateArrivalProbs[1] + t*K, alpha + (t-1)*DStates, alpha + t*DStates);

    if(!missingGenotypes[t]) {    
      for(int j = 0; j < DStates; ++j){
	alpha[t*DStates +j] *=  LambdaGPI(t, j); 
      } 
    }
  }
  alphaIsBad = false;
}

void HapMixHMM::UpdateBackwardProbsDiploid(){
  if(LambdaGPI.isNull() || !theta || !f)throw ("Error: Call to HiddenMarkovModel when inputs are not set!");
  if(!beta) { // allocate beta array if not already done
    beta =  new double[Transitions*K*K];
  }

  //betasumfactor = 0.0;
  
  for(int j = 0; j < DStates; ++j){
    //set beta(T) = 1
    beta[(Transitions - 1)*DStates + j] = 1.0;
  }

  for( int t = Transitions-2; t >= 0; --t ) {
    const double fsq =                f[2*t + 2] *        f[2*t + 3];
    const double f_1minusf =          f[2*t + 2] * (1.0 - f[2*t + 3]);
    const double oneMinusfSq = (1.0 - f[2*t + 2])* (1.0 - f[2*t + 3]);
    std::vector<double> Sumj(K, 0.0);

    //evaluate the linear sums
    for(int j = 0; j < K; ++j)
      for(int i = 0; i < K; ++i)
	Sumj[j] += theta[(t+1)*K+i]*LambdaGPI(t+1, i*K+j)*beta[(t+1)*DStates + i*K + j];

    //evaluate the sum of quadratic terms. 
    double Sumj1j2 = 0.0;
    for(int j = 0; j < K ; ++j){
      Sumj1j2 += theta[(t+1)*K + j]*Sumj[j];
    }

    //now loop over diploid states to evaluate beta
    double NormSum = 0.0;
    for(int i1 = 0; i1 < K ; ++i1)
      for(int i2 = 0; i2 < K; ++i2){
	const int sindex = i1*K +i2;//state index
	beta[t*DStates + sindex] = fsq * LambdaGPI(t+1, sindex)*beta[(t+1)*DStates + sindex]
	  + f_1minusf * ( Sumj[i1] + Sumj[i2])
	  + oneMinusfSq * Sumj1j2;
	//accumulate sum for normalization
	NormSum += beta[t*DStates + sindex];
      }

    //normalize to avoid underflow
    const double Inv = 1.0 / NormSum;//it's better to multiply by the reciprocal
    for(int j = 0; j < DStates; ++j){
      beta[t*DStates + j] *= Inv;
    }

  }

  //beta is now ok to use
  betaIsBad = false;
}

/*
  Updates forward probs, haploid case only
  Here the last dimensions of f and lambda are 1.
*/
void HapMixHMM::UpdateForwardProbsHaploid(){
  if(LambdaGPI.isNull() || !theta || !f)throw ("Error: Call to HiddenMarkovModel when inputs are not set!");
  sumfactor = 0.0;
  //sumfactor2 = 0.0;
  double Sum = 0.0;
  double scaleFactor = 0.0;
  for(int j = 0; j < K; ++j){
    alpha[j] = theta[j] * LambdaGPI(0, j);
  }
  
  for( int t = 1; t < Transitions; t++ ) {
    if(!missingGenotypes[t-1]) {
      Sum = 0.0;
      //scale previous alpha to sum to 1
      for(int j = 0; j < K; ++j){
	Sum += alpha[(t-1)*K + j];
      }
      scaleFactor = 1.0 / Sum;
      for( int j = 0; j <  K; ++j ) {
	alpha[(t-1)*K +j] *= scaleFactor;
      }
      sumfactor += log(Sum);
      //if(t==1) sumfactor2 += log(Sum);
    }

    for(int j = 0; j < K; ++j){
      alpha[t*K + j] = f[2*t]*alpha[(t-1)*K +j] + (1.0 - f[2*t]) * theta[t*K +j];
      alpha[t*K + j] *= LambdaGPI(t, j);
    }
  }
  alphaIsBad =  false;
}

void HapMixHMM::UpdateBackwardProbsHaploid(){
  if(LambdaGPI.isNull() || !theta || !f)throw ("Error: Call to HiddenMarkovModel when inputs are not set!");
  if(!beta) { // allocate diploid-sized beta array if not already done
    //has to be diploid in case there are any diploid individuals
    beta =  new double[Transitions*K*K];
  }
  if(!LambdaBeta)
    LambdaBeta = new double[K*K];

  double Sum = 0.0;
  //betasumfactor = 0.0;

  for(int j = 0; j < K; ++j){
    beta[(Transitions-1)*K + j] = 1.0;
  }
  
  for( int t = Transitions-2; t >=0; t-- ){
    Sum = 0.0;
    for(int j = 0; j < K; ++j){
      LambdaBeta[j] = LambdaGPI(t+1, j)*beta[(t+1)*K + j];
      Sum += theta[(t+1)*K + j]*LambdaBeta[j];
    }

    double NormSum = 0.0;
    for(int j = 0; j < K; ++j){
      beta[t*K + j] = f[2*(t+1)]*LambdaBeta[j] + (1.0 - f[2*(t+1)+1])*Sum;
      NormSum += beta[t*K + j];
    }

    //normalize to avoid underflow
    //it's better to multiply by a reciprocal than divide by a sum
    const double Inv = 1.0 / NormSum;
    for(int j = 0; j < K; ++j){
      beta[t*K + j] *= Inv;
    }
  }
  betaIsBad = false;
}
