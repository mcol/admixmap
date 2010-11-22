//=============================================================================
//
// Copyright (C) 2002-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file HiddenMarkovModel.cc
/// Implementation of the HiddenMarkovModel class.
//=============================================================================

#include "HiddenMarkovModel.h"
#include <cmath>
#include "bclib/rand.h"
#include "bclib/exceptions.h"


using namespace std;
using bclib::Rand;

///set pointers to null etc
void HiddenMarkovModel::SetNullValues(){
  alpha = 0;
  beta = 0;
  LambdaBeta = 0;
  p = 0;
  StateArrivalProbs[0] = 0;
  StateArrivalProbs[1] = 0;
  pi = 0;
  piInv = 0;
  colProb = 0;
  Expectation0 = 0;
  Expectation1 = 0;
  rowProb = 0;
  cov = 0;
  f = 0;
  Lambda = 0;
  alphaIsBad = true;
  betaIsBad = true;
  sumfactor = 0.0;
}



//-----------------------------------------------------------------------------
//
/// Constructor with arguments.
///
/// This version of the HMM only works for individuals, not pedigrees, so the
/// assumption is built-in that the number of states (@a _nStates) must be the
/// square of the number of populations (@a _K).  We nonetheless pass both in as
/// separate arguments so that we can have a standardized interface to HMMBase
/// that can be used by both "specialized" (individual) and "generalized"
/// (pedigree) versions.
//
//-----------------------------------------------------------------------------

HiddenMarkovModel::HiddenMarkovModel( int _transitions, size_t _K, size_t _nStates, const double * _f) :
    K( _K ),
    nStates( _nStates ),
    Transitions( _transitions )
  {
  gp_assert_eq( _nStates, _K*_K );

  SetNullValues();
  SetDimensions(_f);
  }


HiddenMarkovModel::~HiddenMarkovModel()
{
  delete[] p;
  delete[] alpha;
  delete[] beta;
  delete[] LambdaBeta;
  delete[] StateArrivalProbs[0];
  delete[] StateArrivalProbs[1];
  delete[] pi;
  delete[] piInv;
  delete[] rowProb;
  delete[] colProb;
  delete[] Expectation0;
  delete[] Expectation1;
  delete[] cov; //free_matrix(cov, K);
}

/// Allocate the arrays and set the pointer of locus correlations f
void HiddenMarkovModel::SetDimensions(const double* const fin) {

  alpha = new double[Transitions* nStates];
  p = new double[Transitions];
  f = fin;
  StateArrivalProbs[0] = new double[Transitions * K ];
  StateArrivalProbs[1] = new double[Transitions * K ];
  pi = new double[ nStates ];
  piInv = new double[ nStates ];

  SetArraysForRecursionProbs();
}

void HiddenMarkovModel::SetArraysForRecursionProbs(){
  if(K>2){
    rowProb = new double[K];
    colProb = new double[K];
    Expectation0 = new double[K];
    Expectation1 = new double[K];
    cov = new double[ nStates ];
  }
}

///set the pointers to genotype probs and missingGenotypes
void HiddenMarkovModel::SetGenotypeProbs(const double* const lambdain, const bool* const missing){
  Lambda = lambdain;
  missingGenotypes = missing;
  alphaIsBad = true;//new input so reset
  betaIsBad = true;
}

/**
   set state arrival probs.
   Set arrival rates and set pointer to ?? probs and , if diploid, calculate pi and piInv. 
   Requires f to have been set.
   \param Theta mixture proportions
   \param Mcol Maternal gamete column (0 if assortative mating, 1 if random mating)
   \param isdiploid indicator for diploidy
*/
void HiddenMarkovModel::SetStateArrivalProbs(const double* const Theta, const int Mcol, const bool isDiploid){
  theta = Theta;
  alphaIsBad = true;//new input so reset
  betaIsBad = true;

  if(isDiploid){
    for(int j0 = 0; j0 < K; ++j0) {
      for(int j1 = 0; j1 < K; ++j1) {
	pi[j0*K + j1] = Theta[j0] * Theta[j1 + K*Mcol];
	piInv[j0*K + j1] = 1.0 / pi[j0*K + j1];
      }
    }
  }

  //required only for diploid updates
  for(int t = 1; t < Transitions; t++ ){        
    for(int j = 0; j < K; ++j){
      StateArrivalProbs[0][t*K + j]    = (1.0 - f[2*t]) * theta[j];
      if(isDiploid)
	StateArrivalProbs[1][t*K + j] = (1.0 - f[2*t + 1]) * theta[K*Mcol +j ];
    }
    p[t] = f[2*t] * f[2*t + 1];
  }
}

void HiddenMarkovModel::SampleJumpIndicators(const int* const LocusAncestry, const unsigned int gametes, 
			       int *SumLocusAncestry, vector<unsigned> &SumNumArrivals, bool SampleArrivals, unsigned startlocus)const {
  //this does not require forward or backward probs, just state arrival probs
  if(!Lambda || !theta || !f)throw string("Error: Call to HiddenMarkovModel::SampleJumpIndicators when StateArrivalProbs are not set!");
  bool xi = false;//jump indicator
  double ProbJump = 0.0; // prob jump indicator is 1
  // first locus not included in loop below
  for( unsigned int g = 0; g < gametes; g++ ){
    SumLocusAncestry[ g*K + LocusAncestry[g*Transitions] ]++;
  }
  for( int t = 1; t < Transitions; t++ ) {
    for( unsigned int g = 0; g < gametes; g++ ){
      xi = true;
      if( LocusAncestry[g*Transitions + t-1] == LocusAncestry[g*Transitions + t] ){
	ProbJump = StateArrivalProbs[g][t*K + LocusAncestry[t + g*Transitions]];  
	xi = (bool)(ProbJump / (ProbJump + f[2*t+g]) > Rand::myrand());
      } 
      if( xi ){ // increment sumlocusancestry if jump indicator is 1
	SumLocusAncestry[ g*K + LocusAncestry[t+g*Transitions] ]++;
	if(SampleArrivals) { // sample number of arrivals where jump indicator is 1
	  double u = Rand::myrand();
	  // sample distance dlast back to last arrival, as dlast = -log[1 - u(1-f)] / rho
	  // then sample number of arrivals before last as Poisson( rho*(d - dlast) )
	  // algorithm does not require rho or d, only u and f
	  unsigned int sample = Rand::genpoi( log( (1 - u*( 1 - f[2*t+g])) / f[2*t+g] ) );
	  SumNumArrivals[2*(startlocus + t) + g] += sample + 1;
	}
      }//end if xi true
    }//end gamete loop
  } // ends loop over intervals
}

/**
  returns log-likelihood.
  This is the sum over states of products of alpha and beta
  and is the same for all t so it is convenient to compute for
  t = T-1 since beta here is 1 so no backward recursions are required. 
*/
double HiddenMarkovModel::getLogLikelihood(const bool isDiploid) { 
  UpdateForwardProbs(isDiploid);

  double sum = 0.0;
  const int NumStates = isDiploid? nStates : K;
  const double *alpha_Tm1 = alpha + (Transitions - 1) * NumStates;
  
  for (int j = 0; j < NumStates; ++j)
    sum += alpha_Tm1[j];

  return( sumfactor+log(sum) );
}

void HiddenMarkovModel::SampleHiddenStates(int *SStates, const bool isDiploid){
  UpdateForwardProbs(isDiploid);

  int j1,j2;
  if(isDiploid) { 
    double* V = new double[nStates]; //probability vector for possible states (haploid or diploid)
    int C = 0; // sampled state (haploid or diploid) coded as integer
    // array Sstates: elements 0 to T-1 represent paternal gamete, elements T to 2T-1 represent maternal gamete 
    int State = 0;
    // sample rightmost locus  
    for( int j = 0; j < nStates; ++j ) V[State++] = alpha[(Transitions - 1)*nStates + j];
    
    C = Rand::SampleFromDiscrete( V, nStates ); 
    SStates[Transitions-1] = (int)(C/K);
    SStates[Transitions - 1 + Transitions] = (C % K);
    
    for( int t =  Transitions - 2; t >= 0; t-- ) { // loop from right to left
      j1 = SStates[t+1];                     // ancestry on gamete 0 at locus t+1
      j2 = SStates[Transitions + t + 1];     // ancestry on gamete 1 at locus t+1
      State = 0;
      for(int i1 = 0; i1 < K; ++i1)for(int i2 = 0; i2 < K; ++i2) {
	V[State] = 
	  ( (i1==j1)*f[2*t+2] + StateArrivalProbs[0][(t+1)*K + j1] ) * 
	  ( (i2==j2)*f[2*t+3] + StateArrivalProbs[1][(t+1)*K + j2] );
	V[State] *= alpha[t*nStates + i1*K + i2];
	++State;
      }
      C = Rand::SampleFromDiscrete( V, nStates );
      SStates[t] = (int)(C/K);//paternal
      SStates[t + Transitions] = (C % K);//maternal
    }
    delete[] V;
  } else {//haploid
    double* V = new double[K]; //probability vector for possible states 
    for( int state = 0; state < K; state++ )V[state] = alpha[(Transitions - 1)*K + state ];
    SStates[Transitions-1] = Rand::SampleFromDiscrete( V, K );
    for( int t =  Transitions - 2; t >= 0; t-- ){
	//for(int j = 0; j < K; j++)V[j] = (j == C[t+1])*f[2*t+1]+theta[C[t+1]]*(1.0 - f[2*t]);
	for(int state = 0; state < K; state++)
	    V[state] = (state == SStates[t+1]) * f[2*t+2] + StateArrivalProbs[0][(t+1)*K + SStates[t+1]  ];
	for( int state = 0; state < K; state++ ) 
	    V[state] *= alpha[t*K + state];
	SStates[t] = Rand::SampleFromDiscrete( V, K );
    }
    delete[] V;
  }
}

/** 
    Returns a vector of conditional probabilities of each hidden state
    at 'time' t.
    If diploid, the vector is really a matrix of probabilities of pairs
    of states.
 */
const bclib::pvector<double>& HiddenMarkovModel::GetHiddenStateProbs(bool isDiploid, int t){
  UpdateForwardProbs(isDiploid);
  UpdateBackwardProbs(isDiploid);
  unsigned States = K;
  if(isDiploid) {
    States=nStates;
  }
  
  if (hiddenStateProbs.size() != (unsigned)States) {
    hiddenStateProbs.resize(States);
  }

  const double *a = alpha + t*States;
  const double *b = beta + t*States;  
  for( unsigned j = 0; j < States; j++ ){
    //    hiddenStateProbs[j] = alpha[t*States + j] * beta[t*States + j];
    hiddenStateProbs[j] = (*(a++)) * ( *(b++) );
  }

  hiddenStateProbs.normalize();
  return hiddenStateProbs;
}

// ****** End Public Interface *******

/**
 * Update forward probabilities if needed.
 * \param isDiploid indicator for diploidy
 */
void HiddenMarkovModel::UpdateForwardProbs(bool isDiploid)
{
  if(alphaIsBad){
    if(isDiploid) {
      UpdateForwardProbsDiploid();
    } else {
      UpdateForwardProbsHaploid();
    }
  }
}

/**
 * Update backward probabilities if needed.
 * 
 * \param isDiploid indicator for diploidy
 */
void HiddenMarkovModel::UpdateBackwardProbs(bool isDiploid)
{
  if(betaIsBad){
    if(isDiploid) {
      UpdateBackwardProbsDiploid();
    } else {
      UpdateBackwardProbsHaploid();
    }
  }
}

/// Updates Forward probabilities alpha, diploid case only
void HiddenMarkovModel::UpdateForwardProbsDiploid(){
  if(!Lambda || !theta || !f)throw string("Error: Call to HiddenMarkovModel when inputs are not set!");
  // if genotypes missing at locus, skip multiplication by lambda and scaling at next locus   
  sumfactor = 0.0; // accumulates log-likelihood
  
  if(!missingGenotypes[0]) {
    for(int j = 0; j < nStates; ++j) {
      alpha[j] =  pi[j] * Lambda[j];
    } 
  } else {
    for(int j = 0; j < nStates; ++j) {
      alpha[j] = pi[j];
    }
  }

  for( int t = 1; t < Transitions; ++t ){
    if(!missingGenotypes[t-1]) {
      double Sum = 0.0;
      //scale previous alpha to sum to 1
      for( int j = 0; j <  nStates; ++j ) {
	Sum += alpha[(t-1)*nStates +j];
      }
      const double scaleFactor = 1.0 / Sum;
      for( int j = 0; j <  nStates; ++j ) {
	alpha[(t-1)*nStates +j] *= scaleFactor;
      }
      sumfactor += log(Sum);
    }
    
    RecursionProbs(p[t], f + 2*t, StateArrivalProbs[0] + t*K, StateArrivalProbs[1] + t*K, alpha + (t-1)*nStates, alpha + t*nStates);
    for(int j = 0; j < nStates; ++j){
      if(!missingGenotypes[t]) {
	alpha[t*nStates +j] *=  Lambda[t*nStates + j]; //*lam++; 
	//++lam;
      } //else ++lam;
    }

  }
  alphaIsBad = false;
}

/// Updates backard probabilities beta, diploid case only
void HiddenMarkovModel::UpdateBackwardProbsDiploid(){
  if(!Lambda || !theta || !f)throw string("Error: Call to HiddenMarkovModel when inputs are not set!");
  if(!beta) { // allocate beta array if not already done
    beta =  new double[Transitions * nStates];
  }
  if(!LambdaBeta)
    LambdaBeta = new double[ nStates ];

  for(int j = 0; j < nStates; ++j){
    //set beta(T) = 1
    beta[(Transitions - 1)*nStates + j] = 1.0;
  }
  
  for( int t = Transitions-2; t >= 0; --t ) {
    double f2[2] = {f[2*t + 2], f[2*t + 3]};
    double Sum = 0.0;
    for(int j = 0; j < nStates; ++j){
      LambdaBeta[j] = Lambda[(t+1)*nStates + j] * beta[(t+1)*nStates + j] * pi[j];
      Sum += LambdaBeta[j];
    }
    //scale LambdaBeta to sum to 1
    const double scaleFactor = 1.0 / Sum;
    for( int j = 0; j <  nStates; ++j ) {
      LambdaBeta[j] *= scaleFactor;
    }
    
    RecursionProbs(p[t+1], f2, StateArrivalProbs[0]+(t+1)*K, StateArrivalProbs[1]+(t+1)*K, LambdaBeta, beta+ t*nStates);
    for(int j = 0; j < nStates; ++j){ // vectorization successful
      beta[t*nStates + j] *= piInv[j];
    }
  }
  betaIsBad = false;
}

/**
  Updates forward probs, haploid case only.
  Here theta is a column matrix and the last dimensions of f and lambda are 1.
*/
void HiddenMarkovModel::UpdateForwardProbsHaploid(){
  if(!Lambda || !theta || !f)throw string("Error: Call to HiddenMarkovModel when inputs are not set!");
  sumfactor = 0.0;

  //   for(int t = 0; t < Transitions; ++t) {
  //     cout << "t " << t << " ";
  //     for(int j = 0; j < K; ++j) {
  //       cout << Lambda[t*K + j] << " ";
  //     }
  //     cout << "\n";
  //   }
  
  for(int j = 0; j < K; ++j){
    alpha[j] = theta[j] * Lambda[j];
  }
  
  for( int t = 1; t < Transitions; t++ ) {
    if(!missingGenotypes[t-1]) {
      double Sum = 0.0;
      //scale previous alpha to sum to 1
      for(int j = 0; j < K; ++j){
	Sum += alpha[(t-1)*K + j];
      }
      const double scaleFactor = 1.0 / Sum;
      for( int j = 0; j <  K; ++j ) {
	alpha[(t-1)*K +j] *= scaleFactor;
      }
      sumfactor += log(Sum);
    }

    for(int j = 0; j < K; ++j){
      alpha[t*K + j] = f[2*t]*alpha[(t-1)*K +j] + (1.0 - f[2*t]) * theta[j];
      alpha[t*K + j] *= Lambda[t*K + j];
    }
  }
  alphaIsBad =  false;
}

///updates forward probs, haploid case only
void HiddenMarkovModel::UpdateBackwardProbsHaploid(){
  if(!Lambda || !theta || !f)throw string("Error: Call to HiddenMarkovModel when inputs are not set!");
  if(!beta) { // allocate diploid-sized beta array if not already done
    beta =  new double[Transitions*nStates];
  }
  if(!LambdaBeta)
    LambdaBeta = new double[ nStates ];

  for(int j = 0; j < K; ++j){
    beta[(Transitions-1)*K + j] = 1.0;
  }
  
  for( int t = Transitions-2; t >=0; t-- ){
    double Sum = 0.0;
    for(int j = 0; j < K; ++j){
      LambdaBeta[j] = Lambda[(t+1)*K + j]*beta[(t+1)*K + j];
      Sum += theta[j]*LambdaBeta[j];
    }
    for(int j = 0; j < K; ++j){
      beta[t*K + j] = f[2*(t+1)]*LambdaBeta[j] + (1.0 - f[2*(t+1)+1])*Sum;
    }
  }
  betaIsBad = false;
}

/// @param f2        Locus correlations
/// @param oldProbs  Square array of size K * K:
///                  For forward recursions, pass alpha_t and multiply newProbs
///                    by emission probs lambda_t
///                  For backward recursions, pass the array of products
///                    lambda_t+1[jj] * beta_t+1[jj]
void HiddenMarkovModel::RecursionProbs(const double ff, const double f2[2],
                                       const double* const stateArrivalProbs0,
				       const double* const stateArrivalProbs1, 
				       const double* const oldProbs, double *newProbs) {
  if (K==2)
    RecursionProbs2(ff, f2, stateArrivalProbs0, stateArrivalProbs1, oldProbs, newProbs);
  else {
    int j0K = 0;
    for( int j0 = 0; j0 <  K; ++j0 ) {
      rowProb[j0] = 0.0;
      colProb[j0] = 0.0;
      for( int j1 =0; j1 < K; ++j1 ) {
	rowProb[j0] += oldProbs[j0*K + j1];
	colProb[j0] += oldProbs[j1*K + j0];
      }
    }
    // calculate expectations of indicator variables for each ancestry state on each gamete
    for( int j0 = 0; j0 <  K; ++j0 ) {
      Expectation0[j0] = f2[0]*rowProb[j0] + stateArrivalProbs0[j0];
      Expectation1[j0] = f2[1]*colProb[j0] + stateArrivalProbs1[j0];
    }

    // calculate expectation of product as covariance plus product of expectations
    for(int j0 = 0; j0 < K; ++j0) {
      for(int j1 = 0; j1 < K; ++j1) {
        j0K = j0 * K;
	cov[j0K + j1] = ff * ( oldProbs[j0K + j1] - rowProb[j0] * colProb[j1] );
	newProbs[j0K + j1] = cov[j0K + j1] + Expectation0[j0] * Expectation1[j1];
	// newProbs[1] is prob(paternal=1, maternal=0)
      }
    }
  }//end else
}

/// Version of RecursionProbs specialized for 2 populations
void HiddenMarkovModel::RecursionProbs2(const double ff, const double f2[2],
                                        const double* const stateArrivalProbs0,
                                        const double* const stateArrivalProbs1,
                                        const double* const oldProbs,
                                        double *newProbs) {
  // sum row 0 and col 0  
  double row0Prob = oldProbs[0] + oldProbs[2];
  double col0Prob = oldProbs[0] + oldProbs[1];

  // calculate expectations of indicator variables for ancestry=0 on each gamete
  double Exp0 = f2[0]*row0Prob + stateArrivalProbs0[0]; // paternal gamete
  double Exp1 = f2[1]*col0Prob + stateArrivalProbs1[0]; // maternal gamete

  // calculate covariance of indicator variables as ff * deviation from
  // product of row and col probs
  newProbs[0] = Exp0 * Exp1 + ff * ( oldProbs[0] - row0Prob * col0Prob );
  newProbs[1] = Exp0 - newProbs[0]; //prob paternal ancestry=1, maternal=0 
  newProbs[2] = Exp1 - newProbs[0]; //prob paternal ancestry=0, maternal=1 
  newProbs[3] = 1 - Exp0 - Exp1 + newProbs[0];
}
