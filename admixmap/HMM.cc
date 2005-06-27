#include "rand.h"
#include "vector_i.h"
#include "vector_d.h"
#include "matrix_d.h"
#include "HMM.h"
#include "functions.h"

using namespace std;

HMM::HMM()
{
}

//not currently used
//HMM objects are instantiated in Chromosome using default constructor above
//and dimensions set by SetDimensions below
HMM::HMM( int inTransitions, int pops, bool isdiploid )
{
  //inTransitions = #transitions +1 = #Loci 
  //pops = #populations
  K = pops;
  States = isdiploid? K*K : K;//K for haploid, K^2 for diploid

  Transitions = inTransitions;

  alpha = new double[Transitions*States];
  beta =  new double[Transitions*States];
 
  sumfactor=0.0;
  p = new double[Transitions];
  LambdaBeta = new double[States];
}

HMM::~HMM()
{
  //TODO:destroy these properly
  delete[] p;
  delete[] LambdaBeta;
  delete[] alpha;
  delete[] beta;
  delete[] StateArrivalProbs;
  delete[] ThetaThetaPrime;
}

void HMM::SetDimensions( int inTransitions, int pops, bool isdiploid )
{
  //TODO: delete arrays if already allocated
  //this will happen for X chromosome

  //inTransitions = #transitions +1 = #Loci 
  //pops = #populations
  K = pops;
  States = isdiploid? K*K : K;//K for haploid, K^2 for diploid

  Transitions = inTransitions;

  int d = isdiploid? K:1;  
  alpha = new double[Transitions*States];
  beta =  new double[Transitions*States];
  
  sumfactor=0.0;
  p = new double[Transitions];
  LambdaBeta = new double[K*d];

  StateArrivalProbs = new double[Transitions * K * 2];
  ThetaThetaPrime = new double[States];
}

void HMM::SetStateArrivalProbs(double *f[], double *Theta, int Mcol){
  for(int t = 1; t < Transitions; t++ ){        
    for(int j = 0; j < K; ++j){
      StateArrivalProbs[t*K*2 + j*2] = (1.0 - f[0][t]) * Theta[j];
      StateArrivalProbs[t*K*2 + j*2 +1] = (1.0 - f[1][t]) * Theta[K*Mcol +j ];
    }
  }
  for(int j0 = 0; j0 < K; ++j0)for(int j1 = 0; j1 < K; ++j1)
    ThetaThetaPrime[j0*K + j1] = Theta[j0]*Theta[j1 + K*Mcol];
}

/*
  Updates Forward and (if required) Backward probabilities
  diploid case only
  -----------------------------------------------------------
  relates to notation in docs as follows:
  j = j1, j' = j2, i = i1, i' = i2, 
  theta_j = Admixture(j,0),  theta'_j = Admixture(j,Mcol).
*/
void HMM::UpdateForwardProbsDiploid(double *f[], double *lambda)
{
  sumfactor = 0.0;

   for(int j = 0; j < States; ++j)
     //set alpha(0) = StationaryDist * lambda(0) 
     alpha[j] =  ThetaThetaPrime[j] *lambda[j];


  for( int t = 1; t < Transitions; t++ ){        
    p[t] = f[0][t] * f[1][t];
    double f2[2] = {f[0][t], f[1][t]};

    RecursionProbs(p[t], f2, StateArrivalProbs + t*K*2, alpha + (t-1)*States, alpha + t*States);
    for(int j = 0; j < States; ++j){
      alpha[t*States +j] *= lambda[t*States +j];
    }
  }
}

void HMM::UpdateBackwardProbsDiploid(double *f[], double *lambda)
{
  double rec[States];

  for(int j = 0; j < States; ++j){
      //set beta(T) = 1
      beta[(Transitions - 1)*States + j] = 1.0;
      rec[j] = 1.0 / ThetaThetaPrime[j];
    }

  for( int t = Transitions-2; t >=0; t-- ){
    
    double f2[2] = {f[0][t+1], f[1][t+1]};
    
    for(int j = 0; j < States; ++j){
      LambdaBeta[j] = lambda[(t+1)*K*K + j] * beta[(t+1)*States + j] * ThetaThetaPrime[j];
    }
    
    RecursionProbs(p[t+1], f2, StateArrivalProbs+ (t+1)*K*2, LambdaBeta, beta+ t*States);
    for(int j = 0; j < States; ++j){
      beta[t*States + j] *= rec[j];
    }
   }

}


/*
  Updates Forward and (if required) Backward probabilities
  haploid case only
  Here Admixture is a column matrix and the last dimensions of f and lambda are 1.
*/
void HMM::UpdateProbsHaploid(double *f[], double *Admixture, double *lambda, bool CalculateBeta){

  double Sum;

  for(int j=0; j<States;++j){
    alpha[j] = Admixture[j] * lambda[j];
    beta[(Transitions-1)*States + j] = 1.0;
  }

  for( int t = 1; t < Transitions; t++ ){
    Sum = 0.0;
    for(int j=0;j<States;++j){
      Sum += alpha[(t-1)*States + j];
    }
    for(int j=0;j<0;++j){
      alpha[t*States + j] = f[0][t] + (1.0 - f[0][t]) * Admixture[j] * Sum;
      alpha[t*States + j] *= lambda[(t+1)*States + j];
    }

    if(CalculateBeta){
      for( int t = Transitions-2; t >=0; t-- ){
	Sum = 0.0;
	for(int j=0;j<States;++j){
	  Sum += Admixture[j]*lambda[(t+1)*States + j]*beta[(t+1)*States + j];
	}
	for(int j=0;j<States;++j){
	  beta[t*States + j] = f[t+1][0]*lambda[(t+1)*States + j]*beta[(t+1)*States + j] + (1.0 - f[0][t+1])*Sum;
	}
      }
      
    }
  }

}

/*
  computes conditional state probabilities at "time" t
  probs - double array to store state probs
  K = #populations
*/
void HMM::GetStateProbs( double * probs, int t)
{
  double sum = 0.0;
  int State = 0;

  for(int i = 0; i < K; ++i)
   for( int j = 0; j < K; j++ ){
     probs[State++] = alpha[t*States + i*K +j] * beta[t*States + i*K +j];
     sum += probs[State-1];
   }

   for( int j = 0; j < States; j++ ){
     probs[j] /= sum;
   }
}

void HMM::Get3WayStateProbs( int t, double AncestryProbs[][3]){
  double sum = 0.0;
  int State = 0;
  double probs[States];

  for( int j = 0; j < States; j++ ){
    probs[State++] = alpha[t*States + j] * beta[t*States + j];
    sum += probs[State-1];
  }

   for( int j = 0; j < States; j++ ){
     probs[j] /= sum;
   }
   for( int k1 = 0; k1 < K; k1++ ){
     AncestryProbs[k1][2] = probs[ ( K + 1 ) * k1 ];
     AncestryProbs[k1][1] = 0.0;
     for( int k2 = 0 ; k2 < K; k2++ )
       AncestryProbs[k1][1] += probs[k1*K +k2] + probs[k2*K +k1];
     AncestryProbs[k1][1] -= 2.0*AncestryProbs[k1][2];
     AncestryProbs[k1][0] = 1.0 - AncestryProbs[k1][1] - AncestryProbs[k1][2];
   }
}

/*
  returns log-likelihood
  This is the sum over states of products of alpha and beta
  and is the same for all t so it is convenient to compute for
  t = T-1 since beta here is 1. 
*/
double HMM::getLikelihood()
{
   double sum = 0;
   for( int j = 0; j < States; j++ ){
     sum += alpha[(Transitions - 1)*States + j];
   }
   return( sumfactor+log( sum ) );
}

/*
  Samples Hidden States
  ---------------------
  SStates          - an int array to store the sampled states
  isdiploid  - indicator for diploidy
*/
void HMM::Sample(int *SStates, double *Admixture, double *f[], bool isdiploid)
{
  int j1,j2;
  double V[States];
  int C[Transitions];

  if(isdiploid){
    int State = 0;
    for( int j = 0; j < States; j++ )V[State++] = alpha[(Transitions - 1)*States + j];
    C[ Transitions - 1 ] = SampleFromDiscrete3( V, States );
    SStates[Transitions-1] = (int)(C[Transitions-1]/K);
    SStates[Transitions - 1 + Transitions] = (C[Transitions-1] % K);
    
    for( int t =  Transitions - 2; t >= 0; t-- ){
      j1 = (int) (C[t+1]/K);//j
      j2 = C[t+1]-K*j1;     //j'

      State = 0;
      for(int i1 = 0; i1 < K; i1++)for(int i2=0;i2<K;++i2){
	V[State] = 
	  ( (i1==j1)*f[0][t+1] + StateArrivalProbs[(t+1)*K*2 + j1*2] ) * ( (i2==j2)*f[1][t+1] + StateArrivalProbs[(t+1)*K*2 + j2*2 +1] );
	V[State] *= alpha[t*States + i1*K + i2];
	State++;
      }
      C[ t ] = SampleFromDiscrete3( V, States );
      SStates[t] = (int)(C[t]/K);
      SStates[t + Transitions] = (C[t] % K);
    }
   }
  else{//haploid
    for( int j = 0; j < States; j++ )V[j] = alpha[(Transitions - 1)*States + j ];
    C[ Transitions - 1 ] = SampleFromDiscrete3( V, States );
    SStates[Transitions-1] = C[Transitions-1];
    for( int t =  Transitions - 2; t >= 0; t-- ){
      for(int j = 0; j < States; j++)V[j] = (j == C[t+1])*f[0][t+1]+Admixture[C[t+1]]*(1.0 - f[0][t]);
      for( int j = 0; j < States; j++ )	V[j] *= alpha[t*States + j];
      C[ t ] = SampleFromDiscrete3( V, States );
      SStates[t] = C[t];
    }
  }
}

// argument oldProbs is square array of size K, K
// for forward recursions, pass alpha_t and multiply newProbs by observation probs 
// for backward recursions, pass array of products lambda_t+1[jj] * beta_t+1[jj] 
// updates oldProbs (scaled to sum to 1), newProbs and sumfactor if forward = true (for alphas)
void HMM::RecursionProbs(const double ff, const double f[2], 
			 double *stateArrivalProbs, double *oldProbs, double *newProbs) {
  double Sum = 0.0, scaleFactor = 1.0;
  double *rowProb = new double[K];
  double *colProb = new double[K];
  double *Expectation0 = new double[K];
  double *Expectation1 = new double[K];

  double *rowSum = new double[K];
  double *colSum = new double[K];
  double **cov;
  cov = alloc2D_d(K, K);

  // scale array oldProbs so that elements sum to 1, and accumulate row and col sums
  for( int j = 0; j <  States; ++j ) {
    Sum += oldProbs[j];
  }

  scaleFactor = 1.0 / Sum;
  //accumulate sum of logs of scale factor to correct log likelihood
  //sumfactor += log(scaleFactor);

  for( int j0 = 0; j0 <  K; ++j0 ) {
    rowProb[j0] = 0.0;
    colProb[j0] = 0.0;
    for( int j1 =0; j1 < K; ++j1 ) {
      rowProb[j0] += oldProbs[j0*K + j1] * scaleFactor;
      colProb[j0] += oldProbs[j1*K + j0] * scaleFactor;
      //oldProbs[j0][j1] *= scaleFactor; 
    }
  }
  // calculate expectations of indicator variables for each ancestry state on each gamete
  for( int j = 0; j <  K; ++j ) {
    Expectation0[j] = f[0]*rowProb[j] + stateArrivalProbs[j*2];
    Expectation1[j] = f[1]*colProb[j] + stateArrivalProbs[j*2 + 1];
    
  }
  // calculate covariance of ancestry states as ff * deviation from product of row and col probs
  for(int j0 = 0; j0 <  K-1; ++j0) { // leave out last row
    for(int j1 =0; j1 < K-1; ++j1) { // leave out last col
      cov[j0][j1] = ff * ( oldProbs[j0*K + j1]*scaleFactor - rowProb[j0] * colProb[j1] );
    }
  }

  // accumulate sums of covariances over first K-1 rows and K-1 cols
  for(int j0 = 0; j0 <  K-1; ++j0) { // leave out last row
    rowSum[j0] = 0.0;
    colSum[j0] = 0.0;
    for(int j1 =0; j1 < K-1; ++j1) { // leave out last col
      rowSum[j0] += cov[j0][j1];
      colSum[j0] += cov[j1][j0];
    }
  }
  // calculate last row except for last col, by subtracting colSum from 0
  // also accumulate sum of covariances for K th row over first K-1 cols
  rowSum[K-1] = 0.0;
  for( int j = 0; j < K-1; ++j ) {
    cov[K-1][j] =  -colSum[j];
    rowSum[K-1] += cov[K-1][j];
  }
  // calculate last col by subtracting rowSum from 0
  for( int j = 0; j < K; ++j ) {
    cov[j][K-1] =  -rowSum[j];
  }

  // calculate expectation of product as covariance plus product of expectations
  // can speed up with Fourier transform 
  for(int j0 = 0; j0 < K; ++j0) {
    for(int j1 =0; j1 < K; ++j1) {
      
      newProbs[j0*K + j1] = cov[j0][j1] + 
	// newProbs[j0][j1] = ff * (oldProbs[j0*K + j1]*scaleFactor - rowProb[j0] * colProb[j1]) + 
	Expectation0[j0] * Expectation1[j1];
      //	 ( f[0]*rowProb[j0] + stateArrivalProbs[j0*2] ) * ( f[1]*colProb[j1] + stateArrivalProbs[j1*2 + 1] );
      //undo scaling 
      newProbs[j0*K + j1] *= Sum;
    }
  }
  delete[] rowProb;
  delete[] colProb;
  delete[] Expectation0;
  delete[] Expectation1;

  delete[] rowSum;
  delete[] colSum;
  free_matrix(cov, K);
  
}

void HMM::SampleJumpIndicators(int *LocusAncestry, double *f[], const unsigned int gametes, 
			       const Vector &Distances, const int startLocus, 
			       int *sumxi, double *Sumrho0, Matrix_i *SumLocusAncestry, Matrix_i *SumLocusAncestry_X, bool isX, 
			       unsigned int SumN[], unsigned int SumN_X[], bool RhoIndicator){

  int locus;
  double Prob;
  bool xi[2][Transitions];//jump indicators
  xi[0][0] = xi[1][0] = true;

  for( int jj = 1; jj < Transitions; jj++ ){
    locus = startLocus + jj;
    xi[0][jj] = xi[1][jj] = true;    
    for( unsigned int g = 0; g < gametes; g++ ){
      if( LocusAncestry[g*Transitions + jj-1] == LocusAncestry[jj + g*Transitions] ){

	Prob = StateArrivalProbs[jj*K*2 +LocusAncestry[jj + g*Transitions]*2 + g] / 
	  (StateArrivalProbs[jj*K*2 + LocusAncestry[jj + g*Transitions]*2 +g] + f[g][jj] );
	if( Prob > myrand() ){
	  xi[g][jj] = true;
	  sumxi[locus]++;
	} else {
	  xi[g][jj] = false;
	  *Sumrho0 += Distances( jj );
	}
      } else {
	xi[g][jj] = true;
	sumxi[locus]++;
      }
 
      if( xi[g][jj] ){
	// sum ancestry states over loci where jump indicator is 1
	if( !isX )
	  (*SumLocusAncestry)( LocusAncestry[jj + g*Transitions], g )++;
	else
	  (*SumLocusAncestry_X)( LocusAncestry[jj + g*Transitions], g )++;
	//sample number of arrivals where jump indicator is 1
	if(RhoIndicator){
	  double u = myrand();
	  // sample distance dlast back to last arrival, as dlast = -log[1 - u(1-f)] / rho
	  // then sample number of arrivals before last as Poisson( rho*(d - dlast) )
	  // algorithm does not require rho or d, only u and f
	  unsigned int sample = genpoi( log( (1 - u*( 1 - f[g][jj])) / f[g][jj] ) );
	  if( !isX )
	    SumN[g] += sample + 1;
	  else
	    SumN_X[g] += sample + 1;
	}
      }
    }
  }
  //finally for first locus, not include in above loop
    for( unsigned int g = 0; g < gametes; g++ ){
      if( xi[g][0] ){
	if( !isX )
	  (*SumLocusAncestry)( LocusAncestry[g*Transitions], g )++;
	else
	  (*SumLocusAncestry_X)( LocusAncestry[g*Transitions], g )++;
      }
    }
}
