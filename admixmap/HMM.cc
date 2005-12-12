#include "rand.h"
#include "HMM.h"
#include "functions.h"

using namespace std;

HMM::HMM()
{
  alpha = 0;
  beta = 0;
  LambdaBeta = 0;
  p = 0;
  StateArrivalProbs = 0;
  ThetaThetaPrime = 0;
  colProb = 0;
  Expectation0 = 0;
  Expectation1 = 0;
  rowSum = 0;
  colSum = 0;
  rowProb = 0;
  cov = 0;
}

//not currently used
//HMM objects are instantiated in Chromosome using default constructor above
//and dimensions set by SetDimensions below
HMM::HMM( int inTransitions, int pops)
{
  SetDimensions(inTransitions, pops);
}

HMM::~HMM()
{
  delete[] p;
  delete[] LambdaBeta;
  delete[] alpha;
  delete[] beta;
  delete[] StateArrivalProbs;
  delete[] ThetaThetaPrime;
  delete[] rowProb;
  delete[] colProb;
  delete[] Expectation0;
  delete[] Expectation1;
  delete[] rowSum;
  delete[] colSum;
  free_matrix(cov, K);
}

void HMM::SetDimensions( int inTransitions, int pops)
{
  //inTransitions = #transitions +1 = #Loci 
  //pops = #populations
  K = pops;

  Transitions = inTransitions;

  alpha = new double[Transitions*K*K];
  beta =  new double[Transitions*K*K];
  
  sumfactor=0.0;
  p = new double[Transitions];
  LambdaBeta = new double[K*K];

  StateArrivalProbs = new double[Transitions * K * 2];
  ThetaThetaPrime = new double[K*K];
  if(K>2){
    rowProb = new double[K];
    colProb = new double[K];
    Expectation0 = new double[K];
    Expectation1 = new double[K];
    rowSum = new double[K];
    colSum = new double[K];
    cov = alloc2D_d(K,K);
  }
}

void HMM::SetStateArrivalProbs(const double* const f[], const double* const Theta, int Mcol){
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
  Updates Forward probabilities alpha and array p (=f0*f1)
  diploid case only
  -----------------------------------------------------------
  relates to notation in docs as follows:
  j = j1, j' = j2, i = i1, i' = i2, 
  theta_j = Admixture(j,0),  theta'_j = Admixture(j,Mcol).
*/
void HMM::UpdateForwardProbsDiploid(const double* const f[], double* lambda, double coolness)
{
  sumfactor = 0.0;
  double scaleFactor, Sum;
  States = K*K;
  double *alphat;
  alphat = alpha;
  
  for(int j = 0; j < States; ++j) {
    //set alpha[0] = StationaryDist * lambda(0)
    if(coolness > 1.0) 
      alpha[j] =  ThetaThetaPrime[j] * pow(*lambda++, coolness);
    else  {
      *(alphat + j) =  ThetaThetaPrime[j] * *lambda++;
    }
  }
  
  for( int t = 1; t < Transitions; t++ ){        
    Sum = 0.0;
    //scale previous alpha to sum to 1
    for( int j = 0; j <  States; ++j ) {
      Sum  += *(alphat + j);
    }
    scaleFactor = 1.0 / Sum;
    for( int j = 0; j <  States; ++j ) {
      *(alphat + j) *= scaleFactor;
    }
    sumfactor += log(Sum);
    
    double f2[2] = {f[0][t], f[1][t]};
    p[t] = f2[0] * f2[1];
    RecursionProbs(p[t], f2, StateArrivalProbs + t*K*2, alphat, alphat + States);
    alphat += States;
    
    for(int j = 0; j < States; ++j){
      if(coolness > 1.0)
	alpha[t*States +j] *= pow(*lambda++, coolness);
      else {
	*(alphat +j) *= *lambda++;
      }
    }
  }
}
  
void HMM::UpdateBackwardProbsDiploid(const double* const f[], const double* const lambda)
{
  vector<double> rec(States);
  double scaleFactor, Sum;
  
  for(int j = 0; j < States; ++j){
    //set beta(T) = 1
    beta[(Transitions - 1)*States + j] = 1.0;
    rec[j] = 1.0 / ThetaThetaPrime[j];
  }
  
  for( int t = Transitions-2; t >=0; t-- ){
    
    double f2[2] = {f[0][t+1], f[1][t+1]};
    
    Sum = 0.0;    
    for(int j = 0; j < States; ++j){
      LambdaBeta[j] = lambda[(t+1)*K*K + j] * beta[(t+1)*States + j] * ThetaThetaPrime[j];
      Sum += LambdaBeta[j];
    }
    //scale LambdaBeta to sum to 1
    scaleFactor = 1.0 / Sum;
    for( int j = 0; j <  States; ++j ) {
      LambdaBeta[j] *= scaleFactor;
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
void HMM::UpdateForwardProbsHaploid(const double* const f[], const double* const Admixture, const double* const lambda){

  sumfactor = 0.0;
  //double factor = 0.0;
  double Sum;
  States = K;

  for(int j = 0; j < States; ++j){
    alpha[j] = Admixture[j] * lambda[j];
    beta[(Transitions-1)*States + j] = 1.0;
  }

  for( int t = 1; t < Transitions; t++ ){
    Sum = 0.0;
    for(int j = 0; j < States; ++j){
      Sum += alpha[(t-1)*States + j];
    }
    //factor = 0.0;
    for(int j = 0; j < States; ++j){
      alpha[t*States + j] = f[0][t] + (1.0 - f[0][t]) * Admixture[j] * Sum;
      alpha[t*States + j] *= lambda[(t+1)*States + j];
      //factor += alpha[t*States + j];
    }
  }
  //TODO: rescale to avoid underflow
}

void HMM::UpdateBackwardProbsHaploid(const double* const f[], const double* const Admixture, const double* const lambda){
  double Sum;
  for(int j = 0; j < States; ++j){
    beta[(Transitions-1)*States + j] = 1.0;
  }
  
  for( int t = Transitions-2; t >=0; t-- ){
    Sum = 0.0;
    for(int j = 0; j < States; ++j){
      Sum += Admixture[j]*lambda[(t+1)*States + j]*beta[(t+1)*States + j];
    }
    for(int j=0;j<States;++j){
      beta[t*States + j] = f[t+1][0]*lambda[(t+1)*States + j]*beta[(t+1)*States + j] + (1.0 - f[0][t+1])*Sum;
    }
  }
}

/*
  computes conditional state probabilities at "time" t
  probs - double array to store state probs
  K = #populations
*/
// void HMM::GetStateProbs( double * probs, int t)const
// {
//   double sum = 0.0;
//   int State = 0;

//   for(int i = 0; i < K; ++i)
//    for( int j = 0; j < K; j++ ){
//      probs[State++] = alpha[t*States + i*K +j] * beta[t*States + i*K +j];
//      sum += probs[State-1];
//    }

//    for( int j = 0; j < States; j++ ){
//      probs[j] /= sum;
//    }
// }

std::vector<std::vector<double> > HMM::Get3WayStateProbs( int t)const{
  double sum = 0.0;
  int State = 0;
  std::vector<double> probs(States);
  std::vector<std::vector<double> >AncestryProbs(3);

  for( int j = 0; j < States; j++ ){
    probs[State++] = alpha[t*States + j] * beta[t*States + j];
    sum += probs[State-1];
  }

   for( int j = 0; j < States; j++ ){
     probs[j] /= sum;
   }
   for( int k1 = 0; k1 < K; k1++ ){
     AncestryProbs[2].push_back(probs[ ( K + 1 ) * k1 ]);
     AncestryProbs[1].push_back( 0.0 );
     for( int k2 = 0 ; k2 < K; k2++ )
       AncestryProbs[1][k1] += probs[k1*K +k2] + probs[k2*K +k1];
     AncestryProbs[1][k1] -= 2.0*AncestryProbs[2][k1];
     AncestryProbs[0].push_back( 1.0 - AncestryProbs[1][k1] - AncestryProbs[2][k1] );
   }
   return AncestryProbs;
}

/*
  returns log-likelihood
  This is the sum over states of products of alpha and beta
  and is the same for all t so it is convenient to compute for
  t = T-1 since beta here is 1. 
*/
double HMM::getLogLikelihood()const
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
void HMM::Sample(int *SStates, const double* const Admixture, const double* const f[], bool isdiploid)const
{
  int j1,j2;
  double* V = new double[States];
  int* C = new int[Transitions];

  if(isdiploid){
    int State = 0;
    for( int j = 0; j < States; j++ )V[State++] = alpha[(Transitions - 1)*States + j];
    C[ Transitions - 1 ] = SampleFromDiscrete( V, States );
    SStates[Transitions-1] = (int)(C[Transitions-1]/K);
    SStates[Transitions - 1 + Transitions] = (C[Transitions-1] % K);
    
    for( int t =  Transitions - 2; t >= 0; t-- ){
      j1 = (int) (C[t+1]/K);//j
      j2 = C[t+1]-K*j1;     //j'

      State = 0;
      for(int i1 = 0; i1 < K; i1++)for(int i2 = 0; i2 < K; ++i2){
	V[State] = 
	  ( (i1==j1)*f[0][t+1] + StateArrivalProbs[(t+1)*K*2 + j1*2] ) * ( (i2==j2)*f[1][t+1] + StateArrivalProbs[(t+1)*K*2 + j2*2 +1] );
	V[State] *= alpha[t*States + i1*K + i2];
	State++;
      }
      C[ t ] = SampleFromDiscrete( V, States );
      SStates[t] = (int)(C[t]/K);//paternal
      SStates[t + Transitions] = (C[t] % K);//maternal
    }
   }
  else{//haploid
    for( int j = 0; j < States; j++ )V[j] = alpha[(Transitions - 1)*States + j ];
    C[ Transitions - 1 ] = SampleFromDiscrete( V, States );
    SStates[Transitions-1] = C[Transitions-1];
    for( int t =  Transitions - 2; t >= 0; t-- ){
      for(int j = 0; j < States; j++)V[j] = (j == C[t+1])*f[0][t+1]+Admixture[C[t+1]]*(1.0 - f[0][t]);
      for( int j = 0; j < States; j++ )	V[j] *= alpha[t*States + j];
      C[ t ] = SampleFromDiscrete( V, States );
      SStates[t] = C[t];
    }
  }
  delete[] V;
  delete[] C;
}

// argument oldProbs is square array of size K, K
// for forward recursions, pass alpha_t and multiply newProbs by observation probs 
// for backward recursions, pass array of products lambda_t+1[jj] * beta_t+1[jj] 
// updates oldProbs (scaled to sum to 1), newProbs and sumfactor if forward = true (for alphas)
void HMM::RecursionProbs(const double ff, const double f[2], 
			 const double* const stateArrivalProbs, double* oldProbs, double *newProbs) {
  if(K==2) RecursionProbs2(ff, f, stateArrivalProbs, oldProbs, newProbs);
  else
    {
      
      for( int j0 = 0; j0 <  K; ++j0 ) {
	rowProb[j0] = 0.0;
	colProb[j0] = 0.0;
	for( int j1 =0; j1 < K; ++j1 ) {
	  rowProb[j0] += oldProbs[j0*K + j1];
	  colProb[j0] += oldProbs[j1*K + j0];
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
	  cov[j0][j1] = ff * ( oldProbs[j0*K + j1] - rowProb[j0] * colProb[j1] );
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
      for(int j0 = 0; j0 < K; ++j0) {
	for(int j1 =0; j1 < K; ++j1) {
	  newProbs[j0*K + j1] = cov[j0][j1] + Expectation0[j0] * Expectation1[j1];
	  // newProbs[1] is prob(paternal=1, maternal=0)
	}
      }
    }//end else
}

void HMM::RecursionProbs2(const double ff, const double f[2], const double* const stateArrivalProbs, 
			  const double* const oldProbs, double *newProbs) {
  // version for 2 subpopulations
  double row0Prob;
  double col0Prob;
  double Exp0;
  double Exp1;
  double cov0;
  double Product;
  // sum row 0 and col 0  
  row0Prob = ( oldProbs[0] + oldProbs[2] );
  col0Prob = ( oldProbs[0] + oldProbs[1] );
  // calculate expectations of indicator variables for ancestry=0 on each gamete
  Exp0 = f[0]*row0Prob + stateArrivalProbs[0]; // paternal gamete
  Exp1 = f[1]*col0Prob + stateArrivalProbs[1]; // maternal gamete
  Product = Exp0 * Exp1;
  // calculate covariance of indicator variables as ff * deviation from product of row and col probs
  cov0 = ff * ( oldProbs[0] - row0Prob * col0Prob );
  newProbs[0] = Product + cov0; 
  newProbs[1] = Exp0 - newProbs[0]; //prob paternal ancestry=1, maternal=0 
  newProbs[2] = Exp1 - newProbs[0]; //prob paternal ancestry=0, maternal=1 
  newProbs[3] = 1 - Exp0 - Exp1 + newProbs[0];
}

void HMM::SampleJumpIndicators(const int* const LocusAncestry, const double* const f[], const unsigned int gametes, 
			       int *SumLocusAncestry, int *SumLocusAncestry_X, bool isX, 
			       unsigned SumN[], unsigned SumN_X[], bool isGlobalRho)const{

  double Prob;
  vector<bool> xi[2] = {vector<bool>(Transitions), vector<bool>(Transitions)};//jump indicators
  xi[0][0] = xi[1][0] = true;
  for( int jj = 1; jj < Transitions; jj++ ){
    xi[0][jj] = xi[1][jj] = true;    
    for( unsigned int g = 0; g < gametes; g++ ){
      if( LocusAncestry[g*Transitions + jj-1] == LocusAncestry[jj + g*Transitions] ){
	Prob = StateArrivalProbs[jj*K*2 +LocusAncestry[jj + g*Transitions]*2 + g];  
	xi[g][jj] = Prob / (Prob + f[g][jj]) > myrand();
      } else {
	xi[g][jj] = true;
      }
 
      if( xi[g][jj] ){
	// sum ancestry states over loci where jump indicator is 1
	if( !isX )
	  SumLocusAncestry[ LocusAncestry[jj + g*Transitions] +  g*K ]++;
	else
	  SumLocusAncestry_X[ LocusAncestry[jj + g*Transitions] + g*K ]++;
	//sample number of arrivals where jump indicator is 1
	if(!isGlobalRho){
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
	  SumLocusAncestry[ LocusAncestry[g*Transitions] + g*K ]++;
	else
	  SumLocusAncestry_X[ LocusAncestry[g*Transitions] + g*K] ++;
      }
    }
}

