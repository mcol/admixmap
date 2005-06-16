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

  alpha = new double**[Transitions];
  beta =  new double**[Transitions];
  int d = isdiploid? K:1;  
  
  for(int i=0;i < Transitions;++i) {
    alpha[i] = new double*[K];
    beta[i]  = new double*[K];
    for(int j=0;j<K;++j){
      alpha[i][j]= new double[d];
      beta[i][j] = new double[d];
    }
  }
  sumfactor=0.0;
  p = new double[Transitions];
  LambdaBeta = alloc2D_d(K,d);
}

HMM::~HMM()
{
  //TODO:destroy these properly
  delete[] p;
  for(int i = 0; i < K; ++i){
    delete[] LambdaBeta[i];
  }
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
  //TransitionProbs = new double**[ Transitions - 1];
  //for(int i=0;i< Transitions - 1; i++){
  // TransitionProbs[i] = alloc2D_d(States, States);
  //}

  alpha = new double**[Transitions];
  beta =  new double**[Transitions];
  int d = isdiploid? K:1;  
  
  for(int i=0;i < Transitions;++i) {
    alpha[i] = alloc2D_d(K,d);
    beta[i] = alloc2D_d(K,d);
  }
  sumfactor=0.0;
  p = new double[Transitions];
  LambdaBeta = alloc2D_d(K,d);
}



/*
  Updates Forward and (if required) Backward probabilities
  diploid case only
  -----------------------------------------------------------
  Mcol = options->getRandomMatingModel() = 1 for random mating model, 0 otherwise
  K = # populations
  
  relates to notation in docs as follows:
  j = j1, j' = j2, i = i1, i' = i2, 
  theta_j = Admixture(j,0),  theta'_j = Admixture(j,Mcol).
*/
void HMM::UpdateForwardProbsDiploid(double ***StateArrivalProbs, double *f[], Matrix_d& Admixture, double ***lambda, int Mcol)
{
  sumfactor = 0.0;

   for(int j1 = 0; j1 < K; ++j1)
    for(int j2 =0; j2 <K; ++j2){
      //set alpha(0) = StationaryDist * lambda(0) 
      alpha[0][j1][j2] =  Admixture(j1,0)*Admixture(j2,Mcol) *lambda[0][j1][j2];
    }

  for( int t = 1; t < Transitions; t++ ){        
    p[t] = f[0][t] * f[1][t];
    double f2[2] = {f[0][t], f[1][t]};

    RecursionProbs(p[t], f2,StateArrivalProbs[t], alpha[t-1], alpha[t], true);
    for(int j1 = 0; j1 < K; ++j1){
      for(int j2 = 0; j2 < K; ++j2){
	alpha[t][j1][j2] *= lambda[t][j1][j2];
      }
    }
  }
}

void HMM::UpdateBackwardProbsDiploid(double ***StateArrivalProbs,double *f[], double **ThetaThetaPrime, double ***lambda)
{
  double rec[K][K];

  for(int j1 = 0; j1 < K; ++j1)
    for(int j2 =0; j2 <K; ++j2){
      //set beta(T) = 1
      beta[Transitions - 1][j1][j2] = 1.0;
      rec[j1][j2] = 1.0 / ThetaThetaPrime[j1][j2];
    }

  for( int t = Transitions-2; t >=0; t-- ){
    
    double f2[2] = {f[0][t+1], f[1][t+1]};
    
    for(int j1 = 0; j1 < K; ++j1){
	for(int j2 =0; j2 <K; ++j2){
	  LambdaBeta[j1][j2] = lambda[t+1][j1][j2] * beta[t+1][j1][j2] * ThetaThetaPrime[j1][j2];
	}
    }
    
    RecursionProbs(p[t+1], f2,StateArrivalProbs[t+1], LambdaBeta, beta[t], false);
    for(int j1 = 0; j1 < K; ++j1){
      for(int j2 =0; j2 <K; ++j2){
	beta[t][j1][j2] *= rec[j1][j2];
      }
    }
  }

}


/*
  Updates Forward and (if required) Backward probabilities
  haploid case only
  Here Admixture is a column matrix and the last dimensions of f and lambda are 1.
*/
void HMM::UpdateProbsHaploid(double *f[], Matrix_d& Admixture, double ***lambda, bool CalculateBeta){

  double Sum;

  for(int j=0; j<States;++j){
    alpha[0][j][0] = Admixture(j,0) * lambda[0][j][0];
    beta[Transitions-1][j][0] = 1.0;
  }

  for( int t = 1; t < Transitions; t++ ){
    Sum = 0.0;
    for(int j=0;j<States;++j){
      Sum += alpha[t-1][j][0];
    }
    for(int j=0;j<0;++j){
      alpha[t][j][0] = f[0][t] + (1.0 - f[0][t]) * Admixture(0,j) * Sum;
      alpha[t][j][0] *= lambda[t+1][j][0];
    }

    if(CalculateBeta){
      for( int t = Transitions-2; t >=0; t-- ){
	Sum = 0.0;
	for(int j=0;j<States;++j){
	  Sum += Admixture(0,j)*lambda[t+1][j][0]*beta[t+1][j][0];
	}
	for(int j=0;j<0;++j){
	  beta[t][j][0] = f[t+1][0]*lambda[t+1][j][0]*beta[t+1][j][0] + (1.0 - f[0][t+1])*Sum;
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
     probs[State++] = alpha[t][i][j] * beta[t][i][j];
     sum += probs[State-1];
   }

   for( int j = 0; j < States; j++ ){
     probs[j] /= sum;
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
  for(int i = 0; i < K; ++i)
   for( int j = 0; j < K; j++ ){
     sum += alpha[Transitions - 1][i][j];
   }
   return( sumfactor+log( sum ) );
}

// void HMM::PrintLikelihoods(){
//   double oldsum=0.0, newsum=0.0;
//   bool flag = false;
//   cout<<"old lik"<<" "<<"new lik"<<endl;
//   for(int t = 0; t< Transitions;++t){
//     oldsum = newsum = 0.0;
//     for(int i=0;i<K;++i)
//       for(int j=0;j<K;++j){
// 	oldsum += beta[t][i][j]*alpha[t][i][j];
// 	newsum += beta[t][i][j]*alpha[t][i][j];
//       }
//     cout<<"t = "<<t<<": "<<oldsum<<" "<<newsum<<endl;
//     if(oldsum - newsum < -0.000001 || oldsum - newsum > 0.000001)flag=true;
//   }
//   if(flag)system("pause");
// }

/*
  Samples Hidden States
  ---------------------
  C          - an int array to store the sampled states
  isdiploid  - indicator for diploidy
*/
void HMM::Sample(int *C, Matrix_d &Admixture, double *f[], double ***StateArrivalProbs, bool isdiploid)
{
  int j1,j2;
  double V[States];

  if(isdiploid){
    int State = 0;
    for( int i1 = 0; i1 < K; i1++ )for(int i2 = 0; i2<K; ++i2)V[State++] = alpha[Transitions - 1][i1][i2];
    C[ Transitions - 1 ] = SampleFromDiscrete3( V, States );
    
    for( int t =  Transitions - 2; t >= 0; t-- ){
      j1 = (int) (C[t+1]/K);//j
      j2 = C[t+1]-K*j1;     //j'

      State = 0;
      for(int i1 = 0; i1 < K; i1++)for(int i2=0;i2<K;++i2){
	V[State] = 
	  ( (i1==j1)*f[0][t+1] + StateArrivalProbs[t+1][j1][0] ) * ( (i2==j2)*f[1][t+1] + StateArrivalProbs[t+1][j2][1] );
	V[State] *= alpha[t][i1][i2];
	State++;
      }
      C[ t ] = SampleFromDiscrete3( V, States );
    }
  }
  else{//haploid
    for( int j = 0; j < States; j++ )V[j] = alpha[Transitions - 1][j][0];
    C[ Transitions - 1 ] = SampleFromDiscrete3( V, States );
    for( int t =  Transitions - 2; t >= 0; t-- ){
      for(int j = 0; j < States; j++)V[j] = (j == C[t+1])*f[0][t+1]+Admixture(C[t+1],0)*(1.0 - f[0][t]);
      for( int j = 0; j < States; j++ )	V[j] *= alpha[t][j][0];
      C[ t ] = SampleFromDiscrete3( V, States );
   }
  }

}

// argument oldProbs is square array of size K, K
// for forward recursions, pass alpha_t and multiply newProbs by observation probs 
// for backward recursions, pass array of products lambda_t+1[jj] * beta_t+1[jj] 
// updates oldProbs (scaled to sum to 1), newProbs and sumfactor if forward = true (for alphas)
void HMM::RecursionProbs(const double ff, const double f[2], 
			 //const double stateArrivalProbs[][2], double **oldProbs,  
			 double **stateArrivalProbs, double **oldProbs,
			 double **newProbs, bool forward) {
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
  for( int j0 = 0; j0 <  K; ++j0 ) {
    for( int j1 =0; j1 < K; ++j1 ) { 
      Sum += oldProbs[j0][j1];
    }
  }

  scaleFactor = 1.0 / Sum;
  //accumulate sum of logs of scale factor to correct log likelihood
  //sumfactor += log(scaleFactor);

  for( int j0 = 0; j0 <  K; ++j0 ) {
    rowProb[j0] = 0.0;
    colProb[j0] = 0.0;
    for( int j1 =0; j1 < K; ++j1 ) {
      rowProb[j0] += oldProbs[j0][j1] * scaleFactor;
      colProb[j0] += oldProbs[j1][j0] * scaleFactor;
      //oldProbs[j0][j1] *= scaleFactor; 
    }
  }
  // calculate expectations of indicator variables for each ancestry state on each gamete
  for( int j = 0; j <  K; ++j ) {
    Expectation0[j] = f[0]*rowProb[j] + stateArrivalProbs[j][0];
    Expectation1[j] = f[1]*colProb[j] + stateArrivalProbs[j][1];
    
  }
  // calculate covariance of ancestry states as ff * deviation from product of row and col probs
  for(int j0 = 0; j0 <  K-1; ++j0) { // leave out last row
    for(int j1 =0; j1 < K-1; ++j1) { // leave out last col
      cov[j0][j1] = ff * ( oldProbs[j0][j1]*scaleFactor - rowProb[j0] * colProb[j1] );
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
      
      newProbs[j0][j1] = cov[j0][j1] + 
	// newProbs[j0][j1] = ff * (oldProbs[j0][j1]*scaleFactor - rowProb[j0] * colProb[j1]) + 
	Expectation0[j0] * Expectation1[j1];
      //	 ( f[0]*rowProb[j0] + stateArrivalProbs[j0][0] ) * ( f[1]*colProb[j1] + stateArrivalProbs[j1][1] );
      //if(forward)
      //oldProbs[j0][j1] *= scaleFactor;
      //undo scaling 
      newProbs[j0][j1] *= Sum;
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


