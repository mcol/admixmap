#include "rand.h"
#include "vector_i.h"
#include "vector_d.h"
#include "matrix_d.h"
#include "MatrixArray_d.h"
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
  //TransitionProbs = new double**[ Transitions - 1];
  //for(int i=0;i< Transitions - 1; i++){
  // TransitionProbs[i] = new double *[States];
  //for(int j=0;j < States; j++) TransitionProbs[i][j] = new double[States];
  //}
  // forward probabilities
  //alpha = new double*[Transitions];
  // backward probabilities
  //beta = new double*[Transitions];

  newalpha = new double**[Transitions];
  newbeta =  new double**[Transitions];
  int d = isdiploid? K:1;  
  
  for(int i=0;i < Transitions;++i) {
    //alpha[i] = new double[States];
    //beta[i] =  new double[States];
    newalpha[i] = new double*[K];
    newbeta[i]  = new double*[K];
    
    for(int j=0;j<K;++j){
      newalpha[i][j]= new double[d];
      newbeta[i][j] = new double[d];
    }
  }
  sumfactor=0.0;
  p = new double[Transitions];
  q = new double[Transitions];  
  r = new double[Transitions];
  s = new double[Transitions];
  Q = alloc2D_d(Transitions,K);
  R = alloc2D_d(Transitions,K);
  S = new double **[Transitions];
  for(int i=0;i<Transitions;++i)S[i] = alloc2D_d(K,K);
  U = alloc2D_d(K,K);
  V = alloc2D_d(K,K);
  W = alloc2D_d(K,K);
}

HMM::~HMM()
{
  //TODO:destroy these properly
  //delete[] alpha;
  //delete[] beta;
  //delete[] TransitionProbs;
  delete[] p;
  delete[] q;
  delete[] r;
  delete[] s;
  delete[] Q;
  delete[] R;
  for(int i=0;i<K;++i){
    delete[] S[i];
    delete[] U[i];
    delete[] V[i];
    delete[] W[i];
  }
  delete[] S;
  delete[] U;
  delete[] V;
  delete[] W;
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

  //alpha = new double*[Transitions];
  //beta = new double*[Transitions];

  newalpha = new double**[Transitions];
  newbeta =  new double**[Transitions];
  int d = isdiploid? K:1;  
  
  for(int i=0;i < Transitions;++i) {
    //alpha[i] = new double[States];
    //beta[i] =  new double[States];
    newalpha[i] = alloc2D_d(K,d);
    newbeta[i] = alloc2D_d(K,d);
  }
  sumfactor=0.0;
  p = new double[Transitions];
  q = new double[Transitions];  
  r = new double[Transitions];
  s = new double[Transitions];
  Q = alloc2D_d(Transitions,K);
  R = alloc2D_d(Transitions,K);
  S = new double **[Transitions];
  for(int i=0;i<Transitions;++i)S[i] = alloc2D_d(K,K);
  U = alloc2D_d(K,K);
  V = alloc2D_d(K,K);
  W = alloc2D_d(K,K);
}

//Update Forward and (optionally) Backward probs
// void HMM::UpdateFwrdBckwdProbabilities( double *StationaryDist, double **StateDepProbs, bool CalculateBeta)
// {
//   //StationaryDist = delta, double array of length States 
//   //TransitionProbs = Gamma, (States x States) Transition Matrix
//   //StateDepProbs = lambda = (States x States) matrix of State-Dependent Probabilities, p(Hidden State | Observed State)
//   //CalculateBeta indicates whether to calculate backward probs beta. They are only needed for state probs.
//   sumfactor = 0.0;
//    //set alpha(0) = StationaryDist * lambda(0) 
//    for( int j = 0; j < States; j++ )
//      alpha[0][j] =  StationaryDist[j] *StateDepProbs[0][j];

//    for( int t = 1; t < Transitions; t++ ){
//      //set alpha(t) = alpha(t-1) * Gamma * lambda(t)
//      for( int i = 0; i < States; i++ ){
//        alpha[t][i] = 0.0;
//        //explicit element-wise vector-matrix multiplication
//        for( int j = 0; j < States; j++ )alpha[t][i] += alpha[t-1][j] * TransitionProbs[t-1][j][i];
//        alpha[t][i] *= StateDepProbs[t][i];
//      }

//      //This is to avoid underflow
//      //probably needs adjusting
//      factor = 0.0; 
//      //       if( !(t%5) ){
//  //     for(int i=0;i<States;++i)     
// //           factor += alpha[t][i];
// //      sumfactor += log(factor);
// //      for(int i=0;i<States;++i)     
// //                alpha[t][i] /= factor;
//      //      }
//     }

//    if( CalculateBeta ){
//      //set beta(T) = 1
//      for(int i = 0; i < States; ++i)beta[Transitions - 1][i] = 1.0;
     
//      for( int t = Transitions - 2; t >= 0; t-- ){
//        for( int i = 0; i < States; i++ ){
// 	 beta[t][i] = 0.0;
// 	 //set beta(t) to Gamma * lambda(t+1) * beta(t+1) 
// 	 for( int j = 0; j < States; j++ )beta[t][i] += StateDepProbs[t+1][j] * beta[t+1][j] * TransitionProbs[t][i][j];
//        }	
//      }
//    }

// }


/*
  Updates Forward and (if required) Backward probabilities
  diploid case only
  -----------------------------------------------------------
  Mcol = options->getRandomMatingModel() = 1 for random mating model, 0 otherwise
  K = # populations
  StartLocus = number of first locus on chromosome
  
  relates to notation in docs as follows:
  j = j1, j' = j2, i = i1, i' = i2, 
  theta_j = Admixture(j,0),  theta'_j = Admixture(j,Mcol).
*/
void HMM::UpdateProbsDiploid(double *f[],int StartLocus, Matrix_d& Admixture, double ***lambda, int Mcol, bool CalculateBeta)
{
  sumfactor=0.0;
  double Sum1, Sum2, Sum3;

  //First t = 0
  for(int j1 = 0; j1 < K; ++j1)
    for(int j2 =0; j2 <K; ++j2){
      //set alpha(0) = StationaryDist * lambda(0) 
      newalpha[0][j1][j2] =  Admixture(j1,0)*Admixture(j2,Mcol) *lambda[0][j1][j2];
    //set beta(T) = 1
      newbeta[Transitions - 1][j1][j2] = 1.0;
    }
  p[0] = f[0][StartLocus] * f[1][StartLocus]; 
  q[0] = f[1][StartLocus] - p[0];
  r[0] = f[0][StartLocus] - p[0];
  s[0] = 1.0 - p[0] - q[0] - r[0];
  for(int j=0;j<K;++j){
    Q[0][j] = q[0] * Admixture(j,0);
    R[0][j] = r[0] * Admixture(j,Mcol);
    for(int jj=0;jj<K;++jj)S[0][j][jj] = s[0] * Admixture(j,0) * Admixture(jj,Mcol);
  }

  for( int t = 1; t < Transitions; t++ ){           //probs of
    p[t] = f[0][StartLocus+t] * f[1][StartLocus+t]; //no arrivals
    q[t] = f[1][StartLocus+t] - p[t]; //arrival on one gamete
    r[t] = f[0][StartLocus+t] - p[t];//arrival on other gamete
    s[t] = 1.0 - p[t] - q[t] - r[t];//arrival on both gametes

    for(int j=0;j<K;++j){
      Q[t][j] = q[t] * Admixture(j,0);
      R[t][j] = r[t] * Admixture(j,Mcol);
      for(int jj=0;jj<K;++jj)S[t][j][jj] = s[t] * Admixture(j,0) * Admixture(jj,Mcol);
    }

    Sum3 = 0.0;
    for(int i1=0;i1<K;++i1)
      for(int i2=0;i2<K;++i2) Sum3 += newalpha[t-1][i1][i2];

    for(int j1 = 0; j1 < K; ++j1){
      for(int j2 =0; j2 <K; ++j2){
	Sum1 = Sum2 = 0.0;
	for(int i1=0;i1<K;++i1){
	  Sum1 += newalpha[t-1][i1][j2];
	  Sum2 += newalpha[t-1][j1][i1];//using i1 for i' here to save another loop
	}
	U[j1][j2] = Q[t][j1]*Sum1;
	V[j1][j2] = R[t][j2]*Sum2;
	W[j1][j2] = S[t][j1][j2]*Sum3;
      }
    }

    //Forward Probs    
    for(int j1 = 0; j1 < K; ++j1){
      for(int j2 =0; j2 <K; ++j2){
	newalpha[t][j1][j2] = p[t] * newalpha[t-1][j1][j2]  +U[j1][j2] +V[j1][j2] + W[j1][j2];
	newalpha[t][j1][j2] *= lambda[t][j1][j2];
      }
    }//end alpha calculation

     factor = 0.0; 
     //       if( !(t%5) ){
//      for(int i=0;i<K;++i)     for(int j=0;j<K;++j)
//           factor += newalpha[t][i][j];
//      sumfactor += log(factor);
//      for(int i=0;i<K;++i)     for(int j=0;j<K;++j)
//                newalpha[t][i][j] /= factor;
     //      }

  }//end t loop

  //Backward Probs
  if( CalculateBeta ){
    for( int t = Transitions-2; t >=0; t-- ){
      for(int j1 = 0; j1 < K; ++j1){
	for(int j2 =0; j2 <K; ++j2){
	  
	  Sum1 = Sum2 = Sum3 = 0.0;
	  for(int i1=0;i1<K;++i1){
	    Sum1 += lambda[t+1][i1][j2]*newbeta[t+1][i1][j2]*Admixture(i1,0);
	    Sum2 += lambda[t+1][j1][i1]*newbeta[t+1][j1][i1]*Admixture(i1,Mcol);
	    for(int i2=0;i2<K;++i2) Sum3 += lambda[t+1][i1][i2]*newbeta[t+1][i1][i2]*Admixture(i1,0)*Admixture(i2,Mcol);
	  }
	  newbeta[t][j1][j2] = p[t+1] * lambda[t+1][j1][j2] * newbeta[t+1][j1][j2]
	                     + q[t+1] * Sum1
	                     + r[t+1] * Sum2 
	                     + s[t+1] * Sum3;
	}
      }
    }//end beta calculation 
  }

}

/*
  Updates Forward and (if required) Backward probabilities
  haploid case only
  Here Admixture is a column matrix and the last dimensions of f and lambda are 1.
*/
void HMM::UpdateProbsHaploid(double *f[], int StartLocus,Matrix_d& Admixture, double ***lambda, bool CalculateBeta){

  double Sum;

  for(int j=0; j<States;++j){
    newalpha[0][j][0] = Admixture(j,0) * lambda[0][j][0];
    newbeta[Transitions-1][j][0] = 1.0;
  }

  for( int t = 1; t < Transitions; t++ ){
    Sum = 0.0;
    for(int j=0;j<States;++j){
      Sum += newalpha[t-1][j][0];
    }
    for(int j=0;j<0;++j){
      newalpha[t][j][0] = f[0][StartLocus+t] + (1.0 - f[0][StartLocus+t]) * Admixture(0,j) * Sum;
      newalpha[t][j][0] *= lambda[t+1][j][0];
    }

    if(CalculateBeta){
      for( int t = Transitions-2; t >=0; t-- ){
	Sum = 0.0;
	for(int j=0;j<States;++j){
	  Sum += Admixture(0,j)*lambda[t+1][j][0]*newbeta[t+1][j][0];
	}
	for(int j=0;j<0;++j){
	  newbeta[t][j][0] = f[StartLocus+t+1][0]*lambda[t+1][j][0]*newbeta[t+1][j][0] + (1.0 - f[0][StartLocus+t+1])*Sum;
	}
      }
      
    }
  }

}

// void HMM::GetStateProbs( double * probs, int i )
// {
//   //Vector_d probs( States );
//   double sum = 0.0;
//    for( int j = 0; j < States; j++ ){
//      probs[j] = alpha[i][j] * beta[i][j];
//      sum += probs[j];
//    }
//    //if(sum() == 0.0)cout<<alpha<<endl<<beta<<endl;
//    for( int j = 0; j < States; j++ ){
//      probs[j] /= sum;
//    }
//    //return probs;
// }

/*
  computes conditional state probabilities at "time" t
  probs - double array to store state probs
  K = #populations
*/
void HMM::NewGetStateProbs( double * probs, int t)
{
  double sum = 0.0;
  int State = 0;

  for(int i = 0; i < K; ++i)
   for( int j = 0; j < K; j++ ){
     probs[State++] = newalpha[t][i][j] * newbeta[t][i][j];
     sum += probs[State-1];
   }

   for( int j = 0; j < States; j++ ){
     probs[j] /= sum;
   }
}


// double HMM::getLikelihood()
// {
//    double sum = 0;
//    for( int j = 0; j < States; j++ ){
//      sum += alpha[Transitions - 1][j];
//    }
//    return( sumfactor+log( sum ) );
// }

/*
  returns log-likelihood
  This is the sum over states of products of alpha and beta
  and is the same for all t so it is convenient to compute for
  t = T-1 since beta here is 1. 
*/
double HMM::NewgetLikelihood()
{
   double sum = 0;
  for(int i = 0; i < K; ++i)
   for( int j = 0; j < K; j++ ){
     sum += newalpha[Transitions - 1][i][j];
   }
   return( sumfactor+log( sum ) );
}


// void HMM::Sample(int *C)
// {
//   //C - an int array to store the sampled states
//   //   Vector_d V(States);
//   double V[States];

//    for( int j = 0; j < States; j++ )V[j] = alpha[Transitions - 1][j];
//    C[ Transitions - 1 ] = SampleFromDiscrete3( V, States );
//    for( int t =  Transitions - 2; t >= 0; t-- ){
//      for(int j = 0; j < States; j++)V[j] = TransitionProbs[t][j][ C[t+1] ];
//       for( int j = 0; j < States; j++ )
// 	V[j] *= alpha[t][j];
//       C[ t ] = SampleFromDiscrete3( V, States );
//    }
// }

/*
  Samples Hidden States
  ---------------------
  C          - an int array to store the sampled states
  StartLocus - number of first locus on chromosome
  Mcol       - indicator for random mating model
  isdiploid  - indicator for diploidy
*/
void HMM::NewSample(int *C, Matrix_d &Admixture, double *f[], int StartLocus,int Mcol,bool isdiploid)
{
  int j1,j2;
  double V[States];

  if(isdiploid){
    int State = 0;
    for( int i1 = 0; i1 < K; i1++ )for(int i2 = 0; i2<K; ++i2)V[State++] = newalpha[Transitions - 1][i1][i2];
    C[ Transitions - 1 ] = SampleFromDiscrete3( V, States );
    
    for( int t =  Transitions - 2; t >= 0; t-- ){
      j1 = (int) (C[t+1]/K);//j
      j2 = C[t+1]-K*j1;     //j'

      State = 0;
      for(int i1 = 0; i1 < K; i1++)for(int i2=0;i2<K;++i2){
// 	V[State] = (i1==j1)*(i2==j2)*p[t+1] + (i2==j2)*q[t+1]*Admixture(j1,0) 
//                  + (i1==j1)*r[t+1]*Admixture(j2,Mcol) + Admixture(j1,0)*Admixture(j2,Mcol)*s[t+1];
	V[State] = (i1==j1)*(i2==j2)*p[t+1] + (i2==j2)*Q[t+1][j1] 
                 + (i1==j1)*R[t+1][j2] + S[t+1][j1][j2];
	V[State] *= newalpha[t][i1][i2];
	State++;
      }
      C[ t ] = SampleFromDiscrete3( V, States );
    }
  }
  else{//haploid
    for( int j = 0; j < States; j++ )V[j] = newalpha[Transitions - 1][j][0];
    C[ Transitions - 1 ] = SampleFromDiscrete3( V, States );
    for( int t =  Transitions - 2; t >= 0; t-- ){
      for(int j = 0; j < States; j++)V[j] = (j == C[t+1])*f[0][StartLocus+t+1]+Admixture(C[t+1],0)*(1.0 - f[0][StartLocus+t]);
      for( int j = 0; j < States; j++ )	V[j] *= newalpha[t][j][0];
      C[ t ] = SampleFromDiscrete3( V, States );
   }
  }

}

//set a transition probability
//could do with a function to set all at once
//will soon be obsolete
// void HMM::SetTProb(int i, int j, int k, double x){
//   TransitionProbs[i][j][k] = x;
// }

