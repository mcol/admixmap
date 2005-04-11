#include "rand.h"
#include "vector_i.h"
#include "vector_d.h"
#include "matrix_d.h"
#include "MatrixArray_d.h"
#include "HMM.h"

using namespace std;

HMM::HMM()
{
}

//not currently used
//HMM objects are instantiated in Chromosome using default constructor above
//and dimensions set by SetDimensions below
HMM::HMM( int inTransitions, int inStates )
{
  States = inStates;
  Transitions = inTransitions;
  TransitionProbs = new double**[ Transitions - 1];
  for(int i=0;i< Transitions - 1; i++){
    TransitionProbs[i] = new double *[States];
    for(int j=0;j < States; j++) TransitionProbs[i][j] = new double[States];
  }

  // forward probabilities
  alpha = new double*[Transitions];
  // backward probabilities
  beta = new double*[Transitions];
  betarow = new double[States];

  for(int i=0;i < Transitions;++i) {
    alpha[i] = new double[States];
    beta[i] = new double[States];
  }
}

HMM::~HMM()
{
  delete[] alpha;
  delete[] beta;
  delete betarow;
  delete[] TransitionProbs;
}

void HMM::SetDimensions( int inTransitions, int inStates )
{
  States = inStates;
  Transitions = inTransitions;
  TransitionProbs = new double**[ Transitions - 1];
  for(int i=0;i< Transitions - 1; i++){
    TransitionProbs[i] = new double *[States];
    for(int j=0;j < States; j++) TransitionProbs[i][j] = new double[States];
  }
  // forward probabilities
  alpha = new double*[Transitions];
  // backward probabilities
  beta = new double*[Transitions];
  betarow = new double[States];

  for(int i=0;i < Transitions;++i) {
    alpha[i] = new double[States];
    beta[i] = new double[States];
  }
}

//Update Forward and (optionally) Backward probs
void HMM::UpdateFwrdBckwdProbabilities( Vector_d &StationaryDist, MatrixArray_d &Likelihood, bool CalculateBeta)
{
  //StationaryDist = delta, vector of length States 
  //TransitionProbs = Gamma, (States x States) Matrix
  //Likelihood = lambda = diag (States x States) matrix of p(Hidden State | Observed State)
  // CalculateBeta indicates whether to calculate backward probs beta. They are only needed for state probs.

  //Check dimensions of arguments
  if( Transitions > 1 )  if(Likelihood(0).GetNumberOfRows() != States){
    //need to throw an exception here
    cout << "Error in Likelihood passed to HMM::Update"<<endl;
    exit(1);
  }

   //set alpha(0) = StationaryDist * lambda(0) 
   for( int j = 0; j < States; j++ )
     alpha[0][j] =  StationaryDist(j) *Likelihood(0)( j, 0);

   for( int t = 1; t < Transitions; t++ ){
     //set alpha(t) = alpha(t-1) * Gamma * lambda(t)
     for( int i = 0; i < States; i++ ){
       alpha[t][i] = 0.0;
       //explicit element-wise vector-matrix multiplication
       for( int j = 0; j < States; j++ )alpha[t][i] += alpha[t-1][j] * TransitionProbs[t-1][j][i];
       alpha[t][i] *= Likelihood(t)( i, 0);
     }

     //This is to avoid underflow
     //probably needs adjusting
     //factor = 0.0; 
     //       if( !(t%5) ){
     //          factor += log(alpha(t-1).GetRow(0).Sum());
     //          alpha(t) /= alpha(t-1).GetRow(0).Sum();
     //      }
   }

   if( CalculateBeta ){
     //beta( Transitions - 1 ).SetNumberOfElements( States, 1 );
     //beta( Transitions - 1 ).SetElements( 1.0 );
//       for( int t = Transitions - 2; t >= 0; t-- ){
//          for( int j = 0; j < States; j++ )
// 	   beta(t)( 0, j ) = Likelihood( t + 1)(j, 0) * beta( t + 1 )( j, 0 );
//          beta(t) = TransitionProbs(t) * beta(t);
//       }

     //set beta(T) = 1
     for(int i = 0; i < States; ++i)beta[Transitions - 1][i] = 1.0;
      for( int t = Transitions - 2; t >= 0; t-- ){

	//set beta(t) (temporarily) to lambda(t) * beta(t+1) 	
	for( int i = 0; i < States; i++ )
	  beta[t][i] = Likelihood( t + 1)(i, 0) * beta[ t + 1 ][i];

	for( int i = 0; i < States; i++ ){
	  //premultiply beta(t) by Gamma 
	  betarow[i] = 0.0;
	//element-wise matrix-vector multiplication
	  for( int j = 0; j < States; j++ )betarow[i] += beta[t][j] * TransitionProbs[t][i][j];
	  beta[t][i] *= betarow[i];
	}	
      }
   }
   
}

Vector_d HMM::GetStateProbs( int i )
{
   Vector_d probs( States );
   for( int j = 0; j < States; j++ ){
     probs(j) = alpha[i][j] * beta[i][j];
   }
   if(probs.Sum() == 0.0)cout<<alpha<<endl<<beta<<endl;
   probs /= probs.Sum();
   return probs;
}

double HMM::getLikelihood()
{
   double sum = 0;
   for( int j = 0; j < States; j++ ){
     sum += alpha[Transitions - 1][j];
   }
   return( factor+log( sum ) );
}

//Sample Hidden States
Vector_i HMM::Sample()
{
   Vector_i C( Transitions );
   Vector_d V(States);

   for( int j = 0; j < States; j++ )V(j) = alpha[Transitions - 1][j];
   C( Transitions - 1 ) = SampleFromDiscrete2( &V );
   for( int t =  Transitions - 2; t >= 0; t-- ){
     //V = TransitionProbs[ t ].GetColumn( C( t + 1 ) );
     for(int j = 0; j < States; j++)V(j) = TransitionProbs[t][j][ C(t+1) ];
      for( int j = 0; j < States; j++ )
	V(j) *= alpha[t][j];
      C( t ) = SampleFromDiscrete2( &V );
   }
   return( C );
}

//set a transition probability
//could do with a function to set all at once
void HMM::SetTProb(int i, int j, int k, double x){
  TransitionProbs[i][j][k] = x;
}

