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

HMM::HMM( int inTransitions, int inStates )
{
   States = inStates;
   Transitions = inTransitions;
   // forward probabilities
   alpha.SetNumberOfElements( Transitions );
   // backward probabilities
   beta.SetNumberOfElementsWithDimensions( Transitions, States, 1 ); 
}

HMM::~HMM()
{
}

void HMM::SetDimensions( int inTransitions, int inStates )
{
   States = inStates;
   Transitions = inTransitions;
   // forward probabilities
   alpha.SetNumberOfElements( Transitions );
   // backward probabilities
   beta.SetNumberOfElementsWithDimensions( Transitions, States, 1 ); 
}

void HMM::
//UpdateParameters( Matrix_d &inStationaryDist, MatrixArray_d &inTransitionProbs, MatrixArray_d &inLikelihood, bool CalculateBeta )
// name of second parameter changed so no need to allocate new memory for a copy of TransitionProbs 
UpdateParameters( Matrix_d &inStationaryDist, MatrixArray_d &TransitionProbs, MatrixArray_d &inLikelihood, bool CalculateBeta )
{
  StationaryDist = inStationaryDist;
  //   TransitionProbs = inTransitionProbs;
  Likelihood = inLikelihood;
  if( Transitions > 1 )
    CheckArguments(TransitionProbs);
  _CalculateBeta = CalculateBeta;
  UpdateFwrdBckwdProbabilities(TransitionProbs);
}

void HMM::UpdateFwrdBckwdProbabilities(MatrixArray_d &TransitionProbs)
{
// Matrix multiplication replaced by component-wise multiplication.
// ie. wasteful to multiply by diagonal matrix, multiply each column
// in turn by a scalar.
   factor = 0.0;
   alpha(0) = StationaryDist.Transpose();
   for( int j = 0; j < States; j++ )
      alpha(0)( 0, j ) *= Likelihood(0)( j, 0 );
   for( int i = 1; i < Transitions; i++ ){
      alpha(i) = alpha( i - 1 ) * TransitionProbs( i - 1 );
      for( int j = 0; j < States; j++ )
         alpha(i)( 0, j ) *= Likelihood(i)( j, 0 );
//       if( !(i%5) ){
//          factor += log(alpha(i-1).GetRow(0).Sum());
//          alpha(i) /= alpha(i-1).GetRow(0).Sum();
//      }
   }

   if( _CalculateBeta ){
      beta( Transitions - 1 ).SetNumberOfElements( States, 1 );
      beta( Transitions - 1 ).SetElements( 1.0 );
      for( int i = Transitions - 2; i >= 0; i-- ){
         for( int j = 0; j < States; j++ )
            beta(i)( 0, j ) = Likelihood( i + 1 )( j, 0 ) * beta( i + 1 )( j, 0 );
         beta(i) = TransitionProbs(i) * beta(i);
      }
   }
}

Vector_d HMM::GetStateProbs( int i )
{
   Vector_d probs( States );
   for( int j = 0; j < States; j++ ){
      probs(j) = alpha(i)( 0, j ) * beta(i)( j, 0 );
   }
   if(probs.Sum() == 0.0)cout<<alpha<<endl<<beta<<endl;
   probs /= probs.Sum();
   return probs;
}

double HMM::getLikelihood()
{
   double sum = 0;
   for( int j = 0; j < States; j++ ){
      sum += alpha( Transitions - 1 )(0,j);
   }
   return( factor+log( sum ) );
}

Vector_i HMM::Sample(MatrixArray_d &TransitionProbs)
{
   Vector_i C( Transitions );
   Vector_d V;

   V = alpha( Transitions - 1 ).GetRow(0);
   C( Transitions - 1 ) = SampleFromDiscrete2( &V );
   for( int i =  Transitions - 2; i >= 0; i-- ){
      V = TransitionProbs( i ).GetColumn( C( i + 1 ) );
      for( int j = 0; j < States; j++ )
         V(j) *= alpha( i )( 0, j );
      C( i ) = SampleFromDiscrete2( &V );
   }
   return( C );
}

void HMM::CheckArguments(MatrixArray_d &TransitionProbs)
{
   assert( Likelihood(0).GetNumberOfRows() == States );
   assert( TransitionProbs.GetNumberOfElements() == Transitions - 1 );
   assert( TransitionProbs(0).GetNumberOfRows() == States );
}
