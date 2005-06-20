#include "functions.h"
#include <cassert>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_linalg.h>
#include <iostream>

using namespace::std;

double getGammaLogDensity(double alpha, double beta, double x)
{
  return alpha * log(beta) - gsl_sf_lngamma(alpha) + (alpha-1.0)*log(x) - beta*x;
}

double getDirichletLogDensity(const Vector_d& a, const Vector_d& x)
{
  double f;
  Vector_d xdash;    //This adds the implied last (kth) element of x
  xdash = x;                              //when x has length (k-1)
  int k = a.GetNumberOfElements();        //x must sum to 1
  if (k - 1 == x.GetNumberOfElements()) {
    xdash.AddElement( k - 1 );
    xdash( k - 1 ) = 1 - x.Sum();
  }
  assert( k == xdash.GetNumberOfElements() ); //Error in getDirichletLogDensity - lengths of vector arguments do not match.\n";

  f = gsl_sf_lngamma( a.Sum() );
  for( int i = 0; i < k; i++ )
    if( a(i) > 0.0 )
      f += ( a(i) - 1 ) * log( xdash(i) ) - gsl_sf_lngamma( a(i) );

  return f;
}

double AverageOfLogs(const std::vector<double>& vec, double max)
{
  double sum = 0;

  for ( unsigned int i = 0; i < vec.size(); i++ )
    sum += exp( vec[i] - max );

  sum /= vec.size();

  return log(sum) + max;
}

int HH_solve (Matrix_d A, Vector_d b, Vector_d *x)
{
  //Caller for gsl_linalg_HH_solve
  //This function solves the system A x = b directly using Householder transformations. 
  //On output the solution is stored in x and b is not modified. The matrix AA is destroyed by the Householder transformations. 

  gsl_matrix *AA;
  gsl_vector *bb,*xx;

  AA = gsl_matrix_calloc(A.GetNumberOfRows(),A.GetNumberOfCols());
  bb = gsl_vector_calloc(b.GetNumberOfElements());
  xx = gsl_vector_calloc(x->GetNumberOfElements());

  for (int i=0; i < A.GetNumberOfRows(); i++){
    for (int j=0; j < A.GetNumberOfCols(); j++){
      int offset = i * AA->size2 + j;
      AA->data[offset] = A(i,j);
    }
  }
  for(int i=0;i < b.GetNumberOfElements();i++)
    bb->data[i] = b(i);

  int status = gsl_linalg_HH_solve(AA,bb,xx);
  for(int i=0;i < x->GetNumberOfElements();i++)
    (*x)(i) = xx->data[i];

  gsl_matrix_free(AA);
  gsl_vector_free(bb);
  gsl_vector_free(xx);
  return status;
}

int HH_svx (Matrix_d A, Vector_d *x)
{
  //Caller for gsl_linalg_HH_svx
  // This function solves the system A x = b in-place using Householder transformations. 
  // On input x should contain the right-hand side b, which is replaced by the solution on output. 
  // The matrix AA is destroyed by the Householder transformations. 

  gsl_matrix *AA;
  gsl_vector *xx;

  AA = gsl_matrix_calloc(A.GetNumberOfRows(),A.GetNumberOfCols());
  xx = gsl_vector_calloc(x->GetNumberOfElements());

  for (int i=0; i < A.GetNumberOfRows(); i++){
    for (int j=0; j < A.GetNumberOfCols(); j++){
      int offset = i * AA->size2 + j;
      AA->data[offset] = A(i,j);
    }
  }
  for(int i=0;i < x->GetNumberOfElements();i++)
    xx->data[i] = (*x)(i);

  int status = gsl_linalg_HH_svx(AA,xx);
  for(int i=0;i < x->GetNumberOfElements();i++)
    (*x)(i) = xx->data[i];

  gsl_matrix_free(AA);
  gsl_vector_free(xx);

  return status;
}


void CentredGaussianConditional( int kk, Matrix_d mean, Matrix_d var,
					   Matrix_d *newmean, Matrix_d *newvar )
//Computes the conditional mean and variance of a centred subvector of length kk of a zero-mean Multivariate Gaussian vector
//means are matrices to allow for matrix algebra
{
  //some bounds checking would be good here
  Matrix_d mean1, mean2, Vbb, Vab, Vaa,V;
  Vector_d x;
  mean1 = mean.SubMatrix( 0, kk - 1, 0, 0 );
  mean2 = mean.SubMatrix( kk, mean.GetNumberOfRows() - 1, 0, 0 );
  Vaa = var.SubMatrix( 0, kk - 1, 0, kk - 1 );
  Vbb = var.SubMatrix( kk, var.GetNumberOfRows() - 1, kk, var.GetNumberOfCols() - 1 );
  Vab = var.SubMatrix( 0, kk - 1, kk, var.GetNumberOfCols() - 1 );

  x = mean2.GetColumn(0);
  HH_svx(Vbb, &x);
  *newmean = mean1 - Vab * (x.ColumnMatrix()); 

  V.SetNumberOfElements(Vbb.GetNumberOfRows(),Vab.GetNumberOfRows());
  x.SetNumberOfElements(Vab.GetNumberOfCols());

  for(int i =0; i<Vab.GetNumberOfRows();i++){
    x =Vab.GetRow(i);
    HH_svx(Vbb, &x);
    V.SetColumn(i,x);
  }
  *newvar = Vaa - Vab * V;
}

void CentredGaussianConditional1( Vector_d mean, Matrix_d var,
					   double *newmean, double  *newvar )
//As above but for special case of kk = 1 and using vectors and scalars 
{

  Matrix_d mu, mu2, var2;

  mu = mean.ColumnMatrix();

  CentredGaussianConditional(1, mu, var, &mu2, &var2 );
  *newmean = mu2(0,0);
  *newvar = var2(0,0); 
}

//allocate space for 2-way rectangular array of doubles and initialise to zero
double **alloc2D_d(int m,int n)
{
  double **M;
  try{
    M = new double*[m];
    if(M==NULL)throw(0);
    for(int i = 0; i < m; ++i){
      M[i] = new double[n];
      if(M[i] == NULL)throw(0);
      for(int j = 0; j < n; ++j)M[i][j] = 0.0;
    }
  }
  catch(int i){
    cout<<"Unable to allocate space for matrix"<<endl;
    exit(1);
  }
  return M;
}

//allocate space for 2-way rectangular array of ints and initialise to zero
int **alloc2D_i(int m,int n)
{
  int **M;
  try{
    M = new int*[m];
    if(M==NULL)throw(0);
    for(int i = 0; i < m; ++i){
      M[i] = new int[n];
      if(M[i] == NULL)throw(0);
      for(int j = 0; j < n; ++j)M[i][j] = 0;
    }
  }
  catch(int i){
    cout<<"Unable to allocate space for matrix"<<endl;
    exit(1);
  }
  return M;
}

//delete double matrix, even nonrectangular
void free_matrix(double **M, int m){
  if(M){
    for(int i = 0; i < m; ++i) if( M[i] )delete[] M[i];
    delete[] M;
  }
}
//delete int matrix
void free_matrix(int **M, int m){
  if(M){
    for(int i = 0; i < m; ++i) if( M[i] )delete[] M[i];
    delete[] M;
  }
}

void equate_matrix(double **A, double **B, int m, int n){
  for(int i = 0; i < m; ++i)
    for(int j = 0; j < n; ++j)
      A[i][j] = B[i][j];
}

double **MatrixAsArray(Matrix_d &M){
  double **A;
  A = alloc2D_d(M.GetNumberOfRows(), M.GetNumberOfCols());
  for(int row = 0; row < M.GetNumberOfRows(); ++row)
    for(int col = 0; col < M.GetNumberOfCols(); ++col)
      A[row][col] = M(row,col);
  return A;
}
