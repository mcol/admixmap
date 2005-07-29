#include "functions.h"
#include <cassert>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <numeric>
#include <iostream>

using namespace::std;

double getGammaLogDensity(double alpha, double beta, double x)
{
  return alpha * log(beta) - gsl_sf_lngamma(alpha) + (alpha-1.0)*log(x) - beta*x;
}

double getDirichletLogDensity(const Vector_d& a, const Vector_d& x)
{
//   double f;
//   Vector_d xdash;    //This adds the implied last (kth) element of x
//   xdash = x;                              //when x has length (k-1)
//   int k = a.GetNumberOfElements();        //x must sum to 1
//   if (k - 1 == x.GetNumberOfElements()) {
//     xdash.AddElement( k - 1 );
//     xdash( k - 1 ) = 1 - x.Sum();
//   }
//   assert( k == xdash.GetNumberOfElements() ); //Error in getDirichletLogDensity - lengths of vector arguments do not match.\n";

//   f = gsl_sf_lngamma( a.Sum() );
//   for( int i = 0; i < k; i++ )
//     if( a(i) > 0.0 )
//       f += ( a(i) - 1 ) * log( xdash(i) ) - gsl_sf_lngamma( a(i) );

  size_t K = a.GetNumberOfElements();
  double f, xsum = 0.0;
  double theta[K];

  for(size_t k = 0; k < K-1; ++k){
    theta[k] = x(k);
     xsum += x(k);
  }

  theta[K-1] = 1.0 - xsum;


  f = gsl_sf_lngamma( a.Sum() );
  for( unsigned i = 0; i < K; i++ )
    if( a(i) > 0.0 )
      f += ( a(i) - 1 ) * log( theta[i] ) - gsl_sf_lngamma( a(i) );

  return f;
}

double getDirichletLogDensity(const std::vector<double>& a, const Vector_d& x)
{
  size_t K = a.size();
  double f, xsum = 0.0;
  double theta[K];

  for(size_t k = 0; k < K-1; ++k){
    theta[k] = x(k);
     xsum += x(k);
  }

  theta[K-1] = 1.0 - xsum;

  double sum = accumulate(a.begin(), a.end(), 0.0, std::plus<double>());//sum of a
  f = gsl_sf_lngamma( sum );
  for( unsigned i = 0; i < K; i++ )
    if( a[i] > 0.0 )
      f += ( a[i] - 1 ) * log( theta[i] ) - gsl_sf_lngamma( a[i] );

  return f;
}

//calls gsl function for computing the logarithm of the probability density p(x_1, ... , x_K) 
//for a Dirichlet distribution with parameters a[K]. 
//a should be of length K but x can be of length K or K-1 since last element ignored
double getDirichletLogDensity(double *alpha, double *x, size_t K)
{
  double f, xsum = 0.0, sumalpha = alpha[K-1];
  double theta[K];

  for(size_t k = 0; k < K-1; ++k){
    theta[k] = x[k];
    xsum += x[k];
    sumalpha += alpha[k];
  }
  theta[K-1] = 1.0 - xsum;

  f = gsl_sf_lngamma( sumalpha );
  for( size_t i = 0; i < K; i++ )
    if( alpha[i] > 0.0 )// to avoid bug in gsl_sf_lngamma
      f += ( alpha[i] - 1 ) * log( theta[i] ) - gsl_sf_lngamma( alpha[i] );

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

int HH_solve (size_t n, double *A, double *b, double *x)
{
  //Caller for gsl_linalg_HH_solve
  //This function solves the system A x = b directly using Householder transformations. 
  //On output the solution is stored in x and b is not modified. The matrix AA is destroyed by the Householder transformations. 

  gsl_matrix *AA;
  gsl_vector_view bb,xx;

  //create copy of A as gsl matrix; might not be necessary, depending on use
  AA = gsl_matrix_calloc(n,n);
  for (size_t i = 0; i < n*n; i++){
      AA->data[i] = A[i];
    }

  bb = gsl_vector_view_array(b, n);
  xx = gsl_vector_view_array(x, n);

  int status = gsl_linalg_HH_solve(AA, &bb.vector, &xx.vector);

  gsl_matrix_free(AA);
  return status;
}

int HH_svx (size_t n, double *A, double *x)
{
  //Caller for gsl_linalg_HH_svx
  // This function solves the system A x = b in-place using Householder transformations. 
  // On input x should contain the right-hand side b, which is replaced by the solution on output. 
  // The matrix AA is destroyed by the Householder transformations. 

  gsl_matrix *AA; //cannot use a matrix view because of above
  gsl_vector_view xx = gsl_vector_view_array(x,n); 

  //create copy of A as gsl matrix
  AA = gsl_matrix_calloc(n,n);
  for (size_t i = 0; i < n*n; i++){
      AA->data[i] = A[i];
    }

  int status = gsl_linalg_HH_svx(AA, &xx.vector);

  gsl_matrix_free(AA);
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

void CentredGaussianConditional( int kk, double *mean, double *var,
				 double *newmean, double *newvar, size_t dim )
//Computes the conditional mean and variance of a centred subvector of length kk of a zero-mean Multivariate Gaussian vector
//of length dim
{
  //Note that matrix_view's do not allocate new data
  gsl_matrix_view mean_matrix = gsl_matrix_view_array(mean, dim, 1);
  gsl_matrix_view var_matrix = gsl_matrix_view_array(var, dim, dim);

  gsl_matrix_view newmean_view = gsl_matrix_view_array(newmean, kk, 1);
  gsl_matrix_view newvar_view = gsl_matrix_view_array(newvar, kk, kk);

  /*create views of mean and var as:
    mean = [mean1 | mean2]
    var  = [Vaa   | Vab ] kk rows
           [------|-----]                 
           [      | Vbb ]
           <- kk ->
  */
  gsl_matrix *mean1 = gsl_matrix_alloc(kk, 1);
  for(int i = 0; i< kk; ++i){
    mean1->data[i] = mean[i];
    newmean[i] = mean[i];     //copy mean1 into new mean
  }
  gsl_matrix *mean2 = gsl_matrix_alloc(dim-kk, 1);
  for(size_t i = 0; i< dim-kk; ++i) mean2->data[i] = mean[i+kk];

  gsl_matrix *Vaa = gsl_matrix_alloc(kk, kk);
  for(int i = 0; i < kk; ++i)
    for(int j = 0; j < kk; ++j){
      Vaa->data[i*kk +j] = var[i*dim +j];
      newvar[i*kk +j] = var[i*dim +j];  //copy Vaa into newvar
    }

  gsl_matrix *Vbb = gsl_matrix_alloc(dim-kk, dim-kk);
  for(size_t i = 0; i < dim-kk; ++i)
    for(size_t j = 0; j < dim-kk; ++j)Vbb->data[i*(dim-kk) +j] = var[(i+kk)*dim + j+kk];

  gsl_matrix *Vab = gsl_matrix_alloc(kk, dim-kk);
  for(int i = 0; i < kk; ++i)
    for(size_t j = 0; j < dim-kk; ++j)Vab->data[i*(dim-kk) +j] = var[i*dim +j+kk];

  //compute inv(Vbb) * mean2 
  //cannot call gsl function directly as it would destroy Vbb
  HH_svx(dim-kk, Vbb->data, mean2->data);
  //mean2 now holds the solution

  //compute new mean as mean1 - Vab * Vbb^-1 * mean2 = mean1 - Vab * mean2
  gsl_matrix *C = gsl_matrix_alloc(kk,1);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, Vab, mean2, 0, C);
  gsl_matrix_sub(&newmean_view.matrix, C);
  gsl_matrix_free(C);

  //compute new var
  gsl_matrix *V = gsl_matrix_alloc(dim-kk, kk);
  double *x = new double[dim-kk];

  //compute V = Vbb * tr(Vab), column by column
  for(int i = 0; i< kk; i++){
    //copy column of Vba (=row of Vab)into x since we still need Vab and x will be overwritten
    for(size_t j = 0; j < dim-kk; ++j)x[j] = gsl_matrix_get(Vab, i, j);
    HH_svx(dim-kk, Vbb->data, x);//cannot call gsl function directly as it would destroy Vbb

    //set column of V to x
    for(size_t j = 0; j < dim-kk; ++j)gsl_matrix_set(V, j, i, x[j]);//V->data[j*kk +i] = x[j];
  }
  delete[] x;
  //now compute newvar = Vaa - Vab * V

  //compute Vab * V
  gsl_matrix *D = gsl_matrix_alloc(kk,kk);
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, Vab, V, 0, D);
  gsl_matrix_free(V);
  //subtract from newvar
  gsl_matrix_sub(&newvar_view.matrix, D);
  gsl_matrix_free(D);
  gsl_matrix_free(mean1);
  gsl_matrix_free(mean2);
  gsl_matrix_free(Vaa);
  gsl_matrix_free(Vbb);
  gsl_matrix_free(Vab);

}

//allocate space for 2-way rectangular array of doubles and initialise to zero
double **alloc2D_d(int m, int n)
{
  double **M = 0;
  try{
    M = new double*[m];
    if(M==NULL)throw(0);
    else for(int i = 0; i < m; ++i)M[i] = NULL;
    for(int i = 0; i < m; ++i){
      M[i] = new double[n];
      if(M[i] == NULL)throw(0);
      else for(int j = 0; j < n; ++j)M[i][j] = 0.0;
    }
  }
  catch(int i){
    cout<<"Unable to allocate space for matrix"<<endl;
    exit(1);
  }
  return M;
}

//allocate space for 2-way rectangular array of ints and initialise to zero
int **alloc2D_i(int m, int n)
{
  int **M = 0;
  try{
    M = new int*[m];
    if(M==NULL)throw(0);
    else for(int i = 0; i < m; ++i)M[i] = NULL;
    for(int i = 0; i < m; ++i){
      M[i] = new int[n];
      if(M[i] == NULL)throw(0);
      else for(int j = 0; j < n; ++j)M[i][j] = 0;
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
  try{
    if(M){
      for(unsigned i = (unsigned)m; i > 0 ; --i) 
	if(M[i-1])delete[] M[i-1];
      delete[] M;
    }
  }
  catch(...){
    cout<<"unable to delete matrix"<<endl;
    exit(1);
  }
}
//delete int matrix
void free_matrix(int **M, int m){
  try{
    if(M){
      for(unsigned i = (unsigned)m; i > 0; --i) if( M[i-1] )delete[] M[i-1];
      delete[] M;
    }
  }
  catch(...){
    cout<<"unable to delete matrix"<<endl;
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

//adds matrix b to matrix a
//dimensions must both be (d1 x d2)
void add_matrix(double *a, double *b, size_t d1, size_t d2){
  gsl_matrix_view A, B;
  A = gsl_matrix_view_array(a, d1, d2);
  B = gsl_matrix_view_array(b, d1, d2);
  gsl_matrix_add(&A.matrix, &B.matrix);
}


//computes c = a * b, where all args are arrays representing matrices
// and a is (d1 x d2) and b is (d2 x d3)
//
void matrix_product(double *a, double *b, double *c, size_t d1, size_t d2, size_t d3){
  gsl_matrix_view A, B, C;
  A = gsl_matrix_view_array(a, d1, d2);
  B = gsl_matrix_view_array(b, d2, d3);
  C = gsl_matrix_view_array(c, d1, d3);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, &A.matrix, &B.matrix, 0, &C.matrix); 
}
void matrix_product(const double *a, const double *b, double *c, size_t d1, size_t d2, size_t d3){
  gsl_matrix_view A, B, C;
  A = gsl_matrix_view_array(const_cast<double *>(a), d1, d2);
  B = gsl_matrix_view_array(const_cast<double *>(b), d2, d3);
  C = gsl_matrix_view_array(c, d1, d3);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, &A.matrix, &B.matrix, 0, &C.matrix); 
}

//multiplies (d1 x d2) matrix a by c
void scale_matrix(double *a, const double c, size_t d1, size_t d2){
  gsl_matrix_view A = gsl_matrix_view_array(a, d1, d2);
  gsl_matrix_scale(&A.matrix, c);
}
