#include "functions.h"
#include <cassert>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "gsl/gsl_sf_exp.h"
#include <numeric>
#include <iostream>
#include <gsl/gsl_math.h>

using namespace::std;

// ************** Log Densities *******************

double getGammaLogDensity(double alpha, double beta, double x)
{
  return alpha * log(beta) - gsl_sf_lngamma(alpha) + (alpha-1.0)*log(x) - beta*x;
}

//calls gsl function for computing the logarithm of the probability density p(x_1, ... , x_K) 
//for a Dirichlet distribution with parameters a[K]. 
//a should be of length K but x can be of length K or K-1 since last element ignored
double getDirichletLogDensity(const double* const alpha, const double* const x, size_t K)
{
  double f, xsum = 0.0, sumalpha = alpha[K-1];
  vector<double> theta(K);

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

double getDirichletLogDensity(const std::vector<double>& a, const double* const x)
{
  size_t K = a.size();
  double f, xsum = 0.0;
  vector<double> theta(K);

  for(size_t k = 0; k < K-1; ++k){
    theta[k] = x[k];
     xsum += x[k];
  }

  theta[K-1] = 1.0 - xsum;

  double sum = accumulate(a.begin(), a.end(), 0.0, std::plus<double>());//sum of a
  f = gsl_sf_lngamma( sum );
  for( unsigned i = 0; i < K; i++ )
    if( a[i] > 0.0 )
      f += ( a[i] - 1 ) * log( theta[i] ) - gsl_sf_lngamma( a[i] );

  return f;
}
double getDirichletLogDensity(const std::vector<double>& a, const std::vector<double>& x)
{
  size_t K = a.size();
  double f, xsum = 0.0;
  vector<double> theta(K);

  for(size_t k = 0; k < K-1; ++k){
    theta[k] = x[k];
     xsum += x[k];
  }

  theta[K-1] = 1.0 - xsum;

  double sum = accumulate(a.begin(), a.end(), 0.0, std::plus<double>());//sum of a
  f = gsl_sf_lngamma( sum );
  for( unsigned i = 0; i < K; i++ )
    if( a[i] > 0.0 )
      f += ( a[i] - 1 ) * log( theta[i] ) - gsl_sf_lngamma( a[i] );

  return f;
}

double getDirichletLogDensity_Softmax(const std::vector<double>& a, double *x) {
  size_t K = a.size();
  double f = 0.0;
  vector<double> theta(K);
  theta[K-1] = 1.0;
  
  for(size_t k = 0; k < K-1; ++k) {
    theta[k] = x[k];
    theta[K-1] -= x[k];
  }
  
  for( unsigned i = 0; i < K; i++ )
    if( a[i] > 0.0 ) f += ( a[i] ) * log( theta[i] );
  return f;
}

// ********** Misc. functions and transformations ******************
double AverageOfLogs(const std::vector<double>& vec, double max)
{
  double sum = 0;

  for ( unsigned int i = 0; i < vec.size(); i++ )
    sum += exp( vec[i] - max );

  sum /= vec.size();

  return log(sum) + max;
}

void inv_softmax(size_t K, const double* const mu, double *a){
  double sumlogmu = 0.0;
  for(unsigned k = 0; k < K; ++k)sumlogmu += log(mu[k]);
  double logz = -sumlogmu / (double)K;

  for(unsigned k = 0; k< K; ++k)a[k] = log(mu[k]) + logz;
}
void softmax(size_t K, double *mu, const double* a){
  //inverse of softmax transformation above
  double logz = 0.0;
  double amax = a[0];
  for(unsigned k = 1; k < K; ++k)amax = max(amax, a[k]);

  for(unsigned k = 0; k < K; ++k)
    logz += exp(a[k] - amax);
  logz = amax + log(logz);

  for(unsigned k = 0; k < K; ++k)mu[k] = exp(a[k] - logz);
}
void inv_softmax(size_t K, const double* const mu, double *a, const bool* const b){
  //b is an array of indicators
  //transformation is only applied to elements with b=true
  double sumlogmu = 0.0;
  for(unsigned k = 0; k < K; ++k)if(b[k])sumlogmu += log(mu[k]);
  double logz = -sumlogmu / (double)K;

  for(unsigned k = 0; k< K; ++k)if(b[k])a[k] = log(mu[k]) + logz;
}
void softmax(size_t K, double *mu, const double* a, const bool* const b){
  //inverse of softmax transformation
  //b is an array of indicators
  //transformation is only applied to elements with b=true
  double z = 0.0;
  for(unsigned k = 0; k < K; ++k) z += exp(a[k])*b[k];
  for(unsigned k = 0; k < K; ++k)mu[k] = b[k]*exp(a[k]) / z;
}

// ************* Matrix Algebra **************************************
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

double GaussianConditionalQuadraticForm( int kk, double *mean, double *var, size_t dim )
{
  //Note that matrix_view's do not allocate new data
  gsl_matrix *Q = gsl_matrix_alloc(1, 1);

  gsl_matrix *U1 = gsl_matrix_alloc(1, kk);
  for(int i = 0; i< kk; ++i){
    U1->data[i] = mean[i];
  }

  gsl_matrix *V11 = gsl_matrix_alloc(kk, kk);
  for(int i = 0; i < kk; ++i)
    for(int j = 0; j < kk; ++j){
      V11->data[i*kk +j] = var[i*dim +j];
    }

  if(gsl_linalg_LU_det (V11, 1)==0.0) return -1;//potentially dangerous but works


  //compute V = V11^-1 * U1
  gsl_matrix *V = gsl_matrix_alloc(kk, 1);
  gsl_vector *x = gsl_vector_alloc(kk);

  for(int j = 0; j < kk; ++j)gsl_vector_set(x, j, gsl_matrix_get(U1, 0, j));
  gsl_linalg_HH_svx(V11, x);
  
  //set column of V to x
  for(int j = 0; j < kk; ++j)gsl_matrix_set(V, j, 0, gsl_vector_get(x, j));  
  gsl_vector_free(x);

  //compute Q = U1' * V
  gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, U1, V, 0, Q);

  gsl_matrix_free(V);
  gsl_matrix_free(U1);
  gsl_matrix_free(V11);

  double result = gsl_matrix_get(Q, 0,0);
  gsl_matrix_free(Q);
  return result;
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

void submatrix(double **M, double **Sub, int r1, int r2, int c1, int c2){
  //sets Sub as submatrix of M consisting of rows r1 to r2 and cols c1 to c2
  for(int row = r1; row < r2; ++row)
    for(int col = c1; col < c2; ++col)
      Sub[row-r1][col-c1] = M[row][col];
}
void submatrix(double *M, double *Sub, int Mcols, int r1, int r2, int c1, int c2){
  for(int row = r1; row < r2; ++row)
    for(int col = c1; col < c2; ++col)
      Sub[(row-r1)*(r2-r1+1) + (col-c1)] = M[row*Mcols +col];
}
void equate_matrix(double **A, double **B, int m, int n){
  for(int i = 0; i < m; ++i)
    for(int j = 0; j < n; ++j)
      A[i][j] = B[i][j];
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
void matrix_product(double *a, double *c, size_t d1, size_t d2){
  //computes c = a * a'
  gsl_matrix_view A, At, C;
  A = gsl_matrix_view_array(a, d1, d2);
  At = gsl_matrix_view_array(a, d1, d2);
  C = gsl_matrix_view_array(c, d1, d1);

  gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, &A.matrix, &At.matrix, 0, &C.matrix); 
}
void matrix_product(const double* const a, const double* const b, double* c, size_t d1, size_t d2, size_t d3){
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

double determinant(double *a, size_t d){
  gsl_permutation *permutation = gsl_permutation_alloc(d);
  int signum;
  double* aa = new double[d*d];
  copy(a, a+d*d, aa);//make copy as LUdecomp will destroy a
  gsl_matrix_view A  = gsl_matrix_view_array(aa, d, d);

  gsl_linalg_LU_decomp( &A.matrix, permutation, &signum );//LU decomposition
  double det = gsl_linalg_LU_det(&A.matrix, signum); 

  gsl_permutation_free(permutation);
  delete[] aa;
  return det;
}

void matrix_inverse(const double* const a, double* inv, size_t d){
//inverts a matrix using LU decomposition
  gsl_permutation *permutation = gsl_permutation_alloc(d);
  int signum;
  double* aa = new double[d*d];
  copy(a, a+d*d, aa);//make copy as LUdecomp will destroy a
  gsl_matrix_view A  = gsl_matrix_view_array(aa, d, d);

  gsl_linalg_LU_decomp( &A.matrix, permutation, &signum );//LU decomposition
  gsl_matrix_view Inv = gsl_matrix_view_array(inv, d, d);
  gsl_linalg_LU_invert (&A.matrix, permutation, &Inv.matrix);

  delete[] aa;
  gsl_permutation_free(permutation);
}
//can use this version to overwrite a with its inverse
void matrix_inverse(double* a, size_t d){
//inverts a matrix using LU decomposition
  gsl_permutation *permutation = gsl_permutation_alloc(d);
  int signum;
  double* aa = new double[d*d];
  copy(a, a+d*d, aa);//make copy as LUdecomp will destroy a
  gsl_matrix_view A  = gsl_matrix_view_array(aa, d, d);

  gsl_linalg_LU_decomp( &A.matrix, permutation, &signum );//LU decomposition
  gsl_matrix_view Inv = gsl_matrix_view_array(a, d, d);
  gsl_linalg_LU_invert (&A.matrix, permutation, &Inv.matrix);

  delete[] aa;
  gsl_permutation_free(permutation);
}

//Cholesky decomposition, Crout algorithm
void cholDecomp(const double* const a, double *L, int n){
  
  for(int i = 0; i < n; ++i){
    for(int j = i+1; j < n; ++j)L[i*n +j] = 0.0;
 
    for(int j = i; j < n; ++j){
  
      double sum = a[i*n + j];
      for(int k = 0; k < i ; ++k) sum-= L[i*n +k]*L[j*n+k];
  
      if(i == j){
	if(sum <= 0.0) {cerr<<"Cholesky decomposition failed"<<endl;system("pause");exit(1);}
	L[i*n +j] = sqrt(sum);
      }
      else                          {
	L[j*n + i] = sum / L[i*n + i];
	
      }
    }
  }
}

// void invert_pds_matrix(const double* const a, double *Inv, int n){
//   //inverts a positive-definite symmetric matrix using Cholesky decomposition


// }

//useful for stl functions
double xlog(double x){
  return log(x);
}
double xexp(double x){
  double res = 0.0;
  if(x < GSL_LOG_DBL_MIN)res = exp(x);
  else res = gsl_sf_exp(x);
  return res;
}

void print_vector(std::vector<double> a){
  copy(a.begin(), a.end(), ostream_iterator<double>(cout, " "));
  cout<<endl;
}
void print_vector(std::vector<std::string> a){
  copy(a.begin(), a.end(), ostream_iterator<std::string>(cout, " "));
  cout<<endl;
}
void print_vector(std::vector<int> a){
  copy(a.begin(), a.end(), ostream_iterator<int>(cout, " "));
  cout<<endl;
}
