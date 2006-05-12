/** 
 *   ADMIXMAP
 *   functions.cc 
 *   Miscellaneous functions for admixmap, not belonging to any class
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "functions.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "gsl/gsl_sf_exp.h"
#include "gsl/gsl_sf_log.h"
#include <numeric>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <sstream>

using namespace::std;

// ************** Log Densities *******************

double getGammaLogDensity(double alpha, double beta, double x) {
  return alpha * log(beta) - gsl_sf_lngamma(alpha) + (alpha-1.0)*log(x) - beta*x;
}

double getGammaLogDensity_LogBasis(double alpha, double beta, double x) {
  return - gsl_sf_lngamma(alpha) + alpha*(log(beta) + log(x)) - beta*x;
}


double getDirichletLogDensity(const double* const alpha, const double* const x, size_t K) {
  // version with array arguments
  // calls gsl function for computing the logarithm of the probability density p(x_1, ... , x_K) 
  // for a Dirichlet distribution with parameters a[K]. 
  // a should be of length K but x can be of length K or K-1 since K th element is recalculated by subtraction
  vector<double> theta(K);
  double sumalpha = alpha[K-1];
  theta[K-1] = 1.0;
  for(size_t k = 0; k < K-1; ++k) {
    theta[k] = x[k];
    theta[K-1] -= x[k];
    sumalpha += alpha[k];
  }
  double f = gsl_sf_lngamma( sumalpha );
  for( size_t i = 0; i < K; ++i ) {
    if( alpha[i] > 0.0 ) { // to avoid bug in gsl_sf_lngamma
      f += ( alpha[i] - 1 ) * log( theta[i] ) - gsl_sf_lngamma( alpha[i] );
    }
  }
  return f;
}

double getDirichletLogDensity(const std::vector<double>& a, const double* const x) {
  // version with parameters as std vector, proportions as array 
  size_t K = a.size();
  vector<double> theta(K);
  theta[K-1] = 1.0;
  for(size_t k = 0; k < K-1; ++k) {
    theta[k] = x[k];
    theta[K-1] -= x[k];
  }
  double sum = accumulate(a.begin(), a.end(), 0.0, std::plus<double>());//sum of a
  double f = gsl_sf_lngamma( sum );
  for( unsigned i = 0; i < K; ++i )
    if( a[i] > 0.0 )
      f += ( a[i] - 1 ) * log( theta[i] ) - gsl_sf_lngamma( a[i] );
  return f;
}

double getDirichletLogDensity(const std::vector<double>& a, const std::vector<double>& x) {
  // version with both arguments as std vectors
  size_t K = a.size();
  vector<double> theta(K);
  theta[K-1] = 1.0;
  for(size_t k = 0; k < K-1; ++k) {
    theta[k] = x[k];
    theta[K-1] -= x[k];
  }
  double sum = accumulate(a.begin(), a.end(), 0.0, std::plus<double>());//sum of a
  double f = gsl_sf_lngamma( sum );
  for( unsigned i = 0; i < K; ++i ) {
    if( a[i] > 0.0 ) {
      f += ( a[i] - 1 ) * log( theta[i] ) - gsl_sf_lngamma( a[i] );
    }
  }
  return f;
}

double getDirichletLogDensity_Softmax(const std::vector<double>& a, const double* const x) {
  // version with parameters as std vector, proportions as array
  size_t K = a.size();
  vector<double> theta(K);
  theta[K-1] = 1.0;
  for(size_t k = 0; k < K-1; ++k) {
    theta[k] = x[k];
    theta[K-1] -= x[k];
  }
  double sum = accumulate(a.begin(), a.end(), 0.0, std::plus<double>());//sum of a
  double f = gsl_sf_lngamma( sum );
  for( unsigned i = 0; i < K; ++i ) {
    if( a[i] > 0.0 ) {
      f += ( a[i] ) * log( theta[i] ) - gsl_sf_lngamma( a[i] );
    }
  }
  return f;
}

void ddigam(  double *X, double *ddgam  )
{
//  FOLLOWS CLOSELY ALG AS 103 APPL.STATS.(1976) VOL 25
//  CALCS DIGAMMA(X)=D(LOG(GAMMA(X)))/DX
//
//  SET CONSTANTS.SN=NTH STIRLING COEFFICIENT,D1=DIGAMMA(1.)
//
   double S, C, S3, S4, S5, D1, Y, R;

   S = 1.0e-5;
   C = 8.5e0;
   S3 = 8.333333333e-2; 
   S4 = 8.333333333e-3;
   S5 = 3.968253968e-3;
   D1 = -0.5772156649;

//      DATA S,C,S3,S4,S5,D1/1.0D-5,8.5D0,8.333333333D-2,
//    1  8.333333333D-3,3.968253968D-3,-0.5772156649D0/


//  CHECK ARGUMENT IS POSITIVE

   *ddgam=0.0;
   Y=*X;
   if(Y < 0.0){
     throw string("Negative value passed as argument to digamma function ddigam");
   }

//  USE APPROXIMATION IF ARGUMENT .LE.S

   if(Y > S){

//  REDUCE TO DIGAMMA(X+N),(X+N).GE.C

      while( Y < C ){
         *ddgam=*ddgam-(1.0/Y);
         Y=Y+1.0;}
   
//  USE STIRLING IF ARGUMENT .GE.C

      R=1.0/Y;
      *ddgam=*ddgam+log(Y)-0.5*R;
      R=R*R;
      *ddgam=*ddgam-R*(S3-R*(S4-R*S5));}
   else
      *ddgam=D1-1.0/Y;
}

void trigam( double *x, double *trgam )
{
/*
 * closely follows alg. as 121 appl.stats. (1978) 
 * vol 27, 97-99. (b.e. schneider)
 *
 * calculates trigamma(x)=d**2(log(gamma(x)))/dx**2
 */
   double a=1.0e-4,b=5.0,one=1.0,half=0.5,y,z,trigam1=0.0;
   double b2=0.1666666667,b4=-0.03333333333;
   double b6=0.02380952381,b8=-0.0333333333;
/*
 *  b2,b4,b6,b8 are bernoulli numbers
 *
 *  check that argument is positive
 *
 */ 
   z=*x;
/*
 *  use small value approximation if x.le.a
 */
   if(z<=a){ 
      trigam1=1.0/(z*z);
      *trgam=trigam1;}
   else{
/*
 *  increase argument to (x+i).ge.b
 */
      while(z<b){
         trigam1=trigam1+1.0/(z*z);
         z=z+1.0;
      }
/*
 *  apply asymptotic formula if argument.ge.b
 */
      y=1.0/(z*z);
      trigam1=trigam1+half*y+(one+y*(b2+y*(b4+y*(b6+y*b8))))/z;
      *trgam=trigam1;
   }
}

double MultinomialPDF( const std::vector<int> r, const std::vector<double> theta )
{
  if( r.size() != theta.size() ){
    throw string("Unequal lengths of vector arguments to MultinomialPDF");
  }
  double f = 0.0;
  unsigned K = (int)r.size();
  unsigned* n = new unsigned[ K ];
  double* p = new double[ K ];
  for( unsigned i = 0; i < K; i++ ){
    p[i] = theta[i];
    n[i] = r[i];
  }
  f = gsl_ran_multinomial_pdf( K, p , n );
  delete[] n;
  delete[] p;
  return( f );
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
  // transforms proportions mu to numbers a on real line 
  // elements of a sum to 0
  double logz = 0.0;
  gsl_sf_result result;
  gsl_error_handler_t* old_handler =  gsl_set_error_handler_off();//disable default gsl error handler
  int status = 0;

  for(unsigned k = 0; k < K; ++k) {
    status = gsl_sf_log_e(mu[k], &result);
    if(status){
      stringstream err;
      err << "error in inv_softmax: " << gsl_strerror(status) << ": " << mu[k];
      throw(err.str());
    }
    a[k] = result.val;
    logz -= a[k];
  }
  logz /= (double)K;
  for(unsigned k = 0; k< K; ++k) {
    a[k] += logz;
    if( !gsl_finite(a[k]) )throw string("error in inv_softmax");
  }
  gsl_set_error_handler (old_handler);//restore gsl error handler 
}

void softmax(size_t K, double *mu, const double* a){
  //inverse of inv_softmax transformation above
  // elements of array a need not sum to zero 
  double z = 0.0;
  double amax = a[0];
  double amin= a[0];
  gsl_sf_result result;
  gsl_error_handler_t* old_handler =  gsl_set_error_handler_off();//disable default gsl error handler
  int status = 0;
  // standardize a so that max argument to exp() is 0 
  for(unsigned k = 1; k < K; ++k) {
    amax = max(amax, a[k]);
    amin = min(amin, a[k]);
  }
  for(unsigned k = 0; k < K; ++k) {
    status = gsl_sf_exp_e(a[k] + 0.5*(amax+amin), &result);
      if(status){
	stringstream s;
	s << "error in softmax: ";
	s << gsl_strerror(status)<< "\n";
	for(unsigned t = 0; t < K; ++t) s << a[t] << " ";
	throw s.str();

      }
      mu[k] = result.val;
    z += mu[k];
  }
  gsl_set_error_handler (old_handler);//restore gsl error handler 
  for(unsigned k = 0; k < K; ++k) mu[k] /= z;
}

void inv_softmax(size_t K, const double* const mu, double *a, const bool* const b){
  //transformation is applied only to elements with b=true
  double logz = 0.0;
  for(unsigned k = 0; k < K; ++k) {
    if(b[k]) {
      a[k] = log(mu[k]);
      logz -= a[k];
    }
  }
  logz /= (double)K;
  for(unsigned k = 0; k< K; ++k) if(b[k]) a[k] += logz;
}

void softmax(size_t K, double *mu, const double* a, const bool* const b){
  //transformation is applied only to elements with b=true
  double z = 0.0;
  double amax = a[0];
  // standardize a so that max argument to exp() is 0 
  for(unsigned k = 1; k < K; ++k) if(b[k]) amax = max(amax, a[k]);
  for(unsigned k = 0; k < K; ++k) {
    if(b[k]) {
      mu[k] = exp(a[k] - amax);
      z += mu[k];
    }
  }
  for(unsigned k = 0; k < K; ++k) if(b[k]) mu[k] /= z;
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

  gsl_error_handler_t* old_handler =  gsl_set_error_handler_off();//disable default gsl error handler
  int status = gsl_linalg_HH_solve(AA, &bb.vector, &xx.vector);

  //clean up  
  gsl_set_error_handler (old_handler);//restore gsl error handler 
  gsl_matrix_free(AA);

  //check for success
  if(status){
    std::string errstring = "HH_solve failed, "; errstring.append(gsl_strerror(status));
    throw(errstring) ;//throw error message up and let caller decide what to do
  }
  return 0;//will only get here if successful
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
  gsl_error_handler_t* old_handler =  gsl_set_error_handler_off();//disable default gsl error handler
  int status = gsl_linalg_HH_svx(AA, &xx.vector);

  gsl_set_error_handler (old_handler);//restore gsl error handler 
  gsl_matrix_free(AA);
  //check for success
  if(status){
    std::string errstring = "HH_svx failed, "; errstring.append(gsl_strerror(status));
    throw(errstring) ;//throw error message up and let caller decide what to do
  }
  return 0;//will only get here if successful
}

void CentredGaussianConditional( size_t kk, double *mean, double *var,
				 double *newmean, double *newvar, size_t dim )
//Computes the conditional mean and variance of a centred subvector of length kk of a zero-mean Multivariate Gaussian vector
//of length dim
{
  if(dim == (kk+1)){
    CentredGaussianConditional(mean, var, newmean, newvar, dim);
    return;
  }
  int status = 0;
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
  for(size_t i = 0; i< kk; ++i){
    mean1->data[i] = mean[i];
    newmean[i] = mean[i];     //copy mean1 into new mean
  }
  gsl_matrix *mean2 = gsl_matrix_alloc(dim-kk, 1);
  for(size_t i = 0; i< dim-kk; ++i) mean2->data[i] = mean[i+kk];

  gsl_matrix *Vaa = gsl_matrix_alloc(kk, kk);
  for(size_t i = 0; i < kk; ++i)
    for(size_t j = 0; j < kk; ++j){
      Vaa->data[i*kk +j] = var[i*dim +j];
      newvar[i*kk +j] = var[i*dim +j];  //copy Vaa into newvar
    }

  gsl_matrix *Vbb = gsl_matrix_alloc(dim-kk, dim-kk);
  for(size_t i = 0; i < dim-kk; ++i)
    for(size_t j = 0; j < dim-kk; ++j)Vbb->data[i*(dim-kk) +j] = var[(i+kk)*dim + j+kk];

  gsl_matrix *Vab = gsl_matrix_alloc(kk, dim-kk);
  for(size_t i = 0; i < kk; ++i)
    for(size_t j = 0; j < dim-kk; ++j)Vab->data[i*(dim-kk) +j] = var[i*dim +j+kk];

  //compute inv(Vbb) * mean2 
  //cannot call gsl function directly as it would destroy Vbb
  gsl_error_handler_t* old_handler =  gsl_set_error_handler_off();//disable default gsl error handler
  status = HH_svx(dim-kk, Vbb->data, mean2->data);
  //mean2 now holds the solution

  if(!status){
    //compute new mean as mean1 - Vab * Vbb^-1 * mean2 = mean1 - Vab * mean2
    gsl_matrix *C = gsl_matrix_alloc(kk,1);
    gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, Vab, mean2, 0, C);
    gsl_matrix_sub(&newmean_view.matrix, C);
    gsl_matrix_free(C);
    
    //compute new var
    gsl_matrix *V = gsl_matrix_alloc(dim-kk, kk);
    double *x = new double[dim-kk];
    
    //compute V = Vbb * tr(Vab), column by column
    for(size_t i = 0; i< kk; i++){
      //copy column of Vba (=row of Vab)into x since we still need Vab and x will be overwritten
      for(size_t j = 0; j < dim-kk; ++j)x[j] = gsl_matrix_get(Vab, i, j);
      status = HH_svx(dim-kk, Vbb->data, x);//cannot call gsl function directly as it would destroy Vbb
      
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
  }
  gsl_set_error_handler (old_handler);//restore gsl error handler 
  gsl_matrix_free(mean1);
  gsl_matrix_free(mean2);
  gsl_matrix_free(Vaa);
  gsl_matrix_free(Vbb);
  gsl_matrix_free(Vab);
  if(status){
    std::string errstring = "CentredGaussianConditional failed, "; errstring.append(gsl_strerror(status));
    throw(errstring) ;//throw error message up and let caller decide what to do
  }

}

void CentredGaussianConditional( double *mean, double *var,
				 double *newmean, double *newvar, size_t dim ){
  //special case of above with kk=dim-1
  for(unsigned i = 0; i < dim-1; ++i){
    newmean[i] = mean[i] - var[i*dim + dim-1] * mean[dim-1]/var[dim*dim-1];
    for(unsigned j = 0; j < dim-1; ++j)
      newvar[i*(dim-1)+j] = var[i*dim+j] - var[i*dim + dim-1]*var[(dim-1)*dim + j] / var[dim*dim-1];
  }
}
double GaussianMarginalQuadraticForm( int kk, double *mean, double *var, size_t dim )
{
  //returns the quadratic form in the density of the marginal distribution of the subvector of length kk of a zero-mean Gaussian vector of length dim. Useful for score tests. 
  int status = 0;
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

  //if(gsl_linalg_LU_det (V11, 1)==0.0) throw;//matrix is rank-deficient


  //compute V = V11^-1 * U1
  gsl_matrix *V = gsl_matrix_alloc(kk, 1);
  gsl_vector *x = gsl_vector_alloc(kk);

  for(int j = 0; j < kk; ++j)gsl_vector_set(x, j, gsl_matrix_get(U1, 0, j));
  gsl_error_handler_t* old_handler =  gsl_set_error_handler_off();//disable default gsl error handler
  status = gsl_linalg_HH_svx(V11, x);

  if(status){
    string error_string = "Error in HH_svx in GaussianConditionalQuadraticForm:\n";
    error_string.append(gsl_strerror(status));
    throw(error_string);
  }

  //set column of V to x
  for(int j = 0; j < kk; ++j)gsl_matrix_set(V, j, 0, gsl_vector_get(x, j));  
  gsl_vector_free(x);

  //compute Q = U1' * V
  status = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, U1, V, 0, Q);

  gsl_set_error_handler (old_handler);//restore gsl error handler 
  if(status){
    string error_string = "Error with matrix multiplication in GaussianConditionalQuadraticForm:\n";
    error_string.append(gsl_strerror(status));
    throw(error_string);
  }

  gsl_matrix_free(V);
  gsl_matrix_free(U1);
  gsl_matrix_free(V11);

  double result = gsl_matrix_get(Q, 0,0);
  gsl_matrix_free(Q);
  return result;
}

double GaussianQuadraticForm(double* mean, double* var, unsigned dim){
  double* VinvU = new double[dim];
  HH_solve(dim, mean, var, VinvU);
  double result = 0.0;
  for(unsigned i = 0; i < dim; ++i){
    result += mean[i] * VinvU[i];
  }
  delete[] VinvU;
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
    string s = "Unable to allocate space for matrix";
    throw(s);
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
    string s = "Unable to allocate space for matrix";
    throw(s);
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
    string s = "Unable to delete matrix";
    throw(s);
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
    string s = "Unable to delete matrix";
    throw(s);
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
  double det = 0.0;
  gsl_permutation *permutation = gsl_permutation_alloc(d);
  int signum;
  double* aa = new double[d*d];
  copy(a, a+d*d, aa);//make copy as LUdecomp will destroy a
  gsl_matrix_view A  = gsl_matrix_view_array(aa, d, d);

  gsl_error_handler_t* old_handler =  gsl_set_error_handler_off();//disable default gsl error handler
  int status = gsl_linalg_LU_decomp( &A.matrix, permutation, &signum );//LU decomposition
  gsl_permutation_free(permutation);

  if(status){
    delete[] aa;
    std::string errstring = "failed to compute determinant, "; errstring.append(gsl_strerror(status));
    throw(errstring) ;//throw error message up and let caller decide what to do
  }

  gsl_set_error_handler (old_handler);//restore gsl error handler 
  det = gsl_linalg_LU_det(&A.matrix, signum); 
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

  gsl_error_handler_t* old_handler =  gsl_set_error_handler_off();//disable default gsl error handler
  int status = gsl_linalg_LU_decomp( &A.matrix, permutation, &signum );//LU decomposition
  if(!status){
    gsl_matrix_view Inv = gsl_matrix_view_array(inv, d, d);
    status = gsl_linalg_LU_invert (&A.matrix, permutation, &Inv.matrix);
  }

  delete[] aa;
  gsl_set_error_handler (old_handler);//restore gsl error handler 
  gsl_permutation_free(permutation);
  if(status){
    std::string errstring = "failed to compute matrix inverse, "; errstring.append(gsl_strerror(status));
    throw(errstring) ;//throw error message up and let caller decide what to do
  }
}
//can use this version to overwrite a with its inverse
void matrix_inverse(double* a, size_t d){
//inverts a matrix using LU decomposition
  gsl_permutation *permutation = gsl_permutation_alloc(d);
  int signum;
  double* aa = new double[d*d];
  copy(a, a+d*d, aa);//make copy as LUdecomp will destroy a
  gsl_matrix_view A  = gsl_matrix_view_array(aa, d, d);

  gsl_error_handler_t* old_handler =  gsl_set_error_handler_off();//disable default gsl error handler
  int status =  gsl_linalg_LU_decomp( &A.matrix, permutation, &signum );//LU decomposition
  if(!status){
    gsl_matrix_view Inv = gsl_matrix_view_array(a, d, d);
   status =  gsl_linalg_LU_invert (&A.matrix, permutation, &Inv.matrix);
  }

  delete[] aa;
  gsl_set_error_handler (old_handler);//restore gsl error handler 
  gsl_permutation_free(permutation);
  if(status){
    std::string errstring = "failed to compute matrix inverse, "; errstring.append(gsl_strerror(status));
    throw(errstring) ;//throw error message up and let caller decide what to do
  }
}

//Cholesky decomposition, Crout algorithm
int cholDecomp(const double* const a, double *L, int n){
  
  for(int i = 0; i < n; ++i){
    for(int j = i+1; j < n; ++j)L[i*n +j] = 0.0;
 
    for(int j = i; j < n; ++j){
  
      double sum = a[i*n + j];
      for(int k = 0; k < i ; ++k) sum-= L[i*n +k]*L[j*n+k];
  
      if(i == j){
	if(sum <= 0.0) {return 1;}
	L[i*n +j] = sqrt(sum);
      }
      else                          {
	L[j*n + i] = sum / L[i*n + i];
	
      }
    }
  }
  return 0;
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
double mylog(double x){
  //log function with error handling
  gsl_error_handler_t* old_handler =  gsl_set_error_handler_off();//disable default gsl error handler
  gsl_sf_result result;
  int status = gsl_sf_log_e(x, &result);
  gsl_set_error_handler (old_handler);//restore gsl error handler 
  if(status){
    stringstream s;
    s << "Error in log(" << x << "): "<< gsl_strerror(status);
    throw s.str();
  }
  return result.val;
}
double myexp(double x){
  //exp function with error handling
  gsl_error_handler_t* old_handler =  gsl_set_error_handler_off();//disable default gsl error handler
  gsl_sf_result result;
  int status = gsl_sf_exp_e(x, &result);
  gsl_set_error_handler (old_handler);//restore gsl error handler 
  if(status){
    stringstream s;
    s << "Error in exp(" << x << "): "<< gsl_strerror(status);
    throw s.str();
  }
  return result.val;

}
