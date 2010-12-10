//=============================================================================
//
// Copyright (C) 2002-2007  David O'Donnell and Paul McKeigue
// Portions Copyright (C) 2010  Marco Colombo
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file functions.cc
/// Miscellaneous functions not belonging to any class.
//=============================================================================

#include "bclib/dist.h"
#include "bclib/linalg.h"
#include "bclib/misc.h"
#include "bclib/pvector.h"
#include "bclib/Exceptions.h"
#include "bclib/GSLErrorHandler.h"
#include "bclib/GSLExceptions.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_exp.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_log.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>
#include <sstream>

using namespace::std;

BEGIN_BCLIB_NAMESPACE

// ************** Log Densities *******************

double getGammaLogDensity(double alpha, double beta, double x) {
  return alpha * log(beta) - gsl_sf_lngamma(alpha) + (alpha-1.0)*log(x) - beta*x;
}

double getGammaLogDensity_LogBasis(double alpha, double beta, double x) {
  return - gsl_sf_lngamma(alpha) + alpha*(log(beta) + log(x)) - beta*x;
}


/**
Gamma-Gamma density with parameters log(a), log(a0), log(nu) 
generated by mixture x ~ Ga(a, lambda), lambda ~ Ga(a0, nu)
where ~Ga(shape, rate)
log density in log basis is also likelihood function
*/
double getGammaGammaLogDensity_LogBasis(const double a, const double a0, const double nu, 
					const int , const std::vector<double>& x, const double sumlogx) {
  double f = 0.0;
  const int N = x.size();
  for(int i = 0; i < N; ++i) {
    f += log( x[i] + nu );
  }
  f *= -(a + a0);
  f += N * (lngamma(a + a0) - lngamma(a) - lngamma(a0)) + a * sumlogx + N * a0 * log(nu);
  return f;
}

void gradientGammaGammaLogLikelihood(const double a, const double a0, const double nu, 
				     const int , const std::vector<double>& x, 
				     const double sumlogx, double* g) {
  double f1 = 0.0;
  double f2 = 0.0;
  const int N = x.size();
  for(int i = 0; i < N; ++i) {
    f1 += log( x[i] + nu );
    f2 += 1.0 / ( x[i] + nu );
  }
  g[0] = ( N * (digamma(a + a0) - digamma(a))  + sumlogx     - f1) ;
  g[1] = ( N * (digamma(a + a0) - digamma(a0)) + N * log(nu) - f1) ;
  g[2] =                                         N * a0/nu      - (a  + a0) * f2;
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

double getDirichletLogDensity_Softmax(const std::vector<double>& a, const bclib::pvector<double>& x) {
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

double MultinomialPDF(const std::vector<int>& r,
                      const std::vector<double>& theta) {

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



//-------------------------------------------------------------------------
///inverse softmax transformation.
/// Transforms proportions mu to numbers a on real line.
/// elements of a sum to 0
//-------------------------------------------------------------------------

void inv_softmax(size_t K, const double* const mu, double *a){
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



//-------------------------------------------------------------------------
/// softmax transformation
/// Inverse of inv_softmax transformation. 
/// elements of array a need not sum to zero 
//-------------------------------------------------------------------------

void softmax(size_t K, double *mu, const double* a){
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



//-------------------------------------------------------------------------
///partial inverse-softmax transformation.
///transformation is applied only to elements with b=true
//-------------------------------------------------------------------------

void inv_softmax(size_t K, const double* const mu, double *a, const bool* const b){
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



//-------------------------------------------------------------------------
///partial softmax transformation
///transformation is applied only to elements with b=true
//-------------------------------------------------------------------------

void softmax(size_t K, double *mu, const double* a, const bool* const b){
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



//-------------------------------------------------------------------------
// ************* Matrix Algebra **************************************
//-------------------------------------------------------------------------

///Caller for gsl_linalg_HH_solve.
///This function solves the system A x = b directly using Householder transformations. 
///On output the solution is stored in x and b is not modified. The matrix AA is destroyed by the Householder transformations. 
int HH_solve (size_t n, double *A, double *b, double *x)
{
  //create copy of A as gsl matrix; might not be necessary, depending on use
  gsl_matrix *AA = gsl_matrix_alloc(n,n);
  for(size_t i = 0;i < n*n; ++i)
    AA->data[i] = A[i];

  gsl_vector_view bb = gsl_vector_view_array(b, n);
  gsl_vector_view xx = gsl_vector_view_array(x, n);

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

///Caller for gsl_linalg_HH_svx.
/// This function solves the system A x = b in-place using Householder transformations. 
/// On input x should contain the right-hand side b, which is replaced by the solution on output. 
/// The matrix AA is destroyed by the Householder transformations. 
int HH_svx (size_t n, double *A, double *x)
{
  gsl_matrix *AA; //cannot use a matrix view because of above
  gsl_vector_view xx = gsl_vector_view_array(x,n); 

  //create copy of A as gsl matrix
  AA = gsl_matrix_alloc(n,n);
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


/// Compute LU factors for the matrix @a A, which is destroyed in the process.
int LU_decomp(size_t n, double *A, size_t *perm) {

  gsl_matrix_view AA = gsl_matrix_view_array(A, n, n);
  gsl_permutation pp = {n, perm};
  int signum;

  // disable default gsl error handler
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  int status = gsl_linalg_LU_decomp(&AA.matrix, &pp, &signum);

  // restore gsl error handler
  gsl_set_error_handler(old_handler);

  // check for failure
  if (status) {
    std::string errstring = "LU_decomp failed: ";
    errstring.append(gsl_strerror(status));
    throw(errstring);
  }

  return 0;
}

/// Solve a linear system of equations, where @a A and @a perm contain the LU
/// factors of the constraint matrix and the permutation matrix, respectively,
/// produced by LU_decomp().
int LU_solve(size_t n, double *A, size_t *perm, double *b, double *x) {

  gsl_matrix_view AA = gsl_matrix_view_array(A, n, n);
  gsl_vector_view bb = gsl_vector_view_array(b, n);
  gsl_vector_view xx = gsl_vector_view_array(x, n);
  gsl_permutation pp = {n, perm};

  // disable default gsl error handler
  gsl_error_handler_t* old_handler = gsl_set_error_handler_off();
  int status = gsl_linalg_LU_solve(&AA.matrix, &pp, &bb.vector, &xx.vector);

  // restore gsl error handler
  gsl_set_error_handler(old_handler);

  // check for failure
  if (status) {
    std::string errstring = "LU_solve failed: ";
    errstring.append(gsl_strerror(status));
    throw(errstring);
  }

  return 0;
}


///Computes the conditional mean and variance of a centred subvector of length kk of a zero-mean Multivariate Gaussian vector
///of length dim
void CentredGaussianConditional( size_t kk, double *mean, double *var,
				 double *newmean, double *newvar, size_t dim )
{
  if(dim == (kk+1)){//if conditioning on a scalar, call simpler version of function
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
  // gsl_matrix *mean1 = gsl_matrix_alloc(kk, 1);
  for(size_t i = 0; i< kk; ++i){
    // mean1->data[i] = mean[i];
    newmean[i] = mean[i];     //copy mean1 into new mean
  }
  gsl_matrix *mean2 = gsl_matrix_alloc(dim-kk, 1);
  for(size_t i = 0; i< dim-kk; ++i) mean2->data[i] = mean[i+kk];

  // gsl_matrix *Vaa = gsl_matrix_alloc(kk, kk);
  for(size_t i = 0; i < kk; ++i)
    for(size_t j = 0; j < kk; ++j){
      // Vaa->data[i*kk +j] = var[i*dim +j];
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
  // gsl_matrix_free(mean1);
  gsl_matrix_free(mean2);
  // gsl_matrix_free(Vaa);
  gsl_matrix_free(Vbb);
  gsl_matrix_free(Vab);
  if(status){
    std::string errstring = "CentredGaussianConditional failed, "; errstring.append(gsl_strerror(status));
    throw(errstring) ;//throw error message up and let caller decide what to do
  }

}

///special case of above with kk=dim-1
void CentredGaussianConditional( double *mean, double *var,
				 double *newmean, double *newvar, size_t dim ){
  for(unsigned i = 0; i < dim-1; ++i){
    newmean[i] = mean[i] - var[i*dim + dim-1] * mean[dim-1]/var[dim*dim-1];
    for(unsigned j = 0; j < dim-1; ++j)
      newvar[i*(dim-1)+j] = var[i*dim+j] - var[i*dim + dim-1]*var[(dim-1)*dim + j] / var[dim*dim-1];
  }
}

///returns the quadratic form in the density of the marginal distribution of the subvector of length kk of a zero-mean Gaussian vector of length dim. 
///Useful for score tests. 
double GaussianMarginalQuadraticForm( int kk, double *mean, double *var, size_t dim )
{
  int status = 0;
  string error_string;
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

  if(status){//stop if svx failed
    error_string = "Error in HH_svx in GaussianConditionalQuadraticForm:\n";
    error_string.append(gsl_strerror(status));
  }
  else{
    //set column of V to x
    for(int j = 0; j < kk; ++j)gsl_matrix_set(V, j, 0, gsl_vector_get(x, j));  
    
    //compute Q = U1' * V
    status = gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, U1, V, 0, Q);

    if(status){
      string error_string = "Error with matrix multiplication in GaussianConditionalQuadraticForm:\n";
      error_string.append(gsl_strerror(status));
    }
  }

  //clean up
  gsl_set_error_handler (old_handler);//restore gsl error handler 
  gsl_vector_free(x);
  
  gsl_matrix_free(V);
  gsl_matrix_free(U1);
  gsl_matrix_free(V11);

  double result = gsl_matrix_get(Q, 0,0);
  gsl_matrix_free(Q);

  if(status) throw(error_string);
  return result;
}

double GaussianQuadraticForm(double* mean, double* var, unsigned dim){
  double* VinvU = new double[dim];
  double result = 0.0;
  int status = -1;
  string err = "Error in GaussianQuadraticForm: ";
  try{
    status  = HH_solve(dim, var, mean, VinvU);
    for(unsigned i = 0; i < dim; ++i){
      result += mean[i] * VinvU[i];
    }
  }
  catch(string s){
    err += s;
    status = 1;
  }
  delete[] VinvU;
  if(status)throw err;
  return result;
}
///allocates space for 2-way rectangular array of doubles and initialise to zero
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

///allocates space for 2-way rectangular array of ints and initialises to zero
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

///delete double matrix, even nonrectangular
void free_matrix(double **M, int m){
  try{
    if(M){
      for(unsigned i = (unsigned)m; i > 0 ; --i) 
	if(M[i-1])
	  delete[] M[i-1];
      delete[] M;
    }
  }
  catch(...){
    string s = "Unable to delete matrix";
    throw(s);
  }
}

///deletes int matrix
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
///sets Sub as submatrix of M consisting of rows r1 to r2 and cols c1 to c2
void submatrix(double **M, double **Sub, int r1, int r2, int c1, int c2){
  for(int row = r1; row < r2; ++row)
    for(int col = c1; col < c2; ++col)
      Sub[row-r1][col-c1] = M[row][col];
}
void submatrix(double *M, double *Sub, int Mcols, int r1, int r2, int c1, int c2){
  for(int row = r1; row < r2; ++row)
    for(int col = c1; col < c2; ++col)
      Sub[(row-r1)*(r2-r1+1) + (col-c1)] = M[row*Mcols +col];
}
///equates two matrices
void equate_matrix(double **A, double **B, int m, int n){
  for(int i = 0; i < m; ++i)
    for(int j = 0; j < n; ++j)
      A[i][j] = B[i][j];
}
///adds two matrices.
///adds matrix b to matrix a
///dimensions must both be (d1 x d2)
void add_matrix(double *a, double *b, size_t d1, size_t d2){
  gsl_matrix_view A, B;
  A = gsl_matrix_view_array(a, d1, d2);
  B = gsl_matrix_view_array(b, d1, d2);
  gsl_matrix_add(&A.matrix, &B.matrix);
}
///multiplies two matrices.
///computes c = a * b, where all args are arrays representing matrices
/// and a is (d1 x d2) and b is (d2 x d3)
void matrix_product(double *a, double *b, double *c, size_t d1, size_t d2, size_t d3){
  gsl_matrix_view A, B, C;
  A = gsl_matrix_view_array(a, d1, d2);
  B = gsl_matrix_view_array(b, d2, d3);
  C = gsl_matrix_view_array(c, d1, d3);

  gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, &A.matrix, &B.matrix, 0, &C.matrix); 
}
///computes c = a * a'
void matrix_product(double *a, double *c, size_t d1, size_t d2){
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

///multiplies (d1 x d2) matrix a by c
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
///inverts a matrix using LU decomposition
void matrix_inverse(const double* const a, double* inv, size_t d){
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
///inverts a matrix using LU decomposition, overwriting original matrix
///can use this version to overwrite a with its inverse
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

///Cholesky decomposition, Crout algorithm
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
///prints a vector of doubles to screen
void print_vector(std::vector<double> a){
  copy(a.begin(), a.end(), ostream_iterator<double>(cout, " "));
  cout<<endl;
}
///prints a vector of strings to screen
void print_vector(std::vector<std::string> a){
  copy(a.begin(), a.end(), ostream_iterator<std::string>(cout, " "));
  cout<<endl;
}
///prints a vector of ints to screen
void print_vector(std::vector<int> a){
  copy(a.begin(), a.end(), ostream_iterator<int>(cout, " "));
  cout<<endl;
}
///log function with error handling
double eh_log(double x){
  //disable default gsl error handler
  gsl_error_handler_t* old_handler = gsl_set_error_handler(&bclib::GSLErrorHandler);
  const double result = gsl_sf_log(x);
  gsl_set_error_handler (old_handler);//restore gsl error handler 
  return result;
}
///exp function with error handling
double eh_exp(double x){
  gsl_error_handler_t* old_handler =  gsl_set_error_handler(&bclib::GSLErrorHandler);
  double result = 0.0;
  try{
    result = gsl_sf_exp(x);
  }
  catch(bclib::underflow){
    result = 0.0;
  }
  gsl_set_error_handler (old_handler);//restore gsl error handler 
  return result;

}
///lngamma function with error handling
double lngamma(double x){
  if(x <= 0.0)
    throw InfinityException("lngamma", __FILE__);
  //disable default gsl error handler
  gsl_error_handler_t* old_handler = gsl_set_error_handler(&bclib::GSLErrorHandler);
  const double result = gsl_sf_lngamma(x);
  //restore gsl error handler 
  gsl_set_error_handler (old_handler);
  return result;
}
///digamma function with error handling
double digamma(double x){
  if(x <= 0.0)
    throw InfinityException("digamma", __FILE__);
  //disable default gsl error handler
  gsl_error_handler_t* old_handler = gsl_set_error_handler(&bclib::GSLErrorHandler);
  const double result = gsl_sf_psi(x);
  gsl_set_error_handler (old_handler);//restore gsl error handler 
  return result;
}
///trigamma function with error handling
double trigamma(double x){
  if(x <= 0.0)
    throw InfinityException("trigamma", __FILE__);
  //disable default gsl error handler
  gsl_error_handler_t* old_handler = gsl_set_error_handler(&bclib::GSLErrorHandler);
  double result = gsl_sf_psi_n( 1, x);
  gsl_set_error_handler (old_handler);//restore gsl error handler 
  return result;
}
END_BCLIB_NAMESPACE
