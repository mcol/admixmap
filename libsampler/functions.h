// *-*-C++-*-*
/** 
 *   functions.h 
 *   Miscellaneous functions, not belonging to any class
 *   Copyright (c) 2002 - 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef FUNCTIONS_H
#define FUNCTIONS_H 1

#include <math.h>
#include <vector>
#include <string>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

// ************** Densities *******************

/// log of gamma(alpha, beta) density at x
double getGammaLogDensity(double alpha, double beta, double x);
/// log of gamma density, where X is the log of a gamma variable
double getGammaLogDensity_LogBasis(double alpha, double beta, double x);
/// log of Dirichlet density of dimension K with parameters a, x and a as arrays
double getDirichletLogDensity(const double* const a, const double* const x, size_t K);
/// log of Dirichlet density of dimension K with parameters a, a as vector, x as array
double getDirichletLogDensity(const std::vector<double>& a, const double* const x);
/// log of Dirichlet density of dimension K with parameters a, x and a as vectors
double getDirichletLogDensity(const std::vector<double>& a, const std::vector<double>& x);
/// log of Dirichlet density of dimension K with parameters a, where a softmax transformation has been applied to x
double getDirichletLogDensity_Softmax(const std::vector<double>& a, const double* const x);
/// multinomial density qith probabilities r
double MultinomialPDF( const std::vector<int> r, const std::vector<double> theta );

// transformations
double AverageOfLogs(const std::vector<double>& vec, double max);

///inverse-softmax transformation
void inv_softmax(size_t K, const double* const mu, double *a);
///softmax transformation
void softmax(size_t K, double *mu, const double* a);
///partial inverse-softmax transformation
void inv_softmax(size_t K, const double* const mu, double *a, const bool* const b);
///partial softmax transformation
void softmax(size_t K, double *mu, const double* a, const bool* const b);

//matrix algebra
///wrapper for gsl_HH_solve
int HH_solve (size_t n, double *A, double *b, double *x);
///wrapper for gsl_HH_svx
int HH_svx (double *A, double *x);

///adds two matrices
void add_matrix(double *a, double *b, size_t d1, size_t d2);
///multiplies two matrices
void matrix_product(double *a, double *b, double *c, size_t d1, size_t d2, size_t d3);
///multiplies a matrix by its transpose
void matrix_product(double *a, double *c, size_t d1, size_t d2);
///multiplies two matrices, const version
void matrix_product(const double* const a, const double* const b, double* c, size_t d1, size_t d2, size_t d3);
///multiplies a matrix by a scalar
void scale_matrix(double *a, const double c, size_t d1, size_t d2);
///computes the determinant of a square matrix
double determinant(double *a, size_t d);

void CentredGaussianConditional( size_t kk, double *mean, double *var,
				 double *newmean, double *newvar, size_t dim );
void CentredGaussianConditional( double *mean, double *var,
				 double *newmean, double *newvar, size_t dim );
double GaussianMarginalQuadraticForm( int kk, double *mean, double *var, size_t dim );
double GaussianQuadraticForm(double* mean, double* var, unsigned dim);

///inverts a matrix using LU decomposition
void matrix_inverse(const double* const a, double* inv, size_t d);
///inverts a matrix using LU decomposition, overwriting original matrix
void matrix_inverse(double* a, size_t d);
///Cholesky decomposition, Crout algorithm
int cholDecomp(const double* const a, double *L, int n);

double **alloc2D_d(int m, int n);
int **alloc2D_i(int m, int n);
void free_matrix(double **, int);
void free_matrix(int **, int);

void submatrix(double **M, double **Sub, int r1, int r2, int c1, int c2);
void submatrix(double *M, double *Sub, int Mcols, int r1, int r2, int c1, int c2);

void print_vector(std::vector<double> a);
void print_vector(std::vector<int> a);
void print_vector(std::vector<std::string> a);

double xlog(double x);
double xexp(double x);
///log function with error handling
double mylog(double x);
///exp function with error handling
double myexp(double x);
///lngamma function with error handling
double lngamma(double x);
///digamma function with error handling
double digamma(double x);
///trigamma function with error handling
double trigamma(double x);
#endif /* !FUNCTIONS_H */
