// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   functions.h 
 *   Miscellaneous functions for admixmap, not belonging to any class
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
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
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>

// ************** Densities *******************
double getGammaLogDensity(double alpha, double beta, double x);
double getGammaLogDensity_LogBasis(double alpha, double beta, double x);

double getDirichletLogDensity(const double* const a, const double* const x, size_t K);
double getDirichletLogDensity(const std::vector<double>& a, const double* const x);
double getDirichletLogDensity(const std::vector<double>& a, const std::vector<double>& x);
double getDirichletLogDensity_Softmax(const std::vector<double>& a, const double* const x);

double MultinomialPDF( const std::vector<int> r, const std::vector<double> theta );

void ddigam( double *, double * );

void trigam( double *, double * );

// transformations
double AverageOfLogs(const std::vector<double>& vec, double max);

void inv_softmax(size_t K, const double* const mu, double *a);
void softmax(size_t K, double *mu, const double* a);
void inv_softmax(size_t K, const double* const mu, double *a, const bool* const b);
void softmax(size_t K, double *mu, const double* a, const bool* const b);

//matrix algebra
int HH_solve (size_t n, double *A, double *b, double *x);
int HH_svx (double *A, double *x);

void add_matrix(double *a, double *b, size_t d1, size_t d2);
void matrix_product(double *a, double *b, double *c, size_t d1, size_t d2, size_t d3);
void matrix_product(double *a, double *c, size_t d1, size_t d2);
void matrix_product(const double* const a, const double* const b, double* c, size_t d1, size_t d2, size_t d3);
void scale_matrix(double *a, const double c, size_t d1, size_t d2);

double determinant(double *a, size_t d);
void CentredGaussianConditional( int kk, double *mean, double *var,
				 double *newmean, double *newvar, size_t dim );
double GaussianConditionalQuadraticForm( int kk, double *mean, double *var, size_t dim );

void matrix_inverse(const double* const a, double* inv, size_t d);
void matrix_inverse(double* a, size_t d);
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
#endif /* !FUNCTIONS_H */
