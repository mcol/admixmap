// *-*-C++-*-*
/** 
 *   linalg.h 
 *   linear algebra functions
 *   Copyright (c) 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef LINALG_H
#define LINALG_H 1

#include "bclib/bclib.h"

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */



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

void submatrix(double **M, double **Sub, int r1, int r2, int c1, int c2);
void submatrix(double *M, double *Sub, int Mcols, int r1, int r2, int c1, int c2);

double **alloc2D_d(int m, int n);
int **alloc2D_i(int m, int n);
void free_matrix(double **, int);
void free_matrix(int **, int);


/** @} */

END_BCLIB_NAMESPACE

#endif /* !LINALG_H */
