// *-*-C++-*-*
#ifndef FUNCTIONS_H
#define FUNCTIONS_H 1

#include <vector>
#include <gsl/gsl_sf_psi.h>
#include "vector_d.h"
#include "matrix_d.h"

double getGammaLogDensity(double alpha, double beta, double x);

double getDirichletLogDensity(const Vector_d& a, const Vector_d& x);
double getDirichletLogDensity(const std::vector<double>& a, const Vector_d& x);

double getDirichletLogDensity(double *a, double *x, size_t K);
double getDirichletLogDensity(const std::vector<double>& a, double *x);
double getDirichletLogDensity_Softmax(const std::vector<double>& a, double *x);

double AverageOfLogs(const std::vector<double>& vec, double max);

void inv_softmax(size_t K, const double* const mu, double *a);
void softmax(size_t K, double *mu, const double* a);
void inv_softmax(size_t K, const double* const mu, double *a, const bool* const b);
void softmax(size_t K, double *mu, const double* a, const bool* const b);

//matrix algebra
int HH_solve (Matrix_d A, Vector_d b, Vector_d *x);
int HH_solve (size_t n, double *A, double *b, double *x);

int HH_svx (Matrix_d A, Vector_d *x);
int HH_svx (double *A, double *x);

void add_matrix(double *a, double *b, size_t d1, size_t d2);
void matrix_product(double *a, double *b, double *c, size_t d1, size_t d2, size_t d3);
void matrix_product(const double *a, const double *b, double *c, size_t d1, size_t d2, size_t d3);
void scale_matrix(double *a, const double c, size_t d1, size_t d2);

void CentredGaussianConditional( int kk, double *mean, double *var,
				 double *newmean, double *newvar, size_t dim );
double GaussianConditionalQuadraticForm( int kk, double *mean, double *var, size_t dim );

double **alloc2D_d(int m, int n);
int **alloc2D_i(int m, int n);
void free_matrix(double **, int);
void free_matrix(int **, int);
double **MatrixAsArray(Matrix_d &M);
void submatrix(double **M, double **Sub, int r1, int r2, int c1, int c2);
void submatrix(double *M, double *Sub, int Mcols, int r1, int r2, int c1, int c2);

double xlog(double x);
double xexp(double x);
#endif /* !FUNCTIONS_H */
