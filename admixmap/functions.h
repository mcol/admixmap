// *-*-C++-*-*
#ifndef FUNCTIONS_H
#define FUNCTIONS_H 1

#include <vector>
#include "vector_d.h"
#include "matrix_d.h"



double getGammaLogDensity(double alpha, double beta, double x);

double getDirichletLogDensity(const Vector_d& a, const Vector_d& x);

double AverageOfLogs(const std::vector<double>& vec, double max);

int HH_solve (Matrix_d A, Vector_d b, Vector_d *x);

int HH_svx (Matrix_d A, Vector_d *x);

void CentredGaussianConditional( int kk, Matrix_d mean, Matrix_d var, Matrix_d *newmean, Matrix_d *newvar );
void CentredGaussianConditional1( Vector_d mean, Matrix_d var, double *newmean, double  *newvar );

double **alloc2D_d(int m, int n);
int **alloc2D_i(int m, int n);
void free_matrix(double **, int);
void free_matrix(int **, int);
double **MatrixAsArray(Matrix_d &M);
#endif /* !FUNCTIONS_H */
