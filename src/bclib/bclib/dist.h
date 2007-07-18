// *-*-C++-*-*
/** 
 *   dist.h 
 *   density functions
 *   Copyright (c) 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef DIST_H
#define DIST_H 1

#include "bclib/bclib.h"
#include <vector>

BEGIN_BCLIB_NAMESPACE

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
double getGammaGammaLogDensity_LogBasis(const double a, const double a0, const double nu, 
					const int N, const std::vector<double>& x, const double sumlogx);
void gradientGammaGammaLogLikelihood(const double a, const double a0, const double nu, 
				     const int N, const std::vector<double>& x, 
				     const double sumlogx, double* g);

END_BCLIB_NAMESPACE
#endif /* !DIST_H */
