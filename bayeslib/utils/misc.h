// *-*-C++-*-*
/** 
 *   misc.h 
 *   Miscellaneous functions, not belonging to any class
 *   Copyright (c) 2006 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef MISCFUNCTIONS_H
#define MISCFUNCTIONS_H 1

#include <vector>
#include <string>

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

double **alloc2D_d(int m, int n);
int **alloc2D_i(int m, int n);
void free_matrix(double **, int);
void free_matrix(int **, int);
#endif /* !MISCFUNCTIONS_H */
