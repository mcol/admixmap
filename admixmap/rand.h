// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   rand.h 
 *   Serial random number generators for admixmap
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include <iostream>
#include <vector>

///
double myrand();
///
double myrandRange( double Min, double Max );
///
void smyrand( long seed );
///
double gengam(double aa,double bb);
///
double genbet(double aa,double bb);
///
double gennor(double av,double sd);
//
int genbinomial( int n, double p );
//
unsigned int genpoi( double );
//
std::vector<int> genmultinomial2( int n, const std::vector<double> p );

///Poisson generator
long ignpoi(double mu);
///
int SampleFromDiscrete( const double probs[] , int numberofelements);
///
void gendirichlet(const size_t K, const double alpha[], double theta[] );
