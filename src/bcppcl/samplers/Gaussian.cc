/** 
 *   Gaussian.cc
 *   This class represents a multivariate Gaussian distribution, including a sampler
 *   Copyright (c) 2002 - 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <stdlib.h>
#include <algorithm>
#include "rand.h"
#include "utils/linalg.h"//for Cholesky decomp
#include "Gaussian.h"

using namespace::std;

// *************** Constructors *************************************************
Gaussian::Gaussian()
{
  Mean = 0;
  Covariance = 0;
  Dimension = 1;
  SetDimension(1);
}

Gaussian::Gaussian( int NewDimension )
{
  Mean = 0;
  Covariance = 0;
  if ( NewDimension < 1 ){
    Dimension = 1;
    SetDimension(1);
  }
  else SetDimension(NewDimension);
}

Gaussian::Gaussian( int NewDimension, const double* const NewMean, const double* const NewCovariance )
{
  Mean = 0;
  Covariance = 0;
  Dimension = NewDimension;
  SetDimension(NewDimension);

  SetMean(NewMean);
  SetCovariance(NewCovariance);
}

// *************** Destructor *************************************************
Gaussian::~Gaussian()
{
  delete[] Mean;
  delete Covariance;
}

// *************** Setup *************************************************
void Gaussian::SetDimension( int NewDimension )
{

//   if ( Dimension == NewDimension )
//     return;
  if ( NewDimension < 1 )
    NewDimension = 1;
  Dimension = NewDimension;

  // Allocate vector containing mean
  delete[] Mean;
  Mean = new double[ Dimension ];

  delete[] Covariance;
  // Allocate matrix containing covariance matrix
  Covariance = new double[ Dimension * Dimension ];
  CovarianceDeterminantIsDirty = 1;



}

void Gaussian::SetMean(const double* const NewMean){
  for(int i = 0; i < Dimension; ++i)Mean[i] = NewMean[i];
}
void Gaussian::SetCovariance( int row, int col, float NewValue )
{
  Covariance[row * Dimension +col] = NewValue;
  if ( row != col )
    Covariance[col * Dimension +row] = NewValue;
  CovarianceDeterminantIsDirty = 1;
}

void Gaussian::SetCovariance( const double* const NewCovariance )
{
  for(int i = 0; i < Dimension*Dimension; ++i)Covariance[i] = NewCovariance[i];
  CovarianceDeterminantIsDirty = 1;
}

void Gaussian::SetCovariance( double NewValue )
{
  fill(Covariance, Covariance+Dimension*Dimension, NewValue );
  CovarianceDeterminantIsDirty = 1;
}

// *************** Assignment Operator (probably not needed) **********************
Gaussian &Gaussian::operator=( const Gaussian &g )
{
  if ( &g == this )
    return *this;

  SetDimension(g.GetDimension());
  SetMean(g.GetMean());
  SetCovariance(g.GetCovariance());

  return *this;
}
/// Stream Insertion Operator (for printing)
std::ostream& operator<< ( std::ostream& os, const Gaussian& G )
{
  for(int i = 0; i < G.GetDimension(); ++i)os << G.GetMean()[i]<<" ";
  os<<endl<<endl;;
  os << G.GetCovariance();

  return os;
}	

// *************** ACCESSORS **********************
float Gaussian::GetCovariance( int Row, int Col ) const
{
  return( Covariance[ Row*Dimension + Col ] );
}

double* Gaussian::GetCovariance() const
{
  return( Covariance );
}

double *Gaussian::GetMean() const
{
  return( Mean );
}

double &Gaussian::GetMean( int Element ) const
{
  return( Mean[ Element ] );
}

int Gaussian::GetDimension() const
{
  return( Dimension );
}

/// Draw from Gaussian Distribution
void Gaussian::Draw(double *beta)const
{
  int i;
  double* draw = new double[Dimension];
  double* L = new double[Dimension*Dimension];
 
  for ( i = 0; i < Dimension; i++ )
    draw[i] = Rand::gennor( (double)0.0, (double)1.0 );
  
  //cholesky decomposition
  if(cholDecomp(Covariance, L, Dimension)) {
    cerr << "Cholesky decomposition failed in Gaussian::Draw"<<endl;
    exit(1);
  }

  for ( i = 0; i < Dimension; i++ ){
    beta[i] = 0.0;
    for(int j = 0; j < Dimension; ++j)beta[i] += L[i*Dimension + j]*draw[j];
    beta[i] += Mean[i];
  }
  delete[] L;
  delete[] draw;
}
