/** 
 *   ADMIXMAP
 *   Gaussian.cc (formerly gaussian_d.cc)
 *   This class represents a multivariate Gaussian distribution, including a sampler
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <stdlib.h>
#include <algorithm>
#include "rand.h"
#include "functions.h"
#include "Gaussian.h"

using namespace::std;

// *************** Constructors *************************************************
Gaussian::Gaussian()
{
  Mean = 0;
  Covariance = 0;
  SetDimension(1);
}

Gaussian::Gaussian( int NewDimension )
{
  Mean = 0;
  Covariance = 0;
  if ( NewDimension < 1 )SetDimension(1);
  else SetDimension(NewDimension);
}

Gaussian::Gaussian( int NewDimension, const double* const NewMean, const double* const NewCovariance )
{
  Mean = 0;
  Covariance = 0;
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

  if ( Dimension == NewDimension )
    return;

  delete[] Mean;
  delete[] Covariance;

  if ( NewDimension < 1 )
    NewDimension = 1;
  Dimension = NewDimension;

  // Allocate vector containing mean
  Mean = new double[ Dimension ];

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
// *************** Stream Insertion Operator (for printing) **********************
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

// *************** Draw from Gaussian Distribution **********************
void Gaussian::Draw(double *beta)const
{
  int i;
  double draw[Dimension];
  double L[Dimension*Dimension];
 
  for ( i = 0; i < Dimension; i++ )
    draw[i] = gennor( (double)0.0, (double)1.0 );
  
  //cholesky decomposition
  cholDecomp(Covariance, L, Dimension);

  for ( i = 0; i < Dimension; i++ ){
    beta[i] = 0.0;
    for(int j = 0; j < Dimension; ++j)beta[i] += L[i*Dimension + j]*draw[j];
    beta[i] += Mean[i];
  }
}
