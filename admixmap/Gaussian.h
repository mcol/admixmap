// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   Gaussian.h (formerly gaussian_d.h)
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
#ifndef GAUSSIAN_H
#define GAUSSIAN_H 1

#include <fstream>


class Gaussian
{
public:

  Gaussian();
 
  Gaussian( int NewDimension );

  Gaussian( int NewDimension, const double* const NewMean, const double* const NewCovariance );

  ~Gaussian();

  void SetDimension( int NewDimension );

  Gaussian &operator=( const Gaussian &g );

  void Randomize();

  float GetCovariance( int Row, int Col ) const;

  double *GetMean() const;

  double &GetMean( int Element ) const;

  double* GetCovariance() const;

  int GetDimension() const;

  friend std::ostream& operator<< ( std::ostream& os, const Gaussian& G );

  void SetMean( int index, float NewValue );

  void SetMean(const double* const NewMean);

  void SetCovariance( int row, int col, float NewValue );

  void SetCovariance( const double* const NewCovariance );

  void SetCovariance( double NewValue );

  void Draw(double *beta)const;

private:

  int Dimension;

  int CovarianceDeterminantIsDirty;

  float CovarianceDeterminant;

  float NormalisationConstant;

  double *Mean;

  double *Covariance;

};

#endif
