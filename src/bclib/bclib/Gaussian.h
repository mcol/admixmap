// *-*-C++-*-*
/** 
 *   Gaussian.h
 *   This class represents a multivariate Gaussian distribution, including a sampler
 *   Copyright (c) 2002 - 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef GAUSSIAN_H
#define GAUSSIAN_H 1

#include "bclib/bclib.h"
#include <ostream>

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */



///Class to represent a multivariate Gaussian distribution and sample from it
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


/** @} */

END_BCLIB_NAMESPACE

#endif
