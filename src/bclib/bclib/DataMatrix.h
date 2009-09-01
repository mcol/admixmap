// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   DataMatrix.h 
 *   class to represent a matrix of data, possibly read in from file
 *   Copyright (c) 2005, 2006 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef DATAMATRIX_H
#define DATAMATRIX_H 1
#include "bclib/bclib.h"
#include <vector>

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */



///Class to represent a matrix of data, possibly read in from file
class DataMatrix{
public:
  DataMatrix();
  DataMatrix(unsigned, unsigned);
  void clear();//clears vectors
  void setDimensions(unsigned, unsigned);
  bool isMissing(unsigned, unsigned)const;
  void isMissing(unsigned, unsigned, bool);
  bool hasMissing()const{return anyMissing;}
  std::vector<double> getRow(unsigned)const;
  std::vector<double> getCol(unsigned c)const;
  std::vector<double> columnMeans()const;
  double get(unsigned, unsigned) const;
  void set(unsigned, unsigned, double);
  unsigned nRows()const;
  unsigned nCols()const;
  DataMatrix SubMatrix(unsigned, unsigned, unsigned, unsigned);
  void SetMissingValuesToColumnMeans();
  double getSampleVariance(int j, bool na_rm=false)const;
  void Print()const;
  //std::vector<double>::const_iterator getData() const;
  /// allows access to data without danger of changing it
  const double* getData() const;

private:
  std::vector<double > data;
  std::vector<bool > missing;
  unsigned nrows;
  unsigned ncols;
  bool anyMissing;
  void throwBoundsViolation(unsigned, unsigned)const;
  //void operator=(const DataMatrix&);
};

/** @} */

END_BCLIB_NAMESPACE

#endif
