// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   DataMatrix.h 
 *   class to represent a matrix of data, possibly read in from file
 *   Copyright (c) 2005 LSHTM
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
#ifndef DATAMATRIX_H
#define DATAMATRIX_H 1
#include <vector>

class DataMatrix{
public:
  DataMatrix();
  DataMatrix(unsigned, unsigned);
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
  void Print()const;
  //std::vector<double>::const_iterator getData() const;
  const double* const getData() const;// allows access to data without danger of changing it
private:
  std::vector<double > data;
  std::vector<bool > missing;
  unsigned nrows;
  unsigned ncols;
  bool anyMissing;
  class BoundsViolation { };

};
#endif
