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
#include <vector>

class DataMatrix{
public:
  DataMatrix();
  DataMatrix(unsigned, unsigned);
  void setDimensions(unsigned, unsigned);
  bool isMissing(unsigned, unsigned);
  void isMissing(unsigned, unsigned, bool);
  std::vector<double> getRow(unsigned);
  double get(unsigned, unsigned) const;
  void set(unsigned, unsigned, double);
  unsigned nRows()const;
  unsigned nCols()const;
  DataMatrix SubMatrix(unsigned, unsigned, unsigned, unsigned);

private:
  std::vector<double > data;
  std::vector<bool > missing;
  unsigned nrows;
  unsigned ncols;
  class BoundsViolation { };

};
