/** 
 *   ADMIXMAP
 *   DataMatrix.cc 
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
#include "DataMatrix.h"
#include <algorithm>
#include <iostream>
#include <numeric>

DataMatrix::DataMatrix(){
  nrows = 0;
  ncols = 0;
  anyMissing =false;
}
DataMatrix::DataMatrix(unsigned rows, unsigned cols){
  nrows = rows;
  ncols = cols;
  data.resize(nrows*ncols);
  fill(data.begin(), data.end(), 0.0);
  missing.resize(nrows*ncols);
  fill(missing.begin(), missing.end(), false);
  anyMissing = false;
}
void DataMatrix::setDimensions(unsigned rows, unsigned cols){
  nrows = rows;
  ncols = cols;
  data.resize(nrows*ncols);
  fill(data.begin(), data.end(), 0.0);
  missing.resize(nrows*ncols);
  fill(missing.begin(), missing.end(), false);
  anyMissing = false;
}
unsigned DataMatrix::nRows()const{
  return nrows;
}
unsigned DataMatrix::nCols()const{
  return ncols;
}
bool DataMatrix::isMissing(unsigned row, unsigned col)const{
  if (row >= nrows || col >= ncols) throw BoundsViolation();
  return missing[row*ncols +col];
}
void DataMatrix::isMissing(unsigned row, unsigned col, bool b){
  if (row >= nrows || col >= ncols) throw BoundsViolation();
  missing[row*ncols + col] = b;
  if(b)anyMissing = true;
}
void DataMatrix::set(unsigned row, unsigned col, double x){
  if (row >= nrows || col >= ncols) throw BoundsViolation();
  data[row*ncols +col] = x;
}
double DataMatrix::DataMatrix::get(unsigned row, unsigned col)const{
  if (row >= nrows || col >= ncols) throw BoundsViolation();
  return data[row*ncols +col];
}
std::vector<double> DataMatrix::getRow(unsigned r)const{
  std::vector<double> row(ncols);
  std::vector<double>::const_iterator it = data.begin()+r*ncols;
  copy(it, it + ncols, row.begin());
  return row;
}
std::vector<double> DataMatrix::getCol(unsigned c)const{
  std::vector<double> col;
  for(unsigned row = 0; row < nrows; ++row)col.push_back( data[row*ncols + c] );
  return col;
}
std::vector<double> DataMatrix::columnMeans()const{
  std::vector<double> mean(ncols);
  for( unsigned j = 0; j < ncols; j++ ){
    int count = 0;
    mean[j] = 0.0;
    for(unsigned i = 0; i < nrows; i++ )
      if(!isMissing(i,j)){
      mean[j] += get(i,j);
      ++count;
      }
    mean[j] /= (double)count;
    }
  return mean;
}
DataMatrix DataMatrix::SubMatrix(unsigned r1, unsigned r2, unsigned c1, unsigned c2){
  if( r1>r2 || c1>c2 || r2 > nrows-1 || c2 > ncols-1)
    std::cerr<<"Error in DataMatrix::SubMatrix"<<std::endl;
  DataMatrix Sub(r2-r1+1, c2-c1+1);
  for(unsigned i = r1; i<= r2; ++i)
    for(unsigned j = c1; j <= c2;++j){
      Sub.set(i-r1,j-c1, data[i*ncols +j]);
      Sub.isMissing(i-r1,j-c1, missing[i*ncols+j]);
    }
  return Sub;   
}

void DataMatrix::SetMissingValuesToColumnMeans(){
  std::vector<double> mean = columnMeans();

  for(unsigned row = 0; row < nrows; ++row)
    for(unsigned col = 0; col < ncols; ++col)
      if(isMissing(row, col))set(row, col, mean[col]);
}

void DataMatrix::Print()const{
  for(unsigned i = 0; i < nrows; ++i){
    std::vector<double> a = getRow(i);
    copy(a.begin(), a.end(), std::ostream_iterator<double>(std::cout, " "));
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
}

// std::vector<double>::const_iterator DataMatrix::getData()const{
//   return data.begin();
// }
const double* const DataMatrix::getData()const{
  return &(data[0]);
}
