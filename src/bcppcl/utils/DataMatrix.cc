/** 
 *   ADMIXMAP
 *   DataMatrix.cc 
 *   class to represent a matrix of data, possibly read in from file
 *   Copyright (c) 2005, 2006 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "DataMatrix.h"
#include <algorithm>
#include <iostream>
#include <numeric>
#include <sstream>

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
  if (row >= nrows || col >= ncols) throwBoundsViolation(row, col);
  return missing[row*ncols +col];
}
void DataMatrix::isMissing(unsigned row, unsigned col, bool b){
  if (row >= nrows || col >= ncols) throwBoundsViolation(row, col);
  missing[row*ncols + col] = b;
  if(b)anyMissing = true;
}
void DataMatrix::set(unsigned row, unsigned col, double x){
  if (row >= nrows || col >= ncols) throwBoundsViolation(row,col);
  data[row*ncols +col] = x;
}
double DataMatrix::DataMatrix::get(unsigned row, unsigned col)const{
  if (row >= nrows || col >= ncols) throwBoundsViolation(row,col);
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

///returns sample variance of jth col
double DataMatrix::getSampleVariance(int j, bool na_rm)const{
  double sum = 0.0, sumsq = 0.0, x = 0.0;
  int n = 0;
  for(unsigned i = 0; i < nrows; ++i){
    if(!na_rm || !missing[i*ncols+j]){
      ++n;
      x = data[i*ncols +j];
      sum += x;
      sumsq += x*x;
      }
  }
  return (sumsq - sum*sum / (double)n) / (double)n;
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
      if(isMissing(row, col)){
	set(row, col, mean[col]);
	//isMissing(row, col, false);
      }
}

void DataMatrix::Print()const{
  for(unsigned i = 0; i < nrows; ++i){
    //std::vector<double> a = getRow(i);
    //copy(a.begin(), a.end(), std::ostream_iterator<double>(std::cout, " "));
    for(unsigned j = 0; j < ncols; ++j)
      if(isMissing(i,j))std::cout << "# ";
      else std::cout << get(i,j) << " ";
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
}

// std::vector<double>::const_iterator DataMatrix::getData()const{
//   return data.begin();
// }
const double* DataMatrix::getData()const{
  return &(data[0]);
}
void DataMatrix::clear(){
  data.clear();
  missing.clear();
}
void DataMatrix::throwBoundsViolation(unsigned row, unsigned col)const{
  std::stringstream ss;
  ss << "Bounds Violation in DataMatrix: requested element (" << row+1 << ", " << col+1 << ") in a matrix with " << nrows << " rows and " << ncols << " cols"; 
  throw ss.str();
}
