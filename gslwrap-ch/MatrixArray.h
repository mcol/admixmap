#ifndef MATRIX_MATRIXARRAY_H
#define MATRIX_MATRIXARRAY_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include "matrix.h"
#include "matrix_d.h"
#include "matrix_i.h"
#include "MatrixArray_i.h"
#include "MatrixArray_d.h"

class Matrix;
class MatrixArray_d;
class MatrixArray_i;

// array of generic matrices
class MatrixArray
{
 private: // members
  int SizeOfArray;
  Matrix* TheArray;

 public:

  // CONSTRUCTORS

  // default constructor
  MatrixArray();
  // number of elements
  MatrixArray(int NumberOfElements);
  // copy constructor
  MatrixArray(const MatrixArray& M);
  // number of elements, rows and columns
  MatrixArray(int NumberOfElements,int NumberOfRows,int NumberOfCols);

  // DESTRUCTOR
  ~MatrixArray();


  MatrixArray_d Double();
  ///
  MatrixArray_i Integer();
  ///
  void SetNumberOfElements( int NumberOfElements );
  ///
  void SetNumberOfElementsWithDimensions( int NumberOfElements, int NumberOfRows, int NumberOfCols );
  ///
  int GetNumberOfElements() const;
  ///
  void SetElements( const Matrix &TheMatrix );
  ///
  void SetElements( float NewValue );
  ///
  void SetDimensionsOfElements( int NumberOfRows, int NumberOfCols );
  ///
  friend std::ostream& operator<< ( std::ostream& os, const MatrixArray& M);
  ///
  Matrix &operator() ( int ElementNumber ) const;
  ///
  MatrixArray operator/( double f ) const;
  ///
  MatrixArray &operator=( const MatrixArray &m );
  ///
  MatrixArray &operator+=( const MatrixArray &m );
  ///
  void AddElementToEnd( const Matrix &NewElement );
  ///
  void RemoveElement( int ElementNumber );

};

#endif /* !defined MATRIX_MATRIXARRAY_H */
