#ifndef MATRIX_MATRIXARRAY_I_H
#define MATRIX_MATRIXARRAY_I_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>

class Matrix_i;
class MatrixArray;
class MatrixArray_d;

// array of Matrix_i objects
class MatrixArray_i
{
 private: // members
  int SizeOfArray;
  Matrix_i* TheArray;

 public:

  // CONSTRUCTORS

  // default constructor
  MatrixArray_i();
  // number of elements
  MatrixArray_i(int NumberOfElements);
  // copy constructor
  MatrixArray_i(const MatrixArray_i& M);
  // number of elements, rows and columns
  MatrixArray_i(int NumberOfElements,int NumberOfRows,int NumberOfCols);

  // DESTRUCTOR
  ~MatrixArray_i();


  MatrixArray Float();
  ///
  MatrixArray_d Double();
  ///
  void SetNumberOfElements( int NumberOfElements );
  ///
  void SetNumberOfElementsWithDimensions( int NumberOfElements, int NumberOfRows, int NumberOfCols );
  ///
  int GetNumberOfElements() const;
  ///
  void SetElements( const Matrix_i &TheMatrix );
  ///
  void SetElements( int NewValue );
  ///
  void SetDimensionsOfElements( int NumberOfRows, int NumberOfCols );
  ///
  friend std::ostream& operator<< ( std::ostream& os, const MatrixArray_i& M);
  ///
  Matrix_i &operator() ( int ElementNumber ) const;
  ///
  MatrixArray_i &operator=( const MatrixArray_i &m );
  ///
  MatrixArray_i &operator+=( const MatrixArray_i &m );
  ///
  void AddElementToEnd( const Matrix_i &NewElement );
  ///
  void RemoveElement( int ElementNumber );
};

#endif /* !defined MATRIX_MATRIXARRAY_I_H */
