#ifndef MATRIX_MATRIXARRAY_D_H
#define MATRIX_MATRIXARRAY_D_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstring>

class Matrix_d;
class MatrixArray;
class MatrixArray_i;

// an array of Matrix_d objects
class MatrixArray_d
{
 private: // members
  int SizeOfArray;
  Matrix_d* TheArray;

 public:

  // CONSTRUCTORS

  // default constructor
  MatrixArray_d();
  // number of elements
  MatrixArray_d(int NumberOfElements);
  // copy constructor
  MatrixArray_d(const MatrixArray_d& M);
  // number of elements, rows and columns
  MatrixArray_d(int NumberOfElements,int NumberOfRows,int NumberOfCols);

  // DESTRUCTOR
  ~MatrixArray_d();


  MatrixArray Float();
  ///
  MatrixArray_i Integer();
  ///
  void SetNumberOfElements( int NumberOfElements );
  ///
  void SetNumberOfElementsWithDimensions( int NumberOfElements, int NumberOfRows, int NumberOfCols );
  ///
  int GetNumberOfElements() const;
  ///
  void SetElements( const Matrix_d &TheMatrix );
  ///
  void SetElements( float NewValue );
  ///
  void SetDimensionsOfElements( int NumberOfRows, int NumberOfCols );
  ///
  friend std::ostream& operator<< ( std::ostream& os, const MatrixArray_d& M);
  ///
  Matrix_d &operator() ( int ElementNumber ) const;
  ///
  MatrixArray_d &operator=( const MatrixArray_d &m );
  ///
  MatrixArray_d &operator+=( const MatrixArray_d &m );
  ///
  MatrixArray_d operator/( double f ) const;
  ///
  void AddElementToEnd( const Matrix_d &NewElement );
  ///
  void RemoveElement( int ElementNumber );

};

#endif /* !defined MATRIX_MATRIXARRAY_D_H */
