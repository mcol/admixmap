// matrix%TYPE2%.h
#ifndef MATRIX_MATRIX%TYPE2%_H
#define MATRIX_MATRIX%TYPE2%_H 1

#include <iostream>
#include <iomanip>
#include <cassert>

#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_eigen.h>

#include <vector>

class Vector;
class Vector_d;
class Vector_i;

class Matrix;
class Matrix_d;
class Matrix_i;

class Matrix%TYPE2%
{
 private: // members
  std::vector<size_t> _MissingCol;
  std::vector<size_t> _MissingRow;
  gsl_matrix%GSL% * _matrix;
  
  // GLOBAL FRIENDS
  friend Matrix%TYPE2%
    ConcatenateHorizontally( const Matrix%TYPE2%& m1, const Matrix%TYPE2%& m2 );

  friend std::ostream&
    operator<< ( std::ostream& os, const Matrix%TYPE2%& M );

  friend void
    operator>> ( const char * FileName, Matrix%TYPE2%& M );

  friend Matrix%TYPE2%
    operator*( %TYPE%, const Matrix%TYPE2%& m);

 public:
  // CONSTRUCTORS
  // default constructor
  Matrix%TYPE2%();

  // copy constructor
  Matrix%TYPE2%(const Matrix%TYPE2%& M);

  // number of rows and columns
  Matrix%TYPE2%(size_t NewRows, size_t NewCols);


  // DESTRUCTOR
  ~Matrix%TYPE2%();


  // ASSIGNMENT OPERATORS
  Matrix%TYPE2%&
    operator=(const Matrix& m);

  Matrix%TYPE2%&
    operator=(const Matrix_d& m);

  Matrix%TYPE2%&
    operator=(const Matrix_i& m);

  
  // CONVERTERS
  Matrix_d
    Double() const;

  Matrix
    Float() const;

  Matrix_i
    Integer() const;


  // OPERATORS
  %TYPE%&
     operator()(int row, int col);

  %TYPE%&
     operator()(int row, int col) const;

  Matrix%TYPE2%
    operator+(%TYPE% f) const;

  Matrix%TYPE2%
    &operator*=(const %TYPE% f);

  Matrix%TYPE2%
    operator+(const Matrix& m) const;

  Matrix%TYPE2%
    operator+(const Matrix_i& m) const;

  Matrix%TYPE2%
    operator+(const Matrix_d& m) const;

  Matrix%TYPE2%&
    operator+=(const Matrix%TYPE2%& m);

  Matrix%TYPE2%
    operator-(const Matrix%TYPE2%& m) const;

  Matrix%TYPE2%
    operator*(%TYPE% f) const;

  Matrix%TYPE2%
    operator*(const Matrix%TYPE2%& m) const;

  Matrix%TYPE2%
    operator/(%TYPE%) const;

  Matrix%TYPE2%&
    operator/=(%TYPE%);

  // ACCESSORS
  
  Vector%TYPE2%
    GetColumn(int col) const;

  void
    SetColumn(int col, const Vector%TYPE2%& v);

  Vector%TYPE2%
    GetDiagonal() const;

  void
    SetDiagonal(%TYPE%);
  
  void
    SetElements(%TYPE% NewValue);

  void
    SetMissingElement(int row, int col);

  bool
    IsMissingValue(int row, int col) const;

  int
    GetMissingValueCol(int index) const;

  int
    GetMissingValueRow(int index) const;

  int
    GetNumberOfCols() const;

  int
    GetNumberOfRows() const;

  void
    SetNumberOfElements(size_t NewRows, size_t NewCols);

  int
    GetNumberOfMissingValues() const;

  Vector%TYPE2%
    GetRow(int row) const;

  void
    SetRow(int row, const Vector%TYPE2%& v);


  
  // FILE ACCESS
  void
    Load( const char *FileName );



  // OBJECT METHODS
  int
    CholeskyDecomposition(Matrix%TYPE2%*) const;

  Vector%TYPE2%
    ColumnMean() const;

  double
    Determinant() const;

  void
    Eigenvalue2(Vector%TYPE2%*);

  Matrix%TYPE2%
    FindElement( %TYPE% );
  
  void
    InvertUsingLUDecomposition();

  double
    LogDeterminant( int *sign ) const;

  Vector%TYPE2%
    RowMean() const;

  Vector%TYPE2%
    RowSum() const;

  void
    SetMissingValuesToColumnMeans();

  Matrix%TYPE2%
    SubMatrix( int rowstart, int rowfinish,
	       int colstart, int colfinish) const;   
  
  void
    SubMatrix2(int rowstart, int rowfinish,
	       int colstart, int colfinish);   
  
  void
    Symmetrize();
  
  Matrix%TYPE2%
    Transpose() const;
  
  bool
    IsSymmetric() const;

 private:

  // PRIVATE HELPER METHODS
  Vector%TYPE2%
    ColumnSum() const;

  Vector%TYPE2%
    ColumnMeanExcludingMissingValues() const;

  Vector%TYPE2%
    ColumnSumExcludingMissingValues() const;

  double
    LU_lndet();

  int
    LU_decomp (gsl_permutation * p);
};

#endif /* !defined MATRIX_MATRIX%TYPE%_H */
