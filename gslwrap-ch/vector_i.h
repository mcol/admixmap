#ifndef MATRIX_VECTOR_I_H
#define MATRIX_VECTOR_I_H 1

#if !defined( LINELENGTH )
#define LINELENGTH 300000
#endif

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <limits.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas_types.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sf_gamma.h>

#if !defined( MAX_LIN_ACTIVITY )
#define MAX_LIN_ACTIVITY 700
#endif

#if !defined( MIN_LIN_ACTIVITY )
#define MIN_LIN_ACTIVITY -MAX_LIN_ACTIVITY
#endif

class Matrix_i;
class Vector;
class Vector_d;

// vector of integers
class Vector_i
{
 private: // members
  int NumberOfElements;
  int NumberOfMissingValues;
  int *MissingValues;
  gsl_vector_int *vector;

 public:

  // CONSTRUCTORS

  // default constructor
  Vector_i();
  // copy constructor
  Vector_i(const Vector_i& V);
  // initial length
  Vector_i(int NewNumberOfElements);
  // load from file
  Vector_i(const char* FileName);

  // DESTRUCTOR
  ~Vector_i();

  Vector_d Double();
  ///
  Vector Float();
  ///
  int GetElement( int index ) const;
  ///
  int& operator()(int);
  int& operator()(int) const;

  // same as op() and op()const
  int& operator[](int) const;
  int& operator[](int);

  void SetElement( int index, int NewValue );
  ///
  int GetNumberOfElements() const;
  ///
  void SetNumberOfElements( int NewValue );
  ///
  friend std::ostream& operator<< ( std::ostream& os, const Vector_i& V );
  ///
  friend std::istream& operator>> ( std::istream& is, Vector_i& V );
  ///
  Vector_i operator+( const Vector_i &v ) const;
  ///
  Vector_i operator+( int f ) const;
  ///
  friend Vector_i operator+( int f, const Vector_i &v );
  ///
  Vector_i &operator+=( int f );
  ///
  Vector_i &operator+=( const Vector_i &v );
  ///
  Vector_i operator-( const Vector_i &v ) const;
  ///
  Vector_i operator-( int f ) const;
  ///
  friend Vector_i operator-( int f, const Vector_i &v );
  ///
  Vector_i &operator-=( int f );
  ///
  Vector_i &operator-=( const Vector_i &v );
  ///
  Vector_i operator-() const;
  ///
  Vector_i operator*( const Vector_i &v ) const;
  ///
  Vector_i operator*( int f ) const;
  ///
  friend Vector_i operator*( int f, const Vector_i &v );
  ///
  Vector_i &operator*=( int f );
  ///
  Vector_i &operator*=( const Vector_i &v );
  ///
  Vector_i operator/( const Vector_i &v ) const;
  ///
  Vector_i operator/( int f ) const;
  ///
  friend Vector_i operator/( int f, const Vector_i &v );
  ///
  Vector_i &operator/=( int f );
  ///
  Vector_i &operator/=( const Vector_i &v );
  ///
  Vector_i &operator=( const Vector_i &v );
  ///
  int operator==( const Vector_i &v ) const;
  ///
  void SetElements( int );
  ///
  int Sum() const;
  ///
  double Mean() const;
  ///
  double Variance(void) const;
  ///
  Vector_i Absolute(void) const;
  ///
  double Norm();
  ///
  int MinimumElement(void);
  ///
  int MaximumElement(void);
  ///
  int MinimumAbsoluteElement(void);
  ///
  int MaximumAbsoluteElement(void);
  ///
  void Load( const char *FileName );
  ///
  Matrix_i MyTransposeTimesMe() const;
  ///
  Matrix_i RowMatrix() const;
  ///
  Matrix_i ColumnMatrix() const;
  ///
  int GetNumberOfMissingValues() const;
  ///
  int GetMissingValueIndex( int MissingValueNumber ) const;
  ///
  int IsMissingElement( int index ) const;
  ///
  int CumulativeSumDraw() const;
  ///
  void Distinct();
  ///
  void Randomize( int, int );
  ///
  void RemoveElement( int );
  ///
  void AddElement( int );
  ///
  gsl_vector_int *GetGslVector();
};

#endif /* !MATRIX_VECTOR_I_H */

