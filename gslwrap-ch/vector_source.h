// vector%TYPE2%.h
#ifndef MATRIX_VECTOR%TYPE%_H
#define MATRIX_VECTOR%TYPE%_H 1

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cfloat>

#if !defined( LINELENGTH )
#define LINELENGTH 300000
#endif

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

class Vector;
class Vector_d;
class Vector_i;

class Matrix;
class Matrix_d;
class Matrix_i;

class Vector%TYPE2%
{
 private: //members
  int NumberOfElements;
  int NumberOfMissingValues;
  int *MissingValues;
  gsl_vector%GSL% *vector;

  // GLOBAL FRIENDS
  friend class Matrix%TYPE2%;
  ///
  friend Vector%TYPE2% operator+( %TYPE% f, const Vector%TYPE2% &v );
  ///
  friend Vector%TYPE2% operator-( %TYPE% f, const Vector%TYPE2% &v );
  ///
  friend Vector%TYPE2% operator*( %TYPE% f, const Vector%TYPE2% &v );
  ///
  friend Vector%TYPE2% operator/( %TYPE% f, const Vector%TYPE2% &v );
  ///
  friend std::ostream& operator<< ( std::ostream& os, const Vector%TYPE2%& V );
  ///
  friend std::istream& operator>> ( std::istream& is, Vector%TYPE2%& V );
  ///
 public:

  // CONSTRUCTORS

  // default constructor
  Vector%TYPE2%();
  // initial length
  Vector%TYPE2%(int NewNumberOfElements);
  // load from file
  Vector%TYPE2%(const char* FileName);
  // copy constructor
  Vector%TYPE2%(const Vector%TYPE2%& V);

  // DESTRUCTOR
  ~Vector%TYPE2%();

  Vector Float();

  Vector_d Double();
  ///
  %TYPE% GetElement( int index ) const;
  ///
  %TYPE% &operator()(int);
  %TYPE% &operator()(int) const;

  // synonms for op() and op()const
  %TYPE% &operator[](int);
  %TYPE% &operator[](int) const;
  ///
  void SetElement( int index, %TYPE% NewValue );
  ///
  int GetNumberOfElements() const;
  ///
  void SetNumberOfElements( int NewValue );
  Vector%TYPE2% operator+( const Vector%TYPE2% &v ) const;
  ///
  Vector%TYPE2% operator+( %TYPE% f ) const;
  ///
  Vector%TYPE2% &operator+=( %TYPE% f );
  ///
  Vector%TYPE2% &operator+=( const Vector%TYPE2% &v );
  ///
  Vector%TYPE2% operator-( const Vector%TYPE2% &v ) const;
  ///
  Vector%TYPE2% operator-( %TYPE% f ) const;
  ///
  Vector%TYPE2% &operator-=( %TYPE% f );
  ///
  Vector%TYPE2% &operator-=( const Vector%TYPE2% &v );
  ///
  Vector%TYPE2% operator-() const;
  ///
  Vector%TYPE2% operator*( const Vector%TYPE2% &v ) const;
  ///
  Vector%TYPE2% operator*( %TYPE% f ) const;
  ///
  Vector%TYPE2% &operator*=( %TYPE% f );
  ///
  Vector%TYPE2% &operator*=( const Vector%TYPE2% &v );
  ///
  Vector%TYPE2% operator/( const Vector%TYPE2% &v ) const;
  ///
  Vector%TYPE2% operator/(double f) const;
  ///
  Vector%TYPE2% &operator/=( %TYPE% f );
  ///
  Vector%TYPE2% &operator/=( const Vector%TYPE2% &v );
  ///
  Vector%TYPE2% &operator=( const Vector%TYPE2% &v );
  ///
  int operator==( const Vector%TYPE2% &v ) const;
  ///
  void Randomize( %TYPE% LowerLimit, %TYPE% UpperLimit );
  ///
  void RandomizeGaussian(%TYPE%, %TYPE%);
  ///
  void Sigmoid(void);
  ///
  void Exp(void);
  ///
  void Tanh(void);
  ///
  void SetElements( %TYPE% Value );
  ///
  %TYPE% Sum() const;
  ///
  %TYPE% Mean() const;
  ///
  %TYPE% Variance(void) const;
  ///
  Vector%TYPE2% Absolute(void) const;
  ///
  %TYPE% Norm();
  ///
  %TYPE% MinimumElement(void);
  ///
  %TYPE% MaximumElement(void);
  ///
  int WhereMinimumElement();
  ///
  int WhereMaximumElement();
  ///
  %TYPE% MinimumAbsoluteElement(void);
  ///
  %TYPE% MaximumAbsoluteElement(void);
  ///
  void Load( const char *FileName );
  ///
  Matrix%TYPE2% MyTransposeTimesMe() const;
  ///
  Matrix%TYPE2% RowMatrix() const;
  ///
  Matrix%TYPE2% ColumnMatrix() const;
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
  void RemoveElement( int );
  ///
  void AddElement( int );
  ///
  gsl_vector%GSL% *GetGslVector();

};

#endif /* !MATRIX_VECTOR_D_H */
