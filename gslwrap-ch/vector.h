#ifndef MATRIX_VECTOR_H
#define MATRIX_VECTOR_H 1

#if !defined( FLOAT_ACC )
#define FLOAT_ACC float
#endif

#if !defined( LINELENGTH )
#define LINELENGTH 300000
#endif

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cassert>

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

// partial declarations here, because
// of circular dependancies
class Matrix;
class Vector_d;
class Vector_i;

extern "C"{
  //? why C linked?
  int CompareVector( const void *, const void * );
}

// generic vector
class Vector
{
 private: // members
  int NumberOfElements;
  int NumberOfMissingValues;
  int *MissingValues;
  gsl_vector_float *vector;

 public:
  
  // CONSTRUCTORS
  
  // default constructor
  Vector();
  // copy constructor
  Vector(const Vector& V);
  // length
  Vector(int NewNumberOfElements);
  // load from file
  Vector(const char* FileName);

  // DESTRUCTOR
  ~Vector();

  ///
  Vector_d Double();
  ///
  Vector_i Integer();
  ///
  FLOAT_ACC GetElement( int index ) const;
  ///
  float& operator()(int);
  float& operator()(int) const;
  float& operator[](int);
  float& operator[](int) const;
  ///
  void SetElement( int index, FLOAT_ACC NewValue );
  ///
  int GetNumberOfElements() const;
  ///
  void SetNumberOfElements( int NewValue );
  ///
  friend std::ostream& operator<< ( std::ostream& os, const Vector& V );
  ///
  friend std::istream& operator>> ( std::istream& is, Vector& V );
  ///
  Vector operator+( const Vector &v ) const;
  ///
  Vector operator+( double f ) const;
  ///
  friend Vector operator+( double f, const Vector &v );
  ///
  Vector &operator+=( double f );
  ///
  Vector &operator+=( const Vector &v );
  ///
  Vector operator-( const Vector &v ) const;
  ///
  Vector operator-( double f ) const;
  ///
  friend Vector operator-( double f, const Vector &v );
  ///
  Vector &operator-=( double f );
  ///
  Vector &operator-=( const Vector &v );
  ///
  Vector operator-() const;
  ///
  Vector operator*( const Vector &v ) const;
  ///
  Vector operator*( double f ) const;
  ///
  friend Vector operator*( double f, const Vector &v );
  ///
  Vector &operator*=( double f );
  ///
  Vector &operator*=( const Vector &v );
  ///
  Vector operator/( const Vector &v ) const;
  ///
  Vector operator/( double f ) const;
  ///
  friend Vector operator/( double f, const Vector &v );
  ///
  Vector &operator/=( double f );
  ///
  Vector &operator/=( const Vector &v );
  ///
  Vector &operator=( const Vector &v );
  ///
  int operator==( const Vector &v ) const;
  ///
  void Randomize( FLOAT_ACC LowerLimit, FLOAT_ACC UpperLimit );
  ///
  void RandomizeGaussian(float, float);
  ///
  void Sigmoid(void);
  ///
  void Exp(void);
  ///
  void Tanh(void);
  ///
  void SetElements( FLOAT_ACC Value );
  ///
  double Sum() const;
  ///
  double Mean() const;
  ///
  double Variance(void) const;
  ///
  Vector Absolute(void) const;
  ///
  float Norm();
  ///
  float MinimumElement(void);
  ///
  float MaximumElement(void);
  ///
  int WhereMinimumElement(void);
  ///
  int WhereMaximumElement(void);
  ///
  float MinimumAbsoluteElement(void);
  ///
  float MaximumAbsoluteElement(void);
  ///
  void Load( const char *FileName );
  ///
  Matrix MyTransposeTimesMe() const;
  ///
  Matrix RowMatrix() const;
  ///
  Matrix ColumnMatrix() const;
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
  void Sort();
  ///
  gsl_vector_float *GetGslVector();
};

#endif /* !MATRIX_VECTOR_H */
