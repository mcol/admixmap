#ifndef MATRIX_VECTOR_D_H
#define MATRIX_VECTOR_D_H 1

#include <iostream>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cfloat>

#include <gsl/gsl_vector.h>

class Matrix_d;
class Vector;
class Vector_i;

#ifndef LINELENGTH
#define LINELENGTH 300000
#endif

#ifndef MAX_LIN_ACTIVITY
#define MAX_LIN_ACTIVITY 700
#endif

#ifndef MIN_LIN_ACTIVITY 
#define MIN_LIN_ACTIVITY -MAX_LIN_ACTIVITY
#endif

// vector of double-precision
// floating-point numbers
class Vector_d
{
 private: //members
  int NumberOfElements;
  int NumberOfMissingValues;
  int *MissingValues;
  gsl_vector *vector;

  // FRIENDS
  friend std::ostream& operator<<(std::ostream& os,const Vector_d& V);
  friend std::istream& operator>>(std::istream& is,Vector_d& V);

 public:
  // CONSTRUCTORS

  // default constructor
  Vector_d();

  // copy constructor
  Vector_d(const Vector_d& V);

  // initial length
  Vector_d(int NewNumberOfElements);

  // DESTRUCTOR
  ~Vector_d();

  double&
    operator()(int);

  int
    GetNumberOfMissingValues() const;

  int
    GetNumberOfElements() const;

  int
    IsMissingElement(int index) const;

  Vector_d& 
    operator/=(double f);

  Vector_d&
    operator=(const Vector_d& v);

  void
    SetNumberOfElements(int NewValue);

  void
    SetElements(double Value);

  double
    Sum() const;

   double
     Mean() const;

  Vector_d
    operator+(const Vector_d &v) const;

  Vector_d
    operator+(const Vector_i &v) const;

  Vector_d&
    operator+=(const Vector_d &v);

  Vector_d
    operator/(double f) const;

  double
    GetElement(int index) const;

  double&
    operator()(int) const;

  Vector_d
    operator-(const Vector_d &v) const;

  Matrix_d
    ColumnMatrix() const;

  double
    MaximumElement(void);

  void
    AddElement(int);

  void
    RemoveElement(int);

  Vector
    Float();

  Vector_d
    operator*(double f) const;

  Vector_d
    operator/(const Vector_d& v) const;

private:
  int
    GetMissingValueIndex(int MissingValueNumber) const;
};

#endif /* !MATRIX_VECTOR_D_H */
