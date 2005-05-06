// *-*-C++-*-*
#ifndef DARS_H
#define DARS_H 1

#include "vector_d.h"
#include "MatrixArray_i.h"
#include "MatrixArray_d.h"
#include "rand.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <stdlib.h>


class DARS
{
private: // members
  Vector_d parameters;
  int LeftFlag;
  int RightFlag;
  int n;
  int no;
  int lgth;
  int loc;
  MatrixArray_i data_i;
  MatrixArray_d data_d;
  double x0;
  double x1;
  double *x;
  double *f;
  double *df;
  double *psum;
  double *u;
  double *z;
  double newnum;
  double (*function)(Vector_d &, MatrixArray_i &, MatrixArray_d &, double xin);
  double (*dfunction)(Vector_d &, MatrixArray_i &, MatrixArray_d &, double xin);
  double (*ddfunction)(Vector_d &, MatrixArray_i &, MatrixArray_d &, double xin);

  // UNIMPLEMENTED
  // to avoid use
  DARS(const DARS&);
  DARS& operator=(const DARS&);

  // PRIVATE HELPERS
  static double
  fannyexp( double );

public:

  // CONSTRUCTORS

  // default constructor
  DARS();
  // all arguments
  DARS( int, int, double newnum , const double inparameters[],int size,//const Vector_d &inparameters,
        double (*funct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
        double (*dfunct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
        double (*ddfunct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
	const MatrixArray_i&, const MatrixArray_d&);

  // DESTRUCTOR
  ~DARS();
  
  double
  Sample();
  
  double SampleUsingARS();
  
  void SetParameters(int, int, double newnum, const double inparameters[],int size,//const Vector_d &inparameters,
                double (*funct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
                double (*dfunct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
                double (*ddfunct)(Vector_d&, MatrixArray_i&, MatrixArray_d&, double),
		const MatrixArray_i&, const MatrixArray_d&);
  
  void SetLeftTruncation(double);
  
  void SetRightTruncation( double );
  
  //void
  //UpdateParameters(const Vector_d &inparameters);
  void UpdateParameters( const double inparameters[], int size );
  
  void
  UpdateIntegerData(const MatrixArray_i&);
  
  void
  UpdateDoubleData(const MatrixArray_d&);
  
  void
  BeginModeSearch(double);
  
  void
  conspsum();
  
  void
  consz();
  
  void
  SimpleModeSearch(double, double);
  
  void
  NewtonRaphson();
  
};


#endif /* !defined DARS_H */
