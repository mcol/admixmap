// *-*-C++-*-*
#ifndef DARS_H
#define DARS_H 1

#include "vector_d.h"
#include "Matrix_i.h"
#include "Matrix_d.h"
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
  Matrix_i data_i;
  Matrix_d data_d;
  double x0;
  double x1;
  double *x;
  double *f;
  double *df;
  double *psum;
  double *u;
  double *z;
  double newnum;
  double (*function)(Vector_d &, Matrix_i &, Matrix_d &, double xin);
  double (*dfunction)(Vector_d &, Matrix_i &, Matrix_d &, double xin);
  double (*ddfunction)(Vector_d &, Matrix_i &, Matrix_d &, double xin);

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
        double (*funct)(Vector_d&, Matrix_i&, Matrix_d&, double),
        double (*dfunct)(Vector_d&, Matrix_i&, Matrix_d&, double),
        double (*ddfunct)(Vector_d&, Matrix_i&, Matrix_d&, double),
	const Matrix_i &, const Matrix_d &);

  // DESTRUCTOR
  ~DARS();
  
  double
  Sample();
  
  double SampleUsingARS();
  
  void SetParameters(int, int, double newnum, const double inparameters[],int size,//const Vector_d &inparameters,
                double (*funct)(Vector_d&, Matrix_i&, Matrix_d&, double),
                double (*dfunct)(Vector_d&, Matrix_i&, Matrix_d&, double),
                double (*ddfunct)(Vector_d&, Matrix_i&, Matrix_d&, double),
		const Matrix_i &, const Matrix_d &);
  
  void SetLeftTruncation(double);
  
  void SetRightTruncation( double );
  
  //void
  //UpdateParameters(const Vector_d &inparameters);
  void UpdateParameters( const double inparameters[], int size );
  
  void
  UpdateIntegerData(const Matrix_i&);
  
  void
  UpdateDoubleData(const Matrix_d&);
  
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
