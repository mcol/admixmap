// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   DARS.h
 *   header file for DARS class
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef DARS_H
#define DARS_H 1

#include "matrix_i.h"
#include "matrix_d.h"
#include "rand.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <string.h>
#include <stdlib.h>


class DARS
{
public:

  // CONSTRUCTORS
  
  // default constructor
  DARS();
  // all arguments
  DARS( int, int, double newnum , const double inparameters[],
        double (*funct)(const double *, Matrix_i&, Matrix_d&, double),
        double (*dfunct)(const double*, Matrix_i&, Matrix_d&, double),
        double (*ddfunct)(const double*, Matrix_i&, Matrix_d&, double),
	const Matrix_i &, const Matrix_d &);

  // DESTRUCTOR
  ~DARS();
  
  double Sample();
  
  double SampleUsingARS();
  
  void SetParameters(int, int, double newnum, const double inparameters[],
                double (*funct)(const double*, Matrix_i&, Matrix_d&, double),
                double (*dfunct)(const double*, Matrix_i&, Matrix_d&, double),
                double (*ddfunct)(const double*, Matrix_i&, Matrix_d&, double),
		const Matrix_i &, const Matrix_d &);
  
  void SetLeftTruncation(double);
  
  void SetRightTruncation( double );
  
  void UpdateParameters( const double inparameters[]);
  
  void UpdateIntegerData(const Matrix_i&);
  
  void UpdateDoubleData(const Matrix_d&);
  
  void BeginModeSearch(double);
  
  void conspsum();
  
  void consz();
  
  void SimpleModeSearch(double, double);
  
  void NewtonRaphson();

private:
  const double* parameters;
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
  double (*function)(const double*, Matrix_i &, Matrix_d &, double xin);
  double (*dfunction)(const double*, Matrix_i &, Matrix_d &, double xin);
  double (*ddfunction)(const double*, Matrix_i &, Matrix_d &, double xin);

  // UNIMPLEMENTED
  // to avoid use
  DARS(const DARS&);
  DARS& operator=(const DARS&);

  // PRIVATE HELPERS
  static double fannyexp( double );
  
};


#endif /* !defined DARS_H */
