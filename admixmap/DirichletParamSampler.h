// *-*-C++-*-*
#ifndef DirichletParamSampler_H
#define DirichletParamSampler_H 1

#include <gsl/gsl_sf_gamma.h>
#include "matrix_i.h"
#include "matrix_d.h"
#include "vector_d.h"
#include "DARS.h"
#include "TuneRW.h"
#include "rand.h"

class DirichletParamSampler
{
public:
  DirichletParamSampler();
  DirichletParamSampler(unsigned int);
  ~DirichletParamSampler();
  
  void SetSize( unsigned int );
  void SetPriorEta( double, double );
  void SetPriorMu( double* );
  void Sample( unsigned int, double*, double*, double* );
  
private:
  TuneRW TuneEta;
  unsigned int d;
  double etanew;
  double *munew;
  double *gamma;
  double EtaAlpha;
  double EtaBeta;

  DARS** DirParamArray;
  // AlphaParameters is an array with 5 elements
  // element 0 is 
  // element 1 is 
  // 
  // element 4 is 
  double AlphaParameters[5];
  static double
  logf( Vector_d&, Matrix_i&, Matrix_d& , double );
  
  static double
  dlogf( Vector_d&, Matrix_i&, Matrix_d& , double );
  
  static double
  ddlogf( Vector_d&, Matrix_i&, Matrix_d& , double );
};

#endif /* ! DirichletParamSampler_H */
