// *-*-C++-*-*
//RevMCMC.h

#ifndef HMM_H
#define HMM_H 1

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <values.h>

#define FILENAMELENGTH 100

class Vector_i;
class Matrix_d;
class MatrixArray_d;

class HMM
{
public:
  HMM();
  HMM( int, int );
  ~HMM();
   void SetDimensions( int, int );
  Vector_i Sample();
  void UpdateParameters( Matrix_d&, MatrixArray_d&, MatrixArray_d&, bool );
   Vector_d GetStateProbs( int );
   double getLikelihood();
   
private:
  void CheckArguments();
  void UpdateFwrdBckwdProbailities();

  int States;
  int Transitions;
   bool _CalculateBeta;
   double factor;

  Matrix_d StationaryDist;
  MatrixArray_d TransitionProbs;
  MatrixArray_d Likelihood;
   MatrixArray_d alpha;
   MatrixArray_d beta;
};

#endif /* ! HMM_H */
