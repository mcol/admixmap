// *-*-C++-*-*
#ifndef DIPLOID_TRANSITION_MATRIX_H
#define DIPLOID_TRANSITION_MATRIX_H 1

#include "HaploidTransitionMatrix.h"
#include "HaploidProbabilities.h"
#include "Haploid.h"

#include <cassert>

#include "matrix_d.h"

class DiploidTransitionMatrix
{
private: //members
  unsigned int _D;
  double **_T;

  // static helpers
  static unsigned int delta(const unsigned int, const unsigned int);
  
  // UNIMPLIMENTED
  // - to avoid use
  DiploidTransitionMatrix();
  DiploidTransitionMatrix(const DiploidTransitionMatrix&);
  DiploidTransitionMatrix& operator=(const DiploidTransitionMatrix&);

public:

  // CONSTRUCTOR
  DiploidTransitionMatrix(HaploidTransitionMatrix&,HaploidTransitionMatrix&,Vector_d&);

  // DESTRUCTOR
  ~DiploidTransitionMatrix();

  unsigned int size();
  double operator()(const unsigned int,const unsigned int);
  Matrix_d toMatrix();

};

#endif /* !defined DIPLOID_TRANSITION_MATRIX_H */
