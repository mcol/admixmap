// *-*-C++-*-*
#ifndef HAPLOID_TRANSITION_MATRIX_H
#define HAPLOID_TRANSITION_MATRIX_H 1

#include <cassert>

#include "matrix_d.h"
#include "vector_d.h"
//#include "GameteAdmixture.h"

/**
 * The HaploidTransitionMatrix Class.
 *
 * The haploid transition matrix represents the 
 * probability of one haploid state existing at the second
 * time point (t+1), given a particular state at
 * the first time point (t).
 *
 * Generate a HaploidTransitionMatrix from a GameteAdmixture.
 *
 *   double theta = 0.68;
 *   HaploidTransitionMatrix htm(ga,theta);
 *
 * Use the toMatrix() method to convert to a square matrix
 * of doubles.
 */

// David - this class should be edited so that it you can just set elements of the 
// haploid transition matrix without allocating memory for a new object  
class HaploidTransitionMatrix
{
private: // members
  unsigned int _H;
  double **_T;

  // UNIMPLEMENTED
  // - to prevent use
  HaploidTransitionMatrix();
  HaploidTransitionMatrix(const HaploidTransitionMatrix&);
  HaploidTransitionMatrix& operator=(const HaploidTransitionMatrix&);

public:

  // CONSTRUCTOR
  HaploidTransitionMatrix(Vector_d& ,double);

  // DESTRUCTOR
  ~HaploidTransitionMatrix();

  unsigned int size();
  double& operator()(const unsigned int,const unsigned int);
  Matrix_d toMatrix();
};

#endif /* !defined HAPLOID_TRANSITION_MATRIX_H */
