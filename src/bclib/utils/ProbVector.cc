#include "bclib/ProbVector.h"
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>
#include "gsl/gsl_vector.h"

#define PVECTOR_PRECISION 1e-10
#define PVECTOR_SUM 1.0
#define PVECTOR_THRESHOLD 1e-20

using namespace::std;

BEGIN_BCLIB_NAMESPACE

//auxiliary functions

///use to set a number below a given threshold to 0
double SnapToZero(double x){
    if (x < PVECTOR_THRESHOLD)  return 0;
    else return x; 
}
///determine if a number is non-negative
bool IsNonNegative(double x){
  return (x >=0.0);
}

///checks elements sum to 1
bool ProbVector::is_normalized()
{
  // TODO: Accumulate only when necessary, i.e.
  // vector was changed.
  sum = accumulate(this->begin(), this->end(), 0.0);
  return (fabs(sum - PVECTOR_SUM) < PVECTOR_PRECISION);
}

///make probabilities sum to 1 
void ProbVector::normalize()
{
  if (!this->is_normalized()) { // Updates sum
    gsl_vector_view V = gsl_vector_view_array( &(*this)[0], this->size());
    gsl_vector_scale(&V.vector, 1.0/sum);
  }
}

///test that elements are nonnegative and sum to one. No need to check they are <= 1 as that is implied.
bool ProbVector::verify()
{

  if(count_if(this->begin(), this->end(), IsNonNegative) != (int)this->size()) {
    //throw string("pvector::verify(): Some elements are negative.");
    return false;
  }
  if (!is_normalized()) {
    //throw string("pvector::verify(): doesn't sum to one.");
    return false;
  }
  return true;
}

///sets 'small' and negative elements to zero
void ProbVector::snapToZero()
{
  transform(this->begin(), this->end(), this->begin(), SnapToZero);
  this->normalize();
}

///print to an output stream
void ProbVector::print(ostream& os){

  copy(this->begin(), this->end(), ostream_iterator<double>(os, " " ) );
  os << endl;
}
END_BCLIB_NAMESPACE
