#include "HaploidTransitionMatrix.h"

using namespace std;

HaploidTransitionMatrix::HaploidTransitionMatrix
(Vector_d& theta,double f):
  _H(theta.GetNumberOfElements()),
  _T(new double*[theta.GetNumberOfElements()])
{
  double _product[_H];
  for(unsigned int i=0; i<_H; i++){
    _T[i] = new double[_H];
    _product[i] = f * theta(i);
  }

  for(unsigned int i=0; i<_H; i++){
    for(unsigned int j=0; j<i; j++){
      _T[i][j] = theta(j) - _product[j];
    }
    _T[i][i] = theta(i) + f - _product[i];
    for(unsigned int j=i+1; j<_H; j++){
      _T[i][j] = theta(j) - _product[j];
    }
  }
}

HaploidTransitionMatrix::~HaploidTransitionMatrix()
{
  for(unsigned int i=0; i<_H; i++){
    delete _T[i];
  }

  delete _T;
}

double&
HaploidTransitionMatrix::operator()(const unsigned int i, const unsigned int j)
{
  assert(i<_H);
  assert(j<_H);
  return _T[i][j];
}

unsigned int
HaploidTransitionMatrix::size()
{
  return _H;
}

Matrix_d
HaploidTransitionMatrix::toMatrix()
{
  Matrix_d result(_H,_H);
  for(unsigned int i=0;i<_H;i++){
    for(unsigned int j=0;j<_H;j++){
      result(i,j) = _T[i][j];
    }
  }
  return result;
}
