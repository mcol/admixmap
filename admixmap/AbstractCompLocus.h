// *-*-C++-*-*
#ifndef ABSTRACT_COMP_LOCUS_H
#define ABSTRACT_COMP_LOCUS_H 1

#include <string>
#include "matrix.h"
#include "vector_d.h"
#include "matrix_d.h"
#include "vector_i.h"
#include "matrix_i.h"
#include "AdmixOptions.h"
#include "MatrixArray_i.h"
#include <vector>
#include <assert.h>
//#include "AlleleFreqs.h"

class Individual;
class AlleleFreqs;
class AbstractCompLocus
{
public:
  virtual ~AbstractCompLocus() {}

  virtual int GetNumberOfLoci() = 0;
  
  virtual int GetNumberOfStates() = 0 ;
  
  virtual int GetSize() = 0;

  virtual void SetLabel( int, std::string ) = 0;
  
};

#endif /* !ABSTRACT_COMP_LOCUS_H */
