// *-*-C++-*-*
#ifndef INDIVIDUAL_VISITOR_H
#define INDIVIDUAL_VISITOR_H 1
#include <vector>
#include "matrix_d.h"

class Individual;
class IndividualCollection;

class IndividualVisitor
{
public:
  
  virtual void
  visitIndividual(Individual&,std::vector<int>,double, Matrix_d &) = 0;
  
  virtual void
  visitIndividualCollection(IndividualCollection&) = 0;

  virtual ~IndividualVisitor(){}
};

#endif /* !defined INDIVIDUAL_VISITOR_H */
