// *-*-C++-*-*
#ifndef INDIVIDUAL_VISITOR_H
#define INDIVIDUAL_VISITOR_H 1
#include <vector>

class Individual;
class IndividualCollection;

class IndividualVisitor
{
public:
  
  virtual void
  visitIndividual(Individual&,double, std::vector<int>,double) = 0;
  
  virtual void
  visitIndividualCollection(IndividualCollection&) = 0;

  virtual ~IndividualVisitor(){}
};

#endif /* !defined INDIVIDUAL_VISITOR_H */
