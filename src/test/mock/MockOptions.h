#ifndef MOCKOPTIONS_H_
#define MOCKOPTIONS_H_

//#include <exception>
//#include <iostream>

#define MOCKPP_IMPORT_ABBREVIATED
#include <mockpp/mockpp.h>
#include <mockpp/visiting/VisitableMockObject.h>
#include <mockpp/chaining/ChainingMockObjectSupport.h>
#include "../../admixmap/interfaces/IOptions.h"

USING_NAMESPACE_MOCKPP

class MockOptions : public VisitableMockObject,
                    public IOptions
{
public:
	MockOptions()
  : VisitableMockObject(MOCKPP_PCHAR("MockOptions"), 0)
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(getHWTestIndicator)
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(getgenotypesSexColumn)
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(isRandomMatingModel)
//  virtual int getPopulations() const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(getPopulations)
  {}
  virtual ~MockOptions();

  MOCKPP_CONST_VISITABLE0(MockOptions, bool, getHWTestIndicator);
  MOCKPP_CONST_VISITABLE0(MockOptions, unsigned, getgenotypesSexColumn);
  MOCKPP_CONST_VISITABLE0(MockOptions, bool, isRandomMatingModel);
  // bool getHWTestIndicator() const;
  // unsigned int getgenotypesSexColumn() const;
  // bool isRandomMatingModel() const;
  // virtual int getPopulations() const = 0;
  MOCKPP_CONST_VISITABLE0(MockOptions, int, getPopulations);
};

#endif /*MOCKOPTIONS_H_*/
