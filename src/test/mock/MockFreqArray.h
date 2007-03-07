#ifndef MOCKFREQARRAY_H_
#define MOCKFREQARRAY_H_

#include <string>
using std::string;

#define MOCKPP_IMPORT_ABBREVIATED
#include <mockpp/mockpp.h>
#include <mockpp/visiting/VisitableMockObject.h>
#include <mockpp/chaining/ChainingMockObjectSupport.h>
#include "../../admixmap/interfaces/IFreqArray.h"

USING_NAMESPACE_MOCKPP

class MockFreqArray : public VisitableMockObject,
                      public IFreqArray
{
public:
	MockFreqArray()
  : VisitableMockObject(MOCKPP_PCHAR("MockFreqArray"), 0)
  // void dealloc(int) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT1(dealloc, ext)
  {}
	virtual ~MockFreqArray();
  // void dealloc(int) = 0;
  MOCKPP_VOID_VISITABLE_EXT1(MockFreqArray, dealloc, int,
      ext, int);
      
  /*
   * Following methods are unimplemented.
   */
  double* operator[](unsigned int) {throw string("unimplemented");}
  const double* operator[](unsigned int) const {throw string("unimplemented");}
};

#endif /*MOCKFREQARRAY_H_*/
