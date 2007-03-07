#ifndef INDIVIDUALTEST_H_
#define INDIVIDUALTEST_H_

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../admixmap/Individual.h"

class IndividualTest : public CppUnit::TestFixture
{
private:
  Individual *testIndiv1;
  CPPUNIT_TEST_SUITE( IndividualTest );
  CPPUNIT_TEST( testStub );
  CPPUNIT_TEST_SUITE_END();
public:
	IndividualTest();
	virtual ~IndividualTest();
  void setUp();
  void tearDown();
  void testStub();
};

#endif /*INDIVIDUALTEST_H_*/
