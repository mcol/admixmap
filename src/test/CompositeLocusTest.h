#ifndef COMPOSITELOCUSTEST_H_
#define COMPOSITELOCUSTEST_H_

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/Exception.h>

#include "../admixmap/CompositeLocus.h"

class CompositeLocusTest : public CppUnit::TestFixture
{
private:
  CompositeLocus *cl1;

  CPPUNIT_TEST_SUITE( CompositeLocusTest );
  CPPUNIT_TEST( testStub );
  CPPUNIT_TEST( testAddLocus );
  CPPUNIT_TEST_SUITE_END();
public:
	CompositeLocusTest();
	virtual ~CompositeLocusTest();
  void setUp();
  void tearDown();
  
  void testStub();
  void testAddLocus();
};

#endif /*COMPOSITELOCUSTEST_H_*/
