#ifndef LEGENDTEST_H_
#define LEGENDTEST_H_

#include <fstream>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../tools/Legend.h"

class LegendTest : public CppUnit::TestFixture
{
private:
  Legend *legend1;
  CPPUNIT_TEST_SUITE( LegendTest );
  CPPUNIT_TEST( testLocusStruct );
  CPPUNIT_TEST_SUITE_END();
public:
	LegendTest();
	virtual ~LegendTest();
  void setUp();
  void tearDown();
  void testLocusStruct();
  void testInitialized();
};

#endif /*LEGENDTEST_H_*/
