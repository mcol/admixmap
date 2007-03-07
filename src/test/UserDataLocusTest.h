#ifndef USERDATALOCUSTEST_H_
#define USERDATALOCUSTEST_H_

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include "../tools/UserDataLocus.h"
#include "../tools/Legend.h"

class UserDataLocusTest : public CppUnit::TestFixture
{
private:
  Legend *legend;
  UserDataLocus *locus1, *locus2, *locus3, *locus4, *locus5;
  CPPUNIT_TEST_SUITE( UserDataLocusTest );
  CPPUNIT_TEST( testLessSameChromosome );
  CPPUNIT_TEST( testGreaterSameChromosome );
  CPPUNIT_TEST( testLessOtherChromosome );
  CPPUNIT_TEST( testGreaterOtherChromosome );
  CPPUNIT_TEST( testEqual );
  CPPUNIT_TEST( testNotEqual );
  CPPUNIT_TEST( testLocus1 );
  CPPUNIT_TEST( testLocus2and3 );
  CPPUNIT_TEST_SUITE_END();
public:
	UserDataLocusTest();
	virtual ~UserDataLocusTest();
  void setUp();
  void tearDown();
  void testLessSameChromosome();
  void testGreaterSameChromosome();
  void testLessOtherChromosome();
  void testGreaterOtherChromosome();
  void testEqual();
  void testNotEqual();
  void testLocus1();
  void testLocus2and3();
};

#endif /*USERDATALOCUSTEST_H_*/
