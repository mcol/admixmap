#include <cstdlib>
// #include <iostream>
#include <cppunit/TestSuite.h>
// #include <cppunit/TestResult.h>
// #include <cppunit/TestCaller.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/extensions/TestFactoryRegistry.h>

#include "UserDataLocusTest.h"
#include "LegendTest.h"
#include "HapMixIndividualTest.h"
#include "CompositeLocusTest.h"

CPPUNIT_TEST_SUITE_REGISTRATION(UserDataLocusTest);
CPPUNIT_TEST_SUITE_REGISTRATION(LegendTest);
CPPUNIT_TEST_SUITE_REGISTRATION(HapMixIndividualTest);
CPPUNIT_TEST_SUITE_REGISTRATION(CompositeLocusTest);

int main(int arg, char ** argv)
{
  CppUnit::TextUi::TestRunner runner;
  CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();
  runner.addTest( registry.makeTest() );
  bool wasSuccessful = runner.run( "", false );
  return wasSuccessful ? EXIT_SUCCESS : EXIT_FAILURE;
}
