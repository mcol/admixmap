#ifndef HAPMIXINDIVIDUALTEST_H_
#define HAPMIXINDIVIDUALTEST_H_

#include <math.h>

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/Exception.h>
#include "src/admixmap/HapMixIndividual.h"
#include "mock/MockOptions.h"
#include "mock/MockInputData.h"
#include "mock/MockHapMixFreqs.h"
#include "mock/MockGenome.h"
#include "mock/MockFreqArray.h"
#include "mock/convenience-macros.h"

#include <vector>
using std::vector;

class HapMixIndividualTest : public CppUnit::TestFixture
{
private:
  HapMixIndividual *testIndiv1;
  MockOptions *options;
  MockInputData *inputData;
  MockGenome *genome;
  // MockHapMixFreqs *A; // not needed really
  MockFreqArray *haploidGenProbs, *diploidGenProbs;
  
  CPPUNIT_TEST_SUITE( HapMixIndividualTest );
  CPPUNIT_TEST( testUnorderedGenotypeProbs );
  CPPUNIT_TEST_SUITE_END();
public:
  HapMixIndividualTest();
  virtual ~HapMixIndividualTest();
  void setUp();
  void tearDown();
  
  // Tests
  void testUnorderedGenotypeProbs();
};

#endif /*HAPMIXINDIVIDUALTEST_H_*/
