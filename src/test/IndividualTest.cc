#include "IndividualTest.h"

IndividualTest::IndividualTest()
{
}

IndividualTest::~IndividualTest()
{
}

void IndividualTest::setUp()
{
  testIndiv1 = new Individual();
  return;
}

void IndividualTest::tearDown()
{
  delete testIndiv1;
  return;
}

void IndividualTest::testStub()
{
  return;
}
