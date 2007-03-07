#include "LegendTest.h"

LegendTest::LegendTest()
{
}

LegendTest::~LegendTest()
{
}

void LegendTest::setUp()
{
  #include "LegendData.h"
  legend1 = new Legend(legend_data);
}

void LegendTest::tearDown()
{
  delete legend1;
}

void LegendTest::testLocusStruct()
{
  locus_t *locus1 = legend1->getLocusPointerBySnp("rs1000013");
  CPPUNIT_ASSERT(locus1->chromosome == 7);
  CPPUNIT_ASSERT(locus1->position == 40363809);
  CPPUNIT_ASSERT(locus1->a1 == 'C');
  CPPUNIT_ASSERT(locus1->a2 == 'T');
}

void LegendTest::testInitialized()
{
  CPPUNIT_ASSERT(legend1->isInitialized() == true);
}
