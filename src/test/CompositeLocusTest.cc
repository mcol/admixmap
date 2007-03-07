#include "CompositeLocusTest.h"

CompositeLocusTest::CompositeLocusTest()
{
}

CompositeLocusTest::~CompositeLocusTest()
{
}

void CompositeLocusTest::setUp()
{
  cl1 = new CompositeLocus();
}

void CompositeLocusTest::tearDown()
{
  delete cl1;
}

void CompositeLocusTest::testStub()
{
  CPPUNIT_ASSERT( true );
}

void CompositeLocusTest::testAddLocus()
{
  // TODO: Test if the loci were really added.
  cl1->AddLocus(2, "rs01");
  cl1->AddLocus(2, "rs02");
  cl1->AddLocus(2, "rs03");
}
