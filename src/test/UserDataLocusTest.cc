#include "UserDataLocusTest.h"

UserDataLocusTest::UserDataLocusTest()
{
}

UserDataLocusTest::~UserDataLocusTest()
{
}

void UserDataLocusTest::setUp()
{
  #include "LegendData.h"
  legend = new Legend(legend_data);
  locus1 = new UserDataLocus(0, legend->getLocusPointerBySnp("rs1000013"));
  locus2 = new UserDataLocus(1, legend->getLocusPointerBySnp("rs1000058"));
  locus3 = new UserDataLocus(2, legend->getLocusPointerBySnp("rs1000059"));
  locus4 = new UserDataLocus(3, legend->getLocusPointerBySnp("rs1000059"));
  locus5 = new UserDataLocus(4, legend->getLocusPointerBySnp("rs10000007"));
}

void UserDataLocusTest::tearDown()
{
  delete locus1;
  delete locus2;
  delete locus3;
  delete locus4;
  delete locus5;
}

void UserDataLocusTest::testLessSameChromosome()
{
  CPPUNIT_ASSERT(locus3->getChromosome() == locus2->getChromosome());
  CPPUNIT_ASSERT(locus3->getPosition() < locus2->getPosition());
  CPPUNIT_ASSERT(*locus3 < *locus2);
}

void UserDataLocusTest::testGreaterSameChromosome()
{
  CPPUNIT_ASSERT(locus3->getChromosome() == locus2->getChromosome());
  CPPUNIT_ASSERT(locus2->getPosition() > locus3->getPosition());
  CPPUNIT_ASSERT(*locus2 > *locus3);
}

void UserDataLocusTest::testLessOtherChromosome()
{
  CPPUNIT_ASSERT(locus5->getChromosome() < locus2->getChromosome());
  CPPUNIT_ASSERT(*locus5 < *locus2);
}

void UserDataLocusTest::testGreaterOtherChromosome()
{
  CPPUNIT_ASSERT(locus2->getChromosome() > locus5->getChromosome());
  CPPUNIT_ASSERT(*locus2 > *locus5);
}

void UserDataLocusTest::testEqual()
{
  CPPUNIT_ASSERT(*locus2 == *locus2);
  CPPUNIT_ASSERT(*locus3 == *locus4);
}

void UserDataLocusTest::testNotEqual()
{
  CPPUNIT_ASSERT(*locus2 != *locus5);
  CPPUNIT_ASSERT(*locus1 != *locus2);
}

void UserDataLocusTest::testLocus1()
{
  CPPUNIT_ASSERT(locus1->getChromosome() == 7);
  CPPUNIT_ASSERT(locus1->getPosition() == 40363809);
  CPPUNIT_ASSERT(locus1->getAllele1() == 'C');
  CPPUNIT_ASSERT(locus1->getAllele2() == 'T');
}

void UserDataLocusTest::testLocus2and3()
{
  CPPUNIT_ASSERT(locus2->getChromosome() == 7);
  CPPUNIT_ASSERT(locus2->getPosition() == 131554271);
  CPPUNIT_ASSERT(locus3->getChromosome() == 7);
  CPPUNIT_ASSERT(locus3->getPosition() == 131554260);
}
