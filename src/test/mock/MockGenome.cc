#include "MockGenome.h"

//MockGenome::MockGenome()
//{
//}

MockGenome::~MockGenome()
{
  if (chromosome) {
    chromosome->verify();
    delete chromosome;
  }
  if (cl) {
    cl->verify();
    delete cl;
  }
}

void MockGenome::setMockChromosome(MockChromosome *c)
{
  chromosome = c;
}

void MockGenome::setMockCompositeLocus(MockCompositeLocus* cl_tmp)
{
  cl = cl_tmp;
}

ICompositeLocus* MockGenome::operator()(int) const
{
  THROW_IF_EMPTY(cl);
  return cl;
}

IChromosome* MockGenome::getChromosome(unsigned int)
{
  if (not chromosome) {
    throw CppUnit::Exception(CppUnit::Message(
        "IChromosome* MockGenome::getChromosome(): chromosome not set"),
        CPPUNIT_SOURCELINE());
  }
  return chromosome;
}
