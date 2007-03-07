#include "MockInputData.h"

MockInputData::~MockInputData()
{
  delete GetGenotypeDefaultReturnValue;
  delete GetHapMixGenotypeDefaultReturnValue1;
  delete GetHapMixGenotypeDefaultReturnValue2;
}

void MockInputData::GetGenotype(
      int indNumber,
      int sexColumn,
      const IGenome& Loci /* or genome?*/,
      std::vector<genotype>* genotypes,
      bool ** genotypesMissing) const
{
  THROW_IF_EMPTY(GetGenotypeDefaultReturnValue)
  RETURN_VEC_BY_PARAM(GetGenotypeDefaultReturnValue, genotypes);
  return;
}

bool MockInputData::GetHapMixGenotype(
      int i,
      int SexColumn,
      const IGenome &Loci,
      std::vector<unsigned short>* genotypes,
      bool** Missing) const
{
//  genotypes->push_back(2);
  THROW_IF_EMPTY(GetHapMixGenotypeDefaultReturnValue1)
  THROW_IF_EMPTY(GetHapMixGenotypeDefaultReturnValue2)
  RETURN_VEC_BY_PARAM(GetHapMixGenotypeDefaultReturnValue2, genotypes);
  // goes into (HapMix)Individual::isHaploid
  return *GetHapMixGenotypeDefaultReturnValue1;
}

void MockInputData::setGetGenotypeDefaultReturnValue(std::vector<genotype> v)
{
  GetGenotypeDefaultReturnValue = new std::vector<genotype>;
//  copy(GetGenotypeDefaultReturnValue->begin(),
//      GetGenotypeDefaultReturnValue->end(),
//      v.begin());
  GetGenotypeDefaultReturnValue->assign(v.begin(), v.end());
}

void MockInputData::setGetHapMixGenotypeDefaultReturnValues(bool b, std::vector<unsigned short> v)
{
  GetHapMixGenotypeDefaultReturnValue1 = new bool;
  *GetHapMixGenotypeDefaultReturnValue1 = b;
  GetHapMixGenotypeDefaultReturnValue2 = new std::vector<unsigned short>;
  GetHapMixGenotypeDefaultReturnValue2->assign(v.begin(), v.end());
}
