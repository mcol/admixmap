#include "MockCompositeLocus.h"

MockCompositeLocus::~MockCompositeLocus()
{
}

void MockCompositeLocus::getConditionalHapPairProbs(std::vector<double>& Probs,
      const std::vector<hapPair> &PossibleHapPairs,
      const int ancestry[2])const
{
//  Probs.clear();
//  Probs.push_back(0.1);
//  Probs.push_back(0.2);
//  Probs.push_back(0.3);
//  Probs.push_back(0.4);
  if (getConditionalHapPairProbsVals->size() < 1) {
    throw CppUnit::Exception(CppUnit::Message("MockCompositeLocus::getConditionalHapPairProbs(): No values to return"),
      CPPUNIT_SOURCELINE());
  }
  // return getConditionalHapPairProbsVals.pop_front();
  // TODO: implement getting those values.
  vector<double> v = getConditionalHapPairProbsVals->front();
  Probs.assign(v.begin(), v.end());
  getConditionalHapPairProbsVals->pop_front();
  return;
}

//void MockCompositeLocus::getConditionalHapPairProbs(std::vector<double>& Probs,
//      const std::vector<hapPair> &PossibleHapPairs,
//      const int ancestry[2])const
//{
//  getConditionalHapPairProbs_mocker.forward(Probs, PossibleHapPairs, ancestry);
//}

void MockCompositeLocus::getConditionalHapPairProbsAddReturnValue(vector<double> vd)
{
  getConditionalHapPairProbsVals->push_back(vd);
}
