#ifndef MOCKCOMPOSITELOCUS_H_
#define MOCKCOMPOSITELOCUS_H_

//#define MOCKPP_ENABLE_DEFAULT_FORMATTER 
//#include <mockpp/compat/Formatter.h> 

#include <list>
using std::list;

#define MOCKPP_IMPORT_ABBREVIATED
#include <mockpp/mockpp.h>
#include <mockpp/visiting/VisitableMockObject.h>
//#include <mockpp/visiting/CountedVisitableMethod.h>
#include <mockpp/chaining/ChainingMockObjectSupport.h>
#include "../../admixmap/interfaces/ICompositeLocus.h"

#include <cppunit/Exception.h>

USING_NAMESPACE_MOCKPP

/**
 * MockGenome, a class to mimick Genome class behaviour.
 * Implements the same interface, IGenome.
 */

class MockCompositeLocus :  public VisitableMockObject,
                            public ICompositeLocus
{
private:
  vector<double> vecDouble;
  list<vector<double> > *getConditionalHapPairProbsVals;
//  VisitableMockMethod<void, std::vector<double>&, const std::vector<hapPair> &, const int *>   getConditionalHapPairProbs_mocker;
public:
	MockCompositeLocus()
  : VisitableMockObject(MOCKPP_PCHAR("MockGenome"), 0)
  // void AccumulateAlleleProbs() = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE0(AccumulateAlleleProbs)
  // void AddLocus( int, std::string ) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT2(AddLocus, ext)
  // const std::string GetLabel(int)const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(GetLabel, ext)
  // int GetNumberOfAllelesOfLocus( int )const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(GetNumberOfAllelesOfLocus, ext)
  // int GetNumberOfStates()const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(GetNumberOfStates)
  // int GetNumberOfLoci()const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(GetNumberOfLoci)
  // void InitialiseHapPairProbs(const double* const allelefreqs) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT1(InitialiseHapPairProbs, ext)
  // void setAlleleProbsMAP(const double* const Freqs) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT1(setAlleleProbsMAP, ext)
  // void InitialiseHapPairProbsMAP() = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE0(InitialiseHapPairProbsMAP)
  // void SetHapPairProbs() = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE0(SetHapPairProbs)
  // const std::vector<int> getAlleleCounts(int a, const int* happair)const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT2(getAlleleCounts, ext)
  // void getConditionalHapPairProbs(std::vector<double>& Probs, const std::vector<hapPair> &PossibleHapPairs, const int ancestry[2])const = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT3(getConditionalHapPairProbs, ext)
  // void SampleHapPair(hapPair*, const std::vector<hapPair > &PossibleHapPairs, const int ancestry[2])const = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT3(SampleHapPair, ext)
  // const int *GetHapLabels( int ) const = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(GetHapLabels, ext)
  // const vector<int> getHaplotypeCounts(const int* happair) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(getHaplotypeCounts, ext)
  // void getLocusAlleleProbs(double **P, int k)const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT2(getLocusAlleleProbs, ext)
  // void SetDefaultMergeHaplotypes( const double* const alpha) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT1(SetDefaultMergeHaplotypes, ext)
  // int GetNumberOfMergedHaplotypes()const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT0(GetNumberOfMergedHaplotypes, ext)
  // void SetHapPairProbsToPosteriorMeans(int iterations) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT1(SetHapPairProbsToPosteriorMeans, ext)
  // void GetGenotypeProbs(double *Probs, const std::vector<hapPair> &HaplotypePairs, bool chibindicator)const = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT1(AddLocus, ext)
  // void GetHaploidGenotypeProbs(double *Probs, const std::vector<hapPair > &HapPairs, bool chibindicator) const = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(AddLocus, ext)
  {
    getConditionalHapPairProbsVals = new list<vector<double> >;
  }
	virtual ~MockCompositeLocus();
  // void AccumulateAlleleProbs() = 0;
  MOCKPP_VOID_VISITABLE0(MockCompositeLocus, AccumulateAlleleProbs);
  // void AddLocus( int, std::string ) = 0;
  MOCKPP_VOID_VISITABLE_EXT2(MockCompositeLocus, AddLocus, int, std::string,
      ext, int, std::string);
  // const std::string GetLabel(int)const = 0;
  MOCKPP_CONST_VISITABLE_EXT1(MockCompositeLocus, const std::string, GetLabel, int,
      std::string, ext, int);
  // int GetNumberOfAllelesOfLocus( int )const = 0;
  MOCKPP_CONST_VISITABLE_EXT1(MockCompositeLocus, int, GetNumberOfAllelesOfLocus, int,
        int, ext, int);
  // int GetNumberOfStates()const = 0;
  MOCKPP_CONST_VISITABLE0(MockCompositeLocus, int, GetNumberOfStates);
  // int GetNumberOfLoci()const = 0;
  MOCKPP_CONST_VISITABLE0(MockCompositeLocus, int, GetNumberOfLoci);
  // void InitialiseHapPairProbs(const double* const allelefreqs) = 0;
  MOCKPP_VOID_VISITABLE_EXT1(MockCompositeLocus, InitialiseHapPairProbs,
      const double * const, ext, double *);
  // void setAlleleProbsMAP(const double* const Freqs) = 0;
  MOCKPP_VOID_VISITABLE_EXT1(MockCompositeLocus, setAlleleProbsMAP,
      const double * const, ext, double *);
  // void InitialiseHapPairProbsMAP() = 0;
  MOCKPP_VOID_VISITABLE0(MockCompositeLocus, InitialiseHapPairProbsMAP);
  // void SetHapPairProbs() = 0;
  MOCKPP_VOID_VISITABLE0(MockCompositeLocus, SetHapPairProbs);
  // const std::vector<int> getAlleleCounts(int a, const int* happair)const = 0;
  MOCKPP_CONST_VISITABLE_EXT2(MockCompositeLocus, const std::vector<int>, getAlleleCounts, int, const int*,
      std::vector<int>, ext, int, int*);
  // void getConditionalHapPairProbs(std::vector<double>& Probs, const std::vector<hapPair> &PossibleHapPairs, const int ancestry[2])const = 0;
//  MOCKPP_VOID_CONST_VISITABLE_EXT3(MockCompositeLocus, getConditionalHapPairProbs, std::vector<double>&, const std::vector<hapPair> &, const int *,
//      ext, std::vector<double>, std::vector<hapPair> , int *);
  // void SampleHapPair(hapPair*, const std::vector<hapPair> &PossibleHapPairs, const int ancestry[2])const = 0;
//  MOCKPP_VOID_CONST_VISITABLE_EXT3(MockCompositeLocus, SampleHapPair, hapPair*, const std::vector<hapPair> &, const int *,
//      ext, hapPair *, std::vector<hapPair> , int *);
  // const int *GetHapLabels( int ) const = 0;
//  MOCKPP_CONST_VISITABLE_EXT1(MockCompositeLocus, const int*, GetHapLabels, int,
//      int*, ext, int);
  // const vector<int> getHaplotypeCounts(const int* happair) = 0;
  MOCKPP_VISITABLE_EXT1(MockCompositeLocus, const vector<int>, getHaplotypeCounts, const int*,
      vector<int>, ext, int*);
  // void getLocusAlleleProbs(double **P, int k)const = 0;
  MOCKPP_VOID_CONST_VISITABLE_EXT2(MockCompositeLocus, getLocusAlleleProbs,
      double **, int,
      ext, double **, int);
  // void SetDefaultMergeHaplotypes( const double* const alpha) = 0;
  MOCKPP_VOID_VISITABLE_EXT1(MockCompositeLocus, SetDefaultMergeHaplotypes,
      const double *,
      ext, double *);
  // int GetNumberOfMergedHaplotypes()const = 0;
  MOCKPP_CONST_VISITABLE_EXT0(MockCompositeLocus, int, GetNumberOfMergedHaplotypes,
      int, ext);
  // void SetHapPairProbsToPosteriorMeans(int iterations) = 0;
  MOCKPP_VOID_VISITABLE_EXT1(MockCompositeLocus, SetHapPairProbsToPosteriorMeans,
      int, ext, int);
  // void GetGenotypeProbs(double *Probs, const std::vector<hapPair> &HaplotypePairs, bool chibindicator)const = 0;
  // void GetHaploidGenotypeProbs(double *Probs, const std::vector<hapPair > &HapPairs, bool chibindicator) const = 0;

  // Some methods couldn't be mocked with mockpp; mostly the ones which
  // * Return a pointer
  // * Take a vector<SomeClass> as an argument
  
  /// Mockpp couldn't help with this function. It was implemented by hand.
  virtual void getConditionalHapPairProbs(std::vector<double>& Probs,
      const std::vector<hapPair> &PossibleHapPairs,
      const int ancestry[2])const; /*
  {
    throw CppUnit::Exception(
      CppUnit::Message("MockCompositeLocus::getConditionalHapPairProbs(std::vector<double>& Probs, const std::vector<hapPair> &PossibleHapPairs, const int ancestry[2])const is not implemented."),
      CPPUNIT_SOURCELINE());
  } */
  
  // Following functions don't mock anything, they only throw exceptions.
  void SampleHapPair(hapPair*,
      const std::vector<hapPair> &PossibleHapPairs,
      const int ancestry[2])const
  {
    throw CppUnit::Exception(
      CppUnit::Message("MockCompositeLocus::SampleHapPair(hapPair*, const std::vector<hapPair> &PossibleHapPairs, const int ancestry[2])const is not implemented."),
      CPPUNIT_SOURCELINE());
  }
  const int *GetHapLabels( int ) const
  {
    throw CppUnit::Exception(
      CppUnit::Message("const int *MockCompositeLocus::GetHapLabels( int ) const is not implemented"),
      CPPUNIT_SOURCELINE());
  }
  void GetGenotypeProbs(double *Probs,
      const std::vector<hapPair> &HaplotypePairs,
      bool chibindicator)const
  {
    throw CppUnit::Exception(
      CppUnit::Message("MockCompositeLocus::GetGenotypeProbs(double *Probs, const std::vector<hapPair> &HaplotypePairs, bool chibindicator)const is not implemented"),
      CPPUNIT_SOURCELINE());
  }
  void GetHaploidGenotypeProbs(double *Probs,
      const std::vector<hapPair > &HapPairs, 
      bool chibindicator) const
  {
    throw CppUnit::Exception(
      CppUnit::Message("void MockCompositeLocus::GetHaploidGenotypeProbs(double *Probs, const std::vector<hapPair > &HapPairs,  bool chibindicator) const is not implemented"),
      CPPUNIT_SOURCELINE());
  }
  
  void getConditionalHapPairProbsAddReturnValue(vector<double>);
};

//String & operator<< (String &s, const std::vector<hapPair> &pers)
//{
//  s << "std::vector<hapPair>";
//  return s;
//}
//
//String & operator << (String &s, const std::vector<double> &pers)
//{
//  s << "std::vector<double>";
//  return s;
//}

//bool std::vector<hapPair>::operator == (const std::vector<hapPair> &pers) const
//{
//  return false;
//}

#endif /*MOCKCOMPOSITELOCUS_H_*/
