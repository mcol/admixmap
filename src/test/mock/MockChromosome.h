#ifndef MOCKCHROMOSOME_H_
#define MOCKCHROMOSOME_H_

#include <string>
using std::string;

#define MOCKPP_IMPORT_ABBREVIATED
#include <mockpp/mockpp.h>
#include <mockpp/visiting/VisitableMockObject.h>
#include <mockpp/chaining/ChainingMockObjectSupport.h>

USING_NAMESPACE_MOCKPP

#include <cppunit/Exception.h>
#include "../../admixmap/interfaces/IChromosome.h"
//#include "operators.h"

class GenotypeProbIterator;

class MockChromosome :  public VisitableMockObject,
                        public IChromosome
{
public:
	MockChromosome()
  : VisitableMockObject(MOCKPP_PCHAR("MockChromosome"), 0)
  // unsigned int GetSize()const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(GetSize)
  // bool isXChromosome()const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(isXChromosome)
  // void SetDistance(int,double) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT2(SetDistance, ext)
  // void SetLabel(std::string ) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT1(SetLabel, ext)
  // void SetLocusCorrelation(const double rho) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT1(SetLocusCorrelation, ext)
  // void SetLocusCorrelation(
  //    const vector<double> rho_,
  //    bool global,
  //    bool RandomMating) = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT3(SetLocusCorrelation, ext3)
  // const vector<double> getHiddenStateProbs(const bool, int) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT2(getHiddenStateProbs, ext)
  // double getLogLikelihood(const bool isDiploid) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(getLogLikelihood, ext)
  // void SetLocusCorrelation(const std::vector<double>::const_iterator rho) = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT1(SetLocusCorrelation, iter)
  // void SampleLocusAncestry(int *OrderedStates, bool diploid) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT2(SampleLocusAncestry, ext)
  // int GetLocus(int)const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(GetLocus, ext)
  // void SetGenotypeProbs(const GenotypeProbIterator& GenotypeProbs, const bool* const GenotypesMissing) = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT2(SetGenotypeProbs, ext)
  // void SetStateArrivalProbs(bool RandomMating, bool isdiploid) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT2(SetStateArrivalProbs, ext)
  // void SetHMMTheta(const double* const Admixture, bool RandomMating, bool diploid) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT3(SetHMMTheta, ext3)
  // std::vector<std::vector<double> > getAncestryProbs(const bool isDiploid, int) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT2(getAncestryProbs, ext)
  // void SampleJumpIndicators(const int* const LocusAncestry,
  //        const unsigned int gametes, 
  //        int *SumLocusAncestry,
  //        std::vector<unsigned> &SumN, 
  //        bool SampleArrivals)const = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT5(SampleJumpIndicators, ext5)
  {}
  ~MockChromosome();
  
  // unsigned int GetSize()const = 0;
  MOCKPP_CONST_VISITABLE0(MockChromosome, unsigned, GetSize);
  // bool isXChromosome()const = 0;
  MOCKPP_CONST_VISITABLE0(MockChromosome, bool, isXChromosome);
  // void SetDistance(int,double) = 0;
  MOCKPP_VOID_VISITABLE_EXT2(MockChromosome, SetDistance, int, double,
      ext, int, double);
  // void SetLabel(std::string ) = 0;
  MOCKPP_VOID_VISITABLE_EXT1(MockChromosome, SetLabel, std::string,
      ext, std::string);
  // void SetLocusCorrelation(const double rho) = 0;
  MOCKPP_VOID_VISITABLE_EXT1(MockChromosome, SetLocusCorrelation,
      const double, ext, double);
  // void SetLocusCorrelation(
  //    const vector<double> rho_,
  //    bool global,
  //    bool RandomMating) = 0;
//  MOCKPP_VOID_VISITABLE_EXT3(MockChromosome, SetLocusCorrelation,
//      const vector<double>, bool, bool,
//      ext3, vector<double>, bool, bool);
  // const vector<double> getHiddenStateProbs(const bool, int) = 0;
  MOCKPP_VISITABLE_EXT2(MockChromosome, const vector<double>, getHiddenStateProbs,
      const bool, int,
      vector<double>, ext, bool, int);
  // double getLogLikelihood(const bool isDiploid) = 0;
  MOCKPP_VISITABLE_EXT1(MockChromosome, double, getLogLikelihood, const bool,
      double, ext, bool);
  // void SetLocusCorrelation(const std::vector<double>::const_iterator rho) = 0;
//  MOCKPP_VOID_VISITABLE_EXT1(MockChromosome, SetLocusCorrelation,
//      const std::vector<double>::const_iterator,
//      iter, std::vector<double>::const_iterator);
  // void SampleLocusAncestry(int *OrderedStates, bool diploid) = 0;
  MOCKPP_VOID_VISITABLE_EXT2(MockChromosome, SampleLocusAncestry,
      int*, bool, ext, int*, bool);
  // int GetLocus(int)const = 0;
  MOCKPP_CONST_VISITABLE_EXT1(MockChromosome, int, GetLocus, int,
      int, ext, int);
  // void SetGenotypeProbs(const GenotypeProbIterator& GenotypeProbs,
  //    const bool* const GenotypesMissing) = 0;
//  MOCKPP_VOID_VISITABLE_EXT2(MockChromosome, SetGenotypeProbs,
//      const GenotypeProbIterator&, const bool*,
//      ext, GenotypeProbIterator, bool *);
  // void SetStateArrivalProbs(bool RandomMating, bool isdiploid) = 0;
  MOCKPP_VOID_VISITABLE_EXT2(MockChromosome, SetStateArrivalProbs,
      bool, bool, ext, bool, bool);
  // void SetHMMTheta(const double* const Admixture, bool RandomMating, bool diploid) = 0;
  MOCKPP_VOID_VISITABLE_EXT3(MockChromosome, SetHMMTheta,
      const double*, bool, bool, ext3, double *, bool, bool);
  // std::vector<std::vector<double> > getAncestryProbs(
  //    const bool isDiploid, int) = 0;
  MOCKPP_VISITABLE_EXT2(MockChromosome, vector<vector<double> >, getAncestryProbs,
      bool, int,
      vector<vector<double> >, ext, bool, int);
  // void SampleJumpIndicators(
  //        const int* const LocusAncestry,
  //        const unsigned int gametes, 
  //        int *SumLocusAncestry,
  //        std::vector<unsigned> &SumN, 
  //        bool SampleArrivals)const = 0;
//  MOCKPP_VOID_CONST_VISITABLE0(MockChromosome, SampleJumpIndicators,
//      const int*, const unsigned, int *, vector<unsigned>&, bool,
//      ext5, int*, unsigned, int *, vector<unsigned>, bool);

  // Not easy to mock with mockpp. Leaving unimplemented.
 virtual void SetLocusCorrelation(
        const vector<double> rho_,
        bool global,
        bool RandomMating) {
    throw CppUnit::Exception(CppUnit::Message("void MockChromosome::SetLocusCorrelation(const vector<double> rho_, bool global, bool RandomMating) is not implemented."),
        CPPUNIT_SOURCELINE()); }
        
  /*virtual const std::vector<double> getHiddenStateProbs(bool, int);  {
    throw CppUnit::Exception(CppUnit::Message("const std::vector<double> MockChromosome::getHiddenStateProbs(bool, int) is unimplemented"),
        CPPUNIT_SOURCELINE()); } */
        
  void SetLocusCorrelation(const std::vector<double>::const_iterator rho) {
    throw CppUnit::Exception(CppUnit::Message("MockChromosome::SetLocusCorrelation(const std::vector<double>::const_iterator rho) is unimplemented"),
        CPPUNIT_SOURCELINE()); }
        
  void SetGenotypeProbs(const GenotypeProbIterator&, const bool*) {
    throw CppUnit::Exception(CppUnit::Message("MockChromosome::SetGenotypeProbs(const GenotypeProbIterator&, const bool*) is unimplemented"),
        CPPUNIT_SOURCELINE()); }
        
  virtual void SampleJumpIndicators(const int*, unsigned int, int*, std::vector<unsigned int>&, bool) const {
    throw CppUnit::Exception(CppUnit::Message("MockChromosome::SampleJumpIndicators(const int*, unsigned int, int*, std::vector<unsigned int>&, bool) const is unimplemented"),
        CPPUNIT_SOURCELINE()); }
        
  MockChromosome(const IChromosome&);
  MockChromosome& operator=(const IChromosome&) {
    throw CppUnit::Exception(CppUnit::Message("MockChromosome& MockChromosome::operator = (MockChromosome&) is unimplemented"),
        CPPUNIT_SOURCELINE()); }
//    throw string("Thou shalt not use the assignment operator for MockChromosome.");
//  }
};




#endif /*MOCKCHROMOSOME_H_*/
