#ifndef MOCKINPUTDATA_H_
#define MOCKINPUTDATA_H_



#define MOCKPP_IMPORT_ABBREVIATED
#include <mockpp/mockpp.h>
#include <mockpp/visiting/VisitableMockObject.h>
#include <mockpp/chaining/ChainingMockObjectSupport.h>
#include "../../admixmap/interfaces/IGenome.h"
#include "../../admixmap/interfaces/IInputData.h"

#include "convenience-macros.h"

#ifdef HAVE_CPPUNIT
#include <cppunit/Exception.h>
#endif

// class IGenome;

USING_NAMESPACE_MOCKPP

/** MockInputData, a class to mimick InputData class behaviour.
 * 
 * Implements the same interface, IInputData.
 */
class MockInputData : public VisitableMockObject,
                      public IInputData
{
private:
  std::vector<genotype> *GetGenotypeDefaultReturnValue;
  bool *GetHapMixGenotypeDefaultReturnValue1;
  std::vector<unsigned short> *GetHapMixGenotypeDefaultReturnValue2;
public:
	MockInputData()
  : VisitableMockObject(MOCKPP_PCHAR("MockInputData"), 0)
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(GetGenotype, ext1)
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT5(GetGenotype, ext5)
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(isFemale, ext)
  {
    GetGenotypeDefaultReturnValue = NULL;
    GetHapMixGenotypeDefaultReturnValue1 = NULL;
    GetHapMixGenotypeDefaultReturnValue2 = NULL;
  }
	virtual ~MockInputData();
//  vector<unsigned short> GetGenotype(
//      const string genostring) const;
//  bool isFemale(int i) const;
  MOCKPP_CONST_VISITABLE_EXT1(MockInputData, vector<unsigned short>, GetGenotype, const string,
      vector<unsigned short>, ext1, string);
  void GetGenotype(
      int, int, const IGenome&,
      std::vector<genotype>*,
      bool **) const;
//  MOCKPP_VOID_CONST_VISITABLE_EXT5(MockInputData, GetGenotype,
//      int, int, const IGenome&, std::vector<genotype>*, bool **,
//      ext5, int, int, IGenome, std::vector<genotype>*, bool**)
  bool GetHapMixGenotype(
      int i,
      int SexColumn,
      const IGenome &Loci,
      std::vector<unsigned short>* genotypes,
      bool** Missing) const;
  MOCKPP_CONST_VISITABLE_EXT1(MockInputData, bool, isFemale, int,
      bool, ext, int);
  
  void setGetGenotypeDefaultReturnValue(std::vector<genotype>);
  void setGetHapMixGenotypeDefaultReturnValues(bool, std::vector<unsigned short>);
};

#endif /*MOCKINPUTDATA_H_*/
