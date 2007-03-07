#ifndef MOCKGENOME_H_
#define MOCKGENOME_H_

#define MOCKPP_IMPORT_ABBREVIATED
#include <mockpp/mockpp.h>
#include <mockpp/visiting/VisitableMockObject.h>
#include <mockpp/chaining/ChainingMockObjectSupport.h>
#include <cppunit/Exception.h>

#include "../../admixmap/interfaces/IGenome.h"
#include "../../admixmap/interfaces/IChromosome.h"
#include "utils/LogWriter.h"
#include "MockCompositeLocus.h"
#include "MockChromosome.h"

#include "convenience-macros.h"

USING_NAMESPACE_MOCKPP

/** MockGenome, a class to mimick Genome class behaviour.
 * Implements the same interface, IGenome.
 */

class MockGenome :  public VisitableMockObject,
                    public IGenome
{
private:
  MockChromosome *chromosome;
  MockCompositeLocus *cl;
public:
	MockGenome()
  : VisitableMockObject(MOCKPP_PCHAR("MockGenome"), 0)
  // void GetChrmAndLocus(unsigned locus, unsigned* c, unsigned* l) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT3(GetChrmAndLocus, ext3)
  // unsigned GetChrNumOfLocus(unsigned locus) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(GetChrNumOfLocus, ext1)
  //  IChromosome& getChromosome(unsigned) = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(getChromosome, ext)
  //  const int getChromosomeNumber(int) const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(getChromosomeNumber, constint)
  //  unsigned getFirstXLocus()const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(getFirstXLocus)
  //  double GetLengthOfGenome()const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(GetLengthOfGenome)
  //  double GetLengthOfXchrm()const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(GetLengthOfXchrm)
  //  unsigned int GetNumberOfChromosomes()const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(GetNumberOfChromosomes)
  //  unsigned int GetNumberOfCompositeLoci()const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(GetNumberOfCompositeLoci)
  //  int getNumberOfLoci(int)const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(getNumberOfLoci, ext)
  //  int GetNumberOfStates()const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT0(GetNumberOfStates, ext0)
  //  const int getRelativeLocusNumber(int) const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(getRelativeLocusNumber, ext)
  //  unsigned GetSizeOfChromosome(unsigned)const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(GetSizeOfChromosome, ext)
  //  const unsigned int *GetSizesOfChromosomes()const = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT0(GetSizesOfChromosomes, ext)
  //  unsigned int GetTotalNumberOfLoci()const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(GetTotalNumberOfLoci)
  //  bool isX_data()const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE0(isX_data)
  //  unsigned isXChromosome(unsigned)const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(isXChromosome, ext)
  //  bool isXLocus(unsigned j)const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(isXLocus, ext)
  //  int GetNumberOfStates(int) const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(GetNumberOfStates, ext1)
  //  CompositeLocus* operator()(int) const = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE1(GetChromosome)
  //  const Chromosome* const* getChromosomes()const = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT0(getChromosomes, ext)
  //  void Initialise(const InputData* const data_, int populations, LogWriter &Log) = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT3(Initialise, ext)
  //  void PrintLocusTable(const char* filename, const std::vector<double>& Distances)const = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(PrintLocusTable, ext)
  //  void SetLocusCorrelation(double rho) = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT1(SetLocusCorrelation, dbl)
  //  void SetLocusCorrelation(const std::vector<double> rho) = 0;
//  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VOID_VISITABLE_EXT1(SetLocusCorrelation, vecDouble)
  //  double GetDistance(int)const = 0;
  , MOCKPP_CONSTRUCT_MEMBERS_FOR_VISITABLE_EXT1(GetDistance, ext)
//  , chromosome(NULL)
  {
    chromosome = NULL;
    cl = NULL;
  }
	~MockGenome();
  
  void setMockChromosome(MockChromosome*);
  void setMockCompositeLocus(MockCompositeLocus*);
  
  // void GetChrmAndLocus(unsigned locus, unsigned* c, unsigned* l)
  MOCKPP_VOID_VISITABLE_EXT3(MockGenome, GetChrmAndLocus,
            unsigned, unsigned *, unsigned *,
            ext3, unsigned, unsigned *, unsigned *);
  
  // unsigned GetChrNumOfLocus(unsigned locus) = 0;
  MOCKPP_VISITABLE_EXT1(MockGenome, unsigned, GetChrNumOfLocus, unsigned,
      unsigned, ext1, unsigned);
  // IChromosome& getChromosome(unsigned) = 0;
//  MOCKPP_VISITABLE_EXT1(MockGenome, IChromosome&, getChromosome, unsigned,
//      MockChromosome, ext, unsigned);
  // const int getChromosomeNumber(int) const = 0;
  MOCKPP_CONST_VISITABLE_EXT1(MockGenome, const int, getChromosomeNumber, int,
      int, constint, int);
  // unsigned getFirstXLocus()const = 0;
  MOCKPP_CONST_VISITABLE0(MockGenome, unsigned, getFirstXLocus);
  // double GetLengthOfGenome()const = 0;
  MOCKPP_CONST_VISITABLE0(MockGenome, double, GetLengthOfGenome);
  // double GetLengthOfXchrm()const = 0;
  MOCKPP_CONST_VISITABLE0(MockGenome, double, GetLengthOfXchrm);
  // unsigned int GetNumberOfChromosomes()const = 0;
  MOCKPP_CONST_VISITABLE0(MockGenome, unsigned, GetNumberOfChromosomes);
  // unsigned int GetNumberOfCompositeLoci()const = 0;
  MOCKPP_CONST_VISITABLE0(MockGenome, unsigned, GetNumberOfCompositeLoci);
  // int getNumberOfLoci(int)const = 0;
  MOCKPP_CONST_VISITABLE_EXT1(MockGenome, int, getNumberOfLoci, int,
      int, ext, int);
  // int GetNumberOfStates()const = 0;
  MOCKPP_CONST_VISITABLE_EXT0(MockGenome, int, GetNumberOfStates,
    int, ext0);
  // const int getRelativeLocusNumber(int) const = 0;
  MOCKPP_CONST_VISITABLE_EXT1(MockGenome, const int, getRelativeLocusNumber, int,
        int, ext, int);
  // unsigned GetSizeOfChromosome(unsigned)const = 0;
  MOCKPP_CONST_VISITABLE_EXT1(MockGenome, unsigned, GetSizeOfChromosome, unsigned,
        unsigned, ext, unsigned);
  // const unsigned int *GetSizesOfChromosomes()const = 0;
//  MOCKPP_CONST_VISITABLE_EXT0(MockGenome, const unsigned *, GetSizesOfChromosomes,
//      unsigned, ext);
  // unsigned int GetTotalNumberOfLoci()const = 0;
  MOCKPP_CONST_VISITABLE0(MockGenome, unsigned, GetTotalNumberOfLoci);
  // bool isX_data()const = 0;
  MOCKPP_CONST_VISITABLE0(MockGenome, bool, isX_data);
  // unsigned isXChromosome(unsigned)const = 0;
  MOCKPP_CONST_VISITABLE_EXT1(MockGenome, unsigned, isXChromosome, unsigned,
      unsigned, ext, unsigned);
  // bool isXLocus(unsigned j)const = 0;
  MOCKPP_CONST_VISITABLE_EXT1(MockGenome, bool, isXLocus, unsigned,
      bool, ext, unsigned);
  // int GetNumberOfStates(int) const = 0;
  MOCKPP_CONST_VISITABLE_EXT1(MockGenome, int, GetNumberOfStates, int,
      int, ext1, int);
  // CompositeLocus* operator()(int) const = 0;
  // const Chromosome* const* getChromosomes()const = 0;
//  MOCKPP_CONST_VISITABLE_EXT0(MockGenome, const Chromosome * const *, getChromosomes,
//      Chromosome*, ext);
  // void Initialise(const InputData* const data_, int populations, LogWriter &Log) = 0;
//  MOCKPP_VOID_VISITABLE_EXT3(MockGenome, Initialise, const InputData*, int, LogWriter&,
//      ext, InputData *, int, LogWriter);
  // void PrintLocusTable(const char* filename, const std::vector<double>& Distances)const = 0;
  // void SetLocusCorrelation(double rho) = 0;
  MOCKPP_VOID_VISITABLE_EXT1(MockGenome, SetLocusCorrelation, double,
        dbl, double);
  // void SetLocusCorrelation(const std::vector<double> rho) = 0;
  // The following line fails. See:
  // http://sourceforge.net/forum/forum.php?thread_id=1685696&forum_id=235239
//  MOCKPP_VOID_VISITABLE_EXT1(MockGenome, SetLocusCorrelation, const std::vector<double>,
//      vecDouble, std::vector<double>);
  // double GetDistance(int)const = 0;
  MOCKPP_CONST_VISITABLE_EXT1(MockGenome, double, GetDistance, int,
      double, ext, int);
      

  /*
   * Following methods weren't easy to mock with Mockpp, so they are
   * declared here and implemented in the .cc file. Ideally, they
   * should be removed from here and mocked with Mockpp.
   */
  virtual IChromosome* getChromosome(unsigned int);/* {
    throw string("IChromosome& MockGenome::getChromosome(unsigned int) unimplemented");} */
  virtual const unsigned int* GetSizesOfChromosomes() const {
    throw string("const unsigned int* MockGenome::GetSizesOfChromosomes() unimplemented");}
  virtual ICompositeLocus* operator()(int) const;/* {
    throw string("CompositeLocus* MockGenome::operator()(int) const unimplemented");} */
  virtual const IChromosome* const* getChromosomes() const {
    throw string("const Chromosome* const* MockGenome::getChromosomes() const unimplemented");}
  virtual void Initialise(const InputData*, int, LogWriter&) {
    throw string("void MockGenome::Initialise(const InputData*, int, LogWriter&) unimplemented");}
  virtual void PrintLocusTable(const char*, const std::vector<double, std::allocator<double> >&) const {
    throw string("void MockGenome::PrintLocusTable(const char*, const std::vector<double, std::allocator<double> >&) const unimplemented");}
  virtual void SetLocusCorrelation(std::vector<double, std::allocator<double> >) {
    throw string("MockGenome::SetLocusCorrelation(std::vector<double, std::allocator<double> >) unimplemented");}
};

#endif /*MOCKGENOME_H_*/
