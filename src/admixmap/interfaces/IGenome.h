#ifndef IGENOME_H_
#define IGENOME_H_

#include <vector>
//#include "../interfaces/ICompositeLocus.h"
#include "../interfaces/IChromosome.h"
// class ICompositeLocus;
//class IChromosome;
class InputData;
class LogWriter;
class Distances;
class ICompositeLocus;

class IGenome
{
public:
	// IGenome();
	virtual ~IGenome();
  
  virtual void GetChrmAndLocus(unsigned locus, unsigned* c, unsigned* l) = 0;
  virtual unsigned GetChrNumOfLocus(unsigned locus) = 0;
//  virtual IChromosome* getChromosome(unsigned) = 0;
  // Trying to use a reference (to Chromosome) instead of a pointer.
  virtual IChromosome *getChromosome(unsigned) = 0;
  virtual const int getChromosomeNumber(int) const = 0;
  virtual unsigned getFirstXLocus()const = 0;
  virtual double GetLengthOfGenome()const = 0;
  virtual double GetLengthOfXchrm()const = 0;
  virtual unsigned int GetNumberOfChromosomes()const = 0;
  virtual unsigned int GetNumberOfCompositeLoci()const = 0;
  virtual int getNumberOfLoci(int)const = 0;
  virtual int GetNumberOfStates()const = 0;
  virtual const int getRelativeLocusNumber(int) const = 0;
  virtual unsigned GetSizeOfChromosome(unsigned)const = 0;
  virtual const unsigned int *GetSizesOfChromosomes()const = 0;
  virtual unsigned int GetTotalNumberOfLoci()const = 0;
  virtual bool isX_data()const = 0;
  virtual unsigned isXChromosome(unsigned)const = 0;
  virtual bool isXLocus(unsigned j)const = 0;
  virtual int GetNumberOfStates(int) const = 0;
  virtual ICompositeLocus* operator()(int) const = 0;
  virtual const IChromosome* const* getChromosomes()const = 0;
  virtual void Initialise(const InputData* const data_, int populations, LogWriter &Log) = 0;
  virtual void PrintLocusTable(const char* filename, const std::vector<double>& Distances)const = 0;
  virtual void SetLocusCorrelation(double rho) = 0;
  virtual void SetLocusCorrelation(const std::vector<double> rho) = 0;
  virtual double GetDistance(int)const = 0;
//  virtual  = 0;
//  virtual  = 0;
};

#endif /*IGENOME_H_*/
