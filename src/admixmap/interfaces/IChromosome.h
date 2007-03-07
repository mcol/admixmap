#ifndef ICHROMOSOME_H_
#define ICHROMOSOME_H_

#include <string>
#include <vector>

using std::vector;

class GenotypeProbIterator;

class IChromosome
{
public:
//	IChromosome();
	virtual ~IChromosome();
  
  virtual unsigned int GetSize()const = 0;
  virtual bool isXChromosome()const = 0;
  virtual void SetDistance(int,double) = 0;
  virtual void SetLabel(std::string ) = 0;
  virtual void SetLocusCorrelation(const double rho) = 0;
  virtual void SetLocusCorrelation(
      const vector<double> rho_,
      bool global,
      bool RandomMating) = 0;
  virtual const vector<double> getHiddenStateProbs(const bool, int) = 0;
  virtual double getLogLikelihood(const bool isDiploid) = 0;
  virtual void SetLocusCorrelation(const std::vector<double>::const_iterator rho) = 0;
  virtual void SampleLocusAncestry(int *OrderedStates, bool diploid) = 0;
  virtual int GetLocus(int)const = 0;
  virtual void SetGenotypeProbs(const GenotypeProbIterator& GenotypeProbs, const bool* const GenotypesMissing) = 0;
  virtual void SetStateArrivalProbs(bool RandomMating, bool isdiploid) = 0;
  virtual void SetHMMTheta(const double* const Admixture, bool RandomMating, bool diploid) = 0;
  virtual std::vector<std::vector<double> > getAncestryProbs(const bool isDiploid, int) = 0;
  virtual void SampleJumpIndicators(const int* const LocusAncestry, const unsigned int gametes, 
          int *SumLocusAncestry, std::vector<unsigned> &SumN, 
          bool SampleArrivals)const = 0;
//  virtual  = 0;
//  virtual  = 0;
//  virtual  = 0;
//  virtual  = 0;
//  virtual  = 0;
//  virtual  = 0;
};

#endif /*ICHROMOSOME_H_*/
