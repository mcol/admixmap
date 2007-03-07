#ifndef IINPUTDATA_H_
#define IINPUTDATA_H_

#include <vector>
#include <string>
#include "../GeneticDistanceUnit.h"
#include "../interfaces/IGenome.h"

using std::vector;
using std::string;

typedef vector<vector<unsigned short> > genotype;

// class IGenome;

/** Abstract class defining interface for the InputData class.
 * 
 * It has two implementations:
 * InputData
 * MockInputData for testing
 */

class IInputData
{
public:
  // IInputData();
  virtual ~IInputData() = 0;
  virtual vector<unsigned short> GetGenotype(
      const string genostring) const = 0;
  virtual bool isFemale(int i) const = 0;
  
  /** Get genotypes from the IGenome and insert them into genotype vector */
  virtual void GetGenotype(
      int, int, const IGenome&,
      std::vector<genotype>*,
      bool **) const = 0;
  virtual bool GetHapMixGenotype(
      int i,
      int SexColumn,
      const IGenome &Loci,
      std::vector<unsigned short>* genotypes,
      bool** Missing) const = 0;
};

#endif /*IINPUTDATA_H_*/
