#ifndef HAPMIXGENOME_H
#define HAPMIXGENOME_H

#include "Genome.h"

/**
   classes to enable use of polymorphic HMM
   HMM objects are members of Chromosome class and Chromosome objects are members of Genome class.
   So all 3 must be polymorphic.
*/

class HapMixChromosome : public Chromosome{
public:
  HapMixChromosome(int n, int size, int start, int inNumHiddenStates, bool isx);

  void SetLocusCorrelation(const std::vector<double>::const_iterator lambda);
private:
  // UNIMPLEMENTED
  // to avoid use
  HapMixChromosome(const HapMixChromosome&);
  HapMixChromosome& operator=(const HapMixChromosome&);
  HapMixChromosome();


};

class HapMixGenome : public Genome{
public:
  void SetLocusCorrelation(const std::vector<double>& lambda);

  HapMixChromosome* getHapMixChromosome(unsigned c);
private:
  std::vector<HapMixChromosome*> vHapMixChr;

  void CreateChromosome(unsigned i, unsigned size, bool isX, unsigned cstart, int NumHiddenStates );


};

#endif
