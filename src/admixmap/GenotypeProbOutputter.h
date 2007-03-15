// *-*-C++-*-*
#ifndef GPOUTPUTTER_H
#define GPOUTPUTTER_H

#include <fstream>
#include "CompositeLocus.h"
#include "../common/pvector.h"

class GenotypeProbOutputter{
public:

  void Initialise(unsigned Nindivs, unsigned Nloci);
  //TODO: alternative for parallel version since HapPairProbs are not stored and Workers have no CompositeLocus objects
  //update for a single individual at a single locus
  void Update(unsigned i, unsigned j, const CompositeLocus* Locus, const std::vector<hapPair > &HapPairs, const int ancestry[2]);
  void Output(const char* filename);

private:
  unsigned NumMaskedIndivs;
  unsigned NumMaskedLoci;
  unsigned NumIterations;
  pvector<double> Probs;
  std::ofstream outfile;

  std::vector< double > SumProbs;
};



#endif 
