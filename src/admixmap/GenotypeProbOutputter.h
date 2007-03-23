// *-*-C++-*-*
#ifndef GPOUTPUTTER_H
#define GPOUTPUTTER_H

#include <fstream>
#include <iomanip>
#include "CompositeLocus.h"
#include "utils/pvector.h"

class GenotypeProbOutputter{
public:

  void Initialise(unsigned Nindivs, unsigned Nloci);
  //TODO: alternative for parallel version since HapPairProbs are not stored and Workers have no CompositeLocus objects
  //update for a single individual at a single locus
  void Update(unsigned i, unsigned j, const vector<vector<double> >&);
  void Output(const char* filename);

private:
  unsigned NumMaskedIndivs;
  unsigned NumMaskedLoci;
  unsigned NumIterations;
  pvector<double> Probs;
  std::ofstream outfile;

  /** Three-dimensional vector for probabilities.
   * Dimensions are:
   * 1. Individual
   * 2. Locus
   * 3. Genotype
   */
  vector<vector<vector<double> > > sumProbs;
};



#endif 
