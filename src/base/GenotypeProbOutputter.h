// *-*-C++-*-*
#ifndef GPOUTPUTTER_H
#define GPOUTPUTTER_H

#include <fstream>
#include <iomanip>
#include "CompositeLocus.h"
#include "bcppcl/pvector.h"
#include "bcppcl/RObjectWriter.h"

class GenotypeProbOutputter{
public:

  void Initialise(unsigned Nindivs, unsigned Nloci);
  void Update(unsigned i, unsigned j, const std::vector<std::vector<double> >&);
  void Output(const char* filename, const std::vector<std::string>& LocusLabels);

private:
  unsigned NumMaskedIndivs;
  unsigned NumMaskedLoci;
  unsigned NumIterations;
  bcppcl::pvector<double> Probs;
  //std::ofstream outfile;
  RObjectWriter outfile;

  /** Three-dimensional vector for probabilities.
   * Dimensions are:
   * 1. Individual
   * 2. Locus
   * 3. Genotype
   */
  vector<vector<vector<double> > > sumProbs;
};



#endif 
