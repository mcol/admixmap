// *-*-C++-*-*
#ifndef ALLELE_FREQ_OUTPUTTER_H
#define ALLELE_FREQ_OUTPUTTER_H 1

#include "LocusVisitor.h"
#include "CompositeLocus.h"
#include "Chromosome.h"
#include "Genome.h"

#include <iostream>
#include <string>


class AlleleFreqOutputter : public LocusVisitor
{
private:
  AdmixOptions* _options; 
  std::ofstream _out; 
  std::string* _PopulationLabels;
  int _iterations;
  int _numLoci;

  // UNIMPLEMENTED
  // to avoid use
  AlleleFreqOutputter();
  AlleleFreqOutputter(const AlleleFreqOutputter&);
  AlleleFreqOutputter& operator=(const AlleleFreqOutputter&);


public:
  AlleleFreqOutputter( AdmixOptions* options, std::string* PopulationLabels );
  virtual ~AlleleFreqOutputter( );

  virtual void visitCompositeLocus(CompositeLocus&);
  virtual void visitChromosome(Chromosome&);
  virtual void visitGenome(Genome&);

};

#endif /* !defined ALLELE_FREQ_OUTPUTTER_H */
