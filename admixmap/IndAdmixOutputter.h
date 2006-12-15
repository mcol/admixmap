// *-*-C++-*-*
#ifndef IND_ADMIX_OUTPUTTER
#define IND_ADMIX_OUTPUTTER 1

#include "IndividualCollection.h"
#include "AdmixOptions.h"
#include "Genome.h"

#include <vector>
#include <iostream>

class IndividualCollection;
class Individual;

///Class to output individual admixture proportions and sumintensities to file
class IndAdmixOutputter 
{
public:
  IndAdmixOutputter(const AdmixOptions* const, const Genome* const, const Vector_s& PopLabels);
  virtual ~IndAdmixOutputter();
  void visitIndividual(const AdmixedIndividual&, const std::vector<int>);
  void visitIndividualCollection(const IndividualCollection&);


private: 
  std::ofstream _out;

  const AdmixOptions* _options;
  const Genome* _Loci;
  const Vector_s*  _PopulationLabels;

  int _iterations;
  int _totalIndividuals;
  int _currentIndividual;

  bool _RandomMatingModelIndicator;

  // UNIMPLEMENTED
  // to avoid use
  IndAdmixOutputter();
  IndAdmixOutputter(const IndAdmixOutputter&);
  IndAdmixOutputter& operator=(const IndAdmixOutputter&);

};

#endif /* !defined IND_ADMIX_OUTPUTTER */
