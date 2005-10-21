// *-*-C++-*-*
#ifndef IND_ADMIX_OUTPUTTER
#define IND_ADMIX_OUTPUTTER 1

#include "Individual.h"
#include "IndividualCollection.h"
#include "AdmixOptions.h"
#include "Genome.h"
#include "Latent.h"

#include <vector>
#include <iostream>

class IndividualCollection;
class Individual;

class IndAdmixOutputter 
{
public:
  IndAdmixOutputter(const AdmixOptions* const, const Genome* const, const std::string* const);
  virtual ~IndAdmixOutputter();
  void visitIndividual(const Individual&, const std::vector<int>, double);
  void visitIndividualCollection(const IndividualCollection&);


private: 
  std::ofstream _out;

  const AdmixOptions* _options;
  const Genome* _Loci;
  const std::string* _PopulationLabels;

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
