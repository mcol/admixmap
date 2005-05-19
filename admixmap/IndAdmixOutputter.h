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
private: // members
  std::ofstream _out;

  AdmixOptions* _options;
  Genome* _Loci;
  std::string* _PopulationLabels;

  int _iterations;
  int _totalIndividuals;
  int _currentIndividual;

  int _analysisType;
  bool _isLogistic;
   bool _ModelIndicator;

  // UNIMPLEMENTED
  // to avoid use
  IndAdmixOutputter();
  IndAdmixOutputter(const IndAdmixOutputter&);
  IndAdmixOutputter& operator=(const IndAdmixOutputter&);

public:
  IndAdmixOutputter(AdmixOptions*,Genome*,std::string*);
  virtual ~IndAdmixOutputter();
  void visitIndividual(Individual&, std::vector<int>,double);
  void visitIndividualCollection(IndividualCollection&);

};

#endif /* !defined IND_ADMIX_OUTPUTTER */
