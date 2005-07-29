// *-*-C++-*-*
#ifndef STRATIFICATIONTEST
#define STRATIFICATIONTEST 1

#include <vector>
#include <iostream>
#include "IndividualCollection.h"
#include "Genome.h"
#include "AlleleFreqs.h"
#include "LogWriter.h"

class StratificationTest
{
private:
  int T;
  int count;
  int NumberOfTestLoci;
  std::vector<unsigned int> TestLoci;
  bool ModelIndicator;
  std::ofstream outputstream;

  void OpenOutputFile( const char * , LogWriter *);

public:
   StratificationTest();

  void Initialize( AdmixOptions*, Genome& , LogWriter *);

  void calculate( IndividualCollection* individuals, double **AlleleFreqs, std::vector<std::vector<int> > ChrmAndLocus, int Populations );

  void Output();

  std::vector<double>
  GenerateExpectedGenotype( Individual*, double*, int  );

  std::vector<unsigned short>
  GenerateRepGenotype( double*, int ancestry[2] );

  unsigned short **SampleForOrderedSNiP ( double*, int ancestry[2] );

  //float getStatistic();
};

#endif /* !defined STRATIFICATIONTEST_H */
