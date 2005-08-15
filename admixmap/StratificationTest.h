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
public:
   StratificationTest();

  void Initialize( AdmixOptions* options, Genome &Loci, Chromosome **Chr, IndividualCollection *IC, LogWriter *Log );

  void calculate( IndividualCollection* individuals, double **AlleleFreqs, std::vector<std::vector<int> > ChrmAndLocus, int Populations );

  void Output(LogWriter *);

  //float getStatistic();

private:
  int T;
  int count;
  int NumberOfTestLoci;
  std::vector<unsigned int> TestLoci;
  bool ModelIndicator;
  std::ofstream outputstream;

  void OpenOutputFile( const char * , LogWriter *);

  std::vector<double>
  GenerateExpectedGenotype( Individual*, const double*, const int  );

  std::vector<unsigned short>
  SimGenotypeConditionalOnAdmixture( std::vector<double> ProbAllele1 );

  std::vector<unsigned short>
  SimGenotypeConditionalOnAncestry( const double*, const int ancestry[2] );

  std::vector<unsigned short> SampleHeterozygotePhase( const double*, const int ancestry[2] );

  int GetAlleleCounts(int locus, int a, IndividualCollection *IC);

};

#endif /* !defined STRATIFICATIONTEST_H */
