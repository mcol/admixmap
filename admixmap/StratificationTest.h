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

  void Initialize( AdmixOptions* const options, const Genome &Loci,  
				     const IndividualCollection* const IC, LogWriter &Log );

  void calculate( const IndividualCollection* const individuals, const double* const* AlleleFreqs,
		  const std::vector<std::vector<int> > ChrmAndLocus, int Populations );

  void Output(LogWriter &);

  //float getStatistic();

private:
  int T;
  int count;
  int NumberOfTestLoci;
  std::vector<unsigned int> TestLoci;
  bool ModelIndicator;
  std::ofstream outputstream;

  void OpenOutputFile( const char * , LogWriter &);

  std::vector<double>
  GenerateExpectedGenotype( const Individual* const, const double*, const int  );

  std::vector<unsigned short>
  SimGenotypeConditionalOnAdmixture( const std::vector<double> ProbAllele1 );

  std::vector<unsigned short>
  SimGenotypeConditionalOnAncestry( const double*, const int ancestry[2] );

  std::vector<unsigned short> SampleHeterozygotePhase( const double*, const int ancestry[2] );

  int GetAlleleCounts(int locus, int a, const IndividualCollection* const IC);

};

#endif /* !defined STRATIFICATIONTEST_H */
