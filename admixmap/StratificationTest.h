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
  std::ofstream DICstream;

  void Open( const char * , LogWriter *);

public:
   StratificationTest();

  void Initialize( AdmixOptions*, Genome& , LogWriter *);

  void calculate( IndividualCollection* individuals, Matrix_d *AlleleFreqs, std::vector<std::vector<int> > ChrmAndLocus );

  void Output();

  std::vector<double>
  GenerateExpectedGenotype( Individual*, const Matrix_d& );

  std::vector<unsigned short>
  GenerateRepGenotype( const Matrix_d&, const Vector_i& );

  unsigned short **SampleForOrderedSNiP ( const Matrix_d&, const Vector_i& );

  float getStatistic();
};

#endif /* !defined STRATIFICATIONTEST_H */
