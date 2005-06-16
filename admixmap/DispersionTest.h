// *-*-C++-*-*
#ifndef DISPTEST_H
#define DISPTEST_H 1

#include "IndividualCollection.h"
#include "Genome.h"
#include "AlleleFreqs.h"
#include "LogWriter.h"
#include "matrix_i.h"
#include "matrix.h"

class DispersionTest{
 private:
  AdmixOptions *options;
  std::ofstream dispersionoutputstream;
  Matrix_i divergentallelefreqstest;

 public:
  DispersionTest();
  void Initialise(AdmixOptions *,LogWriter *, int);
  void Output(int , Genome &, std::string *PopLabels);

  void TestForDivergentAlleleFrequencies(AlleleFreqs *);
};











#endif /* !defined DISPTEST_H */
