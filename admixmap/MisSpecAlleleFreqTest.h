// *-*-C++-*-*
#ifndef MISSPECALLELEFREQTEST_H
#define MISSPECALLELEFREQTEST_H 1

#include "common.h"
#include <iomanip>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "IndividualCollection.h"
#include "Individual.h"
#include "AlleleFreqs.h"
#include "Genome.h"
#include "AdmixOptions.h"
#include "LogWriter.h"

class MisSpecAlleleFreqTest{

public:
  MisSpecAlleleFreqTest();
  ~MisSpecAlleleFreqTest();

  void Initialise(const AdmixOptions* const options, const Genome* const Loci, LogWriter *Log );
  void Update(const IndividualCollection* const individuals, const AlleleFreqs* const A, const Genome* const Loci);
  void Output(int iteration, const Genome* const Loci, const std::string* const PopLabels);

private:
  bool Test1, Test2;//indicators for the two tests

  double **ScoreGene;
  double **InfoGene;
  double **SumScoreGene;
  double **SumScoreGeneSq;
  double **SumInfoGene;
  int NumTestLoci;     //number of comp loci with a single locus ie those used in scalar score test for misspecified allelefreqs
  int NumCompLoci; //number of composite loci
  int Populations;           //number of populations

  double ***SumNewScore;
  double ***SumNewScoreSq;
  double ***SumNewInfo;
  int dim;

  std::ofstream allelefreqscorestream;
  std::ofstream allelefreqscorestream2;

  void UpdateScoreForMisSpecOfAlleleFreqs(int j, const double* const* phi, const std::vector<std::vector<unsigned short> > x, 
					  const double* const AlleleFreqs);
  void UpdateScoreForMisSpecOfAlleleFreqs2(const int j, const int NumberOfStates, const double* const AlleleFreqs, 
					   const int* const AlleleCounts);

  void OutputTestsForMisSpecifiedAlleleFreqs( int, const Genome* const Loci, const std::string* const PopLabels);
  void OutputTestsForMisSpecifiedAlleleFreqs2( int samples, const Genome* const Loci, const std::string* const PopLabels);
  void R_output3DarrayDimensions(ofstream* stream, const vector<int> dim, const vector<string> labels);
  string double2R( double x );
};

#endif
