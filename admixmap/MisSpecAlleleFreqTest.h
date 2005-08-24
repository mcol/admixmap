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

  void Initialise(AdmixOptions *options, Genome *Loci, LogWriter *Log );
  void Update(IndividualCollection *individuals, AlleleFreqs *A, Genome *Loci);
  void Output(int iteration, Genome *Loci,  std::string * PopLabels);

private:
  bool Test1, Test2;//indicators for the two tests

  Matrix_d *ScoreGene;
  Matrix_d *InfoGene;
  Matrix_d *SumScoreGene;
  Matrix_d *SumScoreGeneSq;
  Matrix_d *SumInfoGene;
  int NumTestLoci;     //number of comp loci with a single locus ie those used in scalar score test for misspecified allelefreqs
  int NumCompLoci; //number of composite loci
  int Populations;           //number of populations

  Matrix_d **SumNewScore;
  Matrix_d **SumNewScoreSq;
  Matrix_d **SumNewInfo;

  std::ofstream allelefreqscorestream;
  std::ofstream allelefreqscorestream2;

  void UpdateScoreForMisSpecOfAlleleFreqs(int j, double** phi, unsigned short **x, double* AlleleFreqs);
  void UpdateScoreForMisSpecOfAlleleFreqs2(const int j, const int NumberOfStates, double* AlleleFreqs, 
					   int *AlleleCounts);

  void OutputTestsForMisSpecifiedAlleleFreqs( int, Genome *Loci,  std::string * PopLabels);
  void OutputTestsForMisSpecifiedAlleleFreqs2( int samples, Genome *Loci, std::string * PopLabels);
  void R_output3DarrayDimensions(ofstream* stream,vector<int> dim,vector<string> labels);
  string double2R( double x );
};

#endif
