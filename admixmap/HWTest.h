// *-*-C++-*-*
#ifndef HWTEST_H
#define HWTEST_H 1

#include "common.h"
#include "AdmixOptions.h"
#include "Genome.h"
#include "IndividualCollection.h"

class LogWriter;

class HWTest{

public:
  HWTest();

  void Initialise(AdmixOptions *options, int nloci, LogWriter *Log);

  void Output(bool IsPedFile, Matrix_s locusdata);

  void Update(IndividualCollection *IC, Chromosome **C, Genome *Loci);

  ~HWTest();

private:
  int NumInd, NumLoci, samples;
  //dmatrix score;
  //dmatrix sumscore;
  //dmatrix sumscore2;
  //dmatrix suminfo;
  double *sumscore;
  double *sumscore2;
  double *suminfo;

  std::ofstream outputfile;

  void R_output3DarrayDimensions(ofstream* stream,vector<int> dim,vector<string> labels);
  string double2R( double x );
  /**
   *  UNIMPLEMENTED: to avoid undesired copying.
   */    
  HWTest(const HWTest&);
  void operator=(const HWTest&);
};






#endif /* !defined HWTEST_H */
