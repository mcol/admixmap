// *-*-C++-*-*
#ifndef REGRESSION_H
#define REGRESSION_H 1

#include "matrix_d.h"
#include "vector_d.h"
#include "vector_i.h"
#include "MatrixArray_d.h"
#include "gaussian_d.h"
#include "AdmixOptions.h"
#include "IndividualCollection.h"
#include "MetropolisHastings.h"

// class Vector_d;
// class Matrix_d;
// class Vector_i;
// class Matrix_i;
// class MatrixArray_d;

class Regression{

public:
  Regression();
  Regression(int);
  ~Regression();
  void Initialise(IndividualCollection *, AdmixOptions *, LogWriter *);
  void Update(int AnalysisTypeIndicator,IndividualCollection *individuals);
  void SumParameters(int);//should be private, part of Update
  void InitializeOutputFile(AdmixOptions *, IndividualCollection *individuals,std::string *PopulationLabels);
  void Output(int iteration, std::ofstream *LogFileStreamPtr, AdmixOptions *, IndividualCollection *individuals);
  void OutputErgodicAvg(int iteration, AdmixOptions *options, IndividualCollection *individuals,std::ofstream *avgstream);
  MatrixArray_d *getbeta();
  Vector_d *getlambda();
  MatrixArray_d *getSumBeta();
  Vector_d *getSumLambda();
  int getNoCovariates();
  double getlambda0();

private:
  int NoCovariates;
  MatrixArray_d beta0; //population
  Matrix_d n0; //population
  Gaussian DrawBeta;
  MetropolisHastings** BetaDrawArray;
  Vector_d BetaParameters;
  Vector_i acceptbeta;
  Matrix_d betan;
  Matrix_d sum;
  double lambda0; //population
  double lambda1; //population
  Vector_d lambda;
  MatrixArray_d beta;//regression parameters
  MatrixArray_d SumBeta; //population
  Vector_d SumLambda; //population
  std::ofstream outputstream;//output to regparamfile

  static double
  lr( Vector_d & , MatrixArray_i&, MatrixArray_d& , double );
  
  static double
  dlr( Vector_d & , MatrixArray_i&, MatrixArray_d& , double );
  
  static double
  ddlr( Vector_d & , MatrixArray_i&, MatrixArray_d& , double );
  

};



























#endif
