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
  void Update(IndividualCollection *individuals);
  void SumParameters(int);//should be private, part of Update
  void InitializeOutputFile(AdmixOptions *, IndividualCollection *individuals,std::string *PopulationLabels);
  void Output(int iteration, std::ofstream *LogFileStreamPtr, AdmixOptions *, IndividualCollection *individuals);
  void OutputErgodicAvg(int iteration, IndividualCollection *individuals,std::ofstream *avgstream);
  MatrixArray_d *getbeta();
  Vector_d *getlambda();
  int getNoCovariates();
  double getlambda0();
  double getDispersion(int);

private:
  int NoCovariates,AnalysisTypeIndicator;

  Matrix_d n0; //
  Gaussian DrawBeta;
  MetropolisHastings** BetaDrawArray;
  Vector_d BetaParameters;
  Vector_i acceptbeta;
  Matrix_d betan;
  Matrix_d sum;

  Vector_d lambda; //precision parameter
  double lambda0; //parameters of
  double lambda1; //prior for lambda

  MatrixArray_d beta;//regression parameters
  MatrixArray_d beta0; //
  MatrixArray_d SumBeta;//running sums (for ergodic averages) 
  Vector_d SumLambda;
 
  std::ofstream outputstream;//output to regparamfile

  static double
  lr( Vector_d & , MatrixArray_i&, MatrixArray_d& , double );
  
  static double
  dlr( Vector_d & , MatrixArray_i&, MatrixArray_d& , double );
  
  static double
  ddlr( Vector_d & , MatrixArray_i&, MatrixArray_d& , double );
  

};



























#endif
