// *-*-C++-*-*
#ifndef REGRESSION_H
#define REGRESSION_H 1


#include "Gaussian.h"
#include "AdmixOptions.h"
#include "IndividualCollection.h"
#include "GaussianProposalMH.h"
#include "LogWriter.h"

class Matrix_d;

class Regression{

public:
  Regression();
  Regression(int);
  ~Regression();
  void Initialise(IndividualCollection *, AdmixOptions *, std::string *Populationlabels, LogWriter *);
  void Update(bool afterBurnIn, IndividualCollection *individuals);
  void InitializeOutputFile(AdmixOptions *, IndividualCollection *individuals,std::string *PopulationLabels);
  void Output(int iteration, AdmixOptions *, IndividualCollection *individuals, LogWriter *Log);
  void OutputErgodicAvg(int iteration, IndividualCollection *individuals,std::ofstream *avgstream);
  double **getbeta();
  double *getlambda();
  int getNumCovariates();
  double getlambda0();
  double getDispersion(int);

private:
  int NumCovariates, NumOutcomeVars, AnalysisTypeIndicator;

  double *lambda; //precision parameter
  double lambda0; //parameters of
  double lambda1; //Gamma prior for lambda
  double *SumLambda;

  double **beta;//regression parameters
  Matrix_d *beta0; //beta prior mean
  Matrix_d n0; //for linear regression, beta prior variance is lambda*n0
  double **SumBeta;//running sums (for ergodic averages)

  Gaussian DrawBeta;//for linear regression
  GaussianProposalMH** BetaDrawArray;//for logistic
  double *BetaParameters;
  int *acceptbeta;
  Matrix_d betan;
  Matrix_d sum;
 
  std::ofstream outputstream;//output to regparamfile

  void SumParameters();

  static double
  lr( const double* , Matrix_i&, Matrix_d& , double );
  
  static double
  dlr( const double* , Matrix_i&, Matrix_d& , double );
  
  static double
  ddlr( const double* , Matrix_i&, Matrix_d& , double );
  

};



























#endif
