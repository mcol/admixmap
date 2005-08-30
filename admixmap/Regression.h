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
   ~Regression();
  void Initialise(unsigned RegNumber, IndividualCollection *, AdmixOptions *, LogWriter *);
  void SetExpectedY(IndividualCollection *IC);
  void Update(bool afterBurnIn, IndividualCollection *individuals);
  static void OpenOutputFile(AdmixOptions *options, IndividualCollection *individuals,std::string *PopulationLabels, LogWriter *Log);  
  static void InitializeOutputFile(AdmixOptions *, IndividualCollection *individuals,std::string *PopulationLabels);
  void Output(int iteration, AdmixOptions *, LogWriter *Log);
  void OutputErgodicAvg(int iteration, std::ofstream *avgstream);
  double *getbeta();
  double getlambda();
  int getNumCovariates();
  double getDispersion();

private:
  int NumCovariates, NumOutcomeVars, AnalysisTypeIndicator;
  RegressionType RegType;
  unsigned RegNumber;

  double *beta;//regression parameters
  double *SumBeta;//running sums (for ergodic averages)
  double lambda; //precision parameter
  double SumLambda;

  // ** Linear Regression Objects
  double lambda0; //parameters of
  double lambda1; //Gamma prior for lambda
  Matrix_d beta0; //beta prior mean
  Matrix_d n0; //for linear regression, beta prior variance is lambda*n0
  Gaussian DrawBeta;//sampler
  Matrix_d betan;

  // ** Logistic Regression Objects
  GaussianProposalMH** BetaDrawArray;
  double *BetaParameters;
  int acceptbeta;
  Matrix_d sum;
  double *aCovariates;
  int *dims;
 
  static std::ofstream outputstream;//output to regparamfile

  void SumParameters();

  static double lr( const double* const, const int* const, const double* const, const double );
  
  static double dlr( const double* const, const int* const, const double* const, const double );
  
  static double ddlr( const double* const, const int* const, const double* const, const double );
  

};



























#endif
