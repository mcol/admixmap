// *-*-C++-*-*
#ifndef REGRESSION_H
#define REGRESSION_H 1


#include "Gaussian.h"
#include "AdmixOptions.h"
#include "IndividualCollection.h"
#include "GaussianProposalMH.h"
#include "LogWriter.h"

//can go once the blas code in update is removed
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

typedef struct{
  int n;//number of individuals
  int d;//numberof covariates
  int index;//index of current parameter
  double beta0;//prior mean
  double lambda;//prior precision 
  double XtY;
  const double* Covariates;
  const double* beta;//regression parameters

}BetaArgs;

class Regression{

public:
  Regression();
   ~Regression();
  void Initialise(unsigned RegNumber, const IndividualCollection* const, LogWriter *);
  void SetExpectedY(IndividualCollection* IC)const;
  void Update(bool sumbeta, IndividualCollection* individuals);
  static void OpenOutputFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
			     const std::string *PopulationLabels, LogWriter *Log);  
  static void InitializeOutputFile(const AdmixOptions* const , const IndividualCollection* const individuals, 
				   const std::string* const PopulationLabels);
  void Output(int iteration, AdmixOptions *, LogWriter *Log)const;
  void OutputErgodicAvg(int iteration, std::ofstream *avgstream)const;
  const double* const getbeta() const;
  double getlambda() const ;
  int getNumCovariates()const;
  double getDispersion()const;
  double getLogLikelihood(const IndividualCollection* const IC)const;
  double getLogLikelihoodAtPosteriorMeans(IndividualCollection *IC, int iterations);

private:
  int NumCovariates, NumOutcomeVars;
  RegressionType RegType;
  unsigned RegNumber;

  double *beta;//regression parameters
  double *SumBeta;//running sums (for ergodic averages)
  double lambda; //precision parameter
  double SumLambda;

  // ** Linear Regression Objects
  double lambda0; //parameters of
  double lambda1; //Gamma prior for lambda
  double *beta0; //beta prior mean
  double *n0; //for linear regression, beta prior variance is lambda*n0
  Gaussian DrawBeta;//sampler
  double *betan;

  // ** Logistic Regression Objects
  GaussianProposalMH** BetaDrawArray;
  BetaArgs BetaParameters;
  int acceptbeta;
  double* XtY;
  const double *X;
  int *dims;
 
  static std::ofstream outputstream;//output to regparamfile

  void SumParameters();

  static double lr( const double beta, const void* const vargs );
  
  static double dlr( const double beta, const void* const vargs );
  
  static double ddlr( const double beta, const void* const vargs );
  

};



























#endif
