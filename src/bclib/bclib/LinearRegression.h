// *-*-C++-*-*
#ifndef LINEARREGRESSION_H
#define LINEARREGRESSION_H 1

#include "bclib/Regression.h"
#include "bclib/Gaussian.h"

///Class to sample the parameters of a linear regression
class LinearRegression : public Regression{
public:
  LinearRegression(unsigned Number, double priorPrecision, const DataMatrix& Covars, const DataMatrix& Outcome, 
		   LogWriter &Log);
  ~LinearRegression();
  void InitializeOutputFile(const std::vector<std::string>& CovariateLabels, unsigned NumOutcomes);

  double getDispersion()const;
  double DerivativeInverseLinkFunction(unsigned i)const;
  void OutputParams(std::ostream* out)const;
  void OutputErgodicAvg(int samples, std::ofstream *avgstream)const;
  void Update(bool sumbeta, const std::vector<double>& Outcome, double coolness);
  double getLogLikelihood(const std::vector<double>& Outcome)const;
  double getLogLikelihoodAtPosteriorMeans(int iterations, const std::vector<double>& Outcome);
private:
  // ** Linear Regression Objects
  double lambda0; //parameters of
  double lambda1; //Gamma prior for lambda
  double *R, *QY, *QX, *V, *betahat;
  Gaussian DrawBeta;//sampler

  void SetExpectedY(const double* const _beta);
  void QRSolve(int dim1, int dim2, const double* a, const double* b, double* x);
  void SamplePrecision(double* lambda, const double* Y, const double* X, int NumIndivs, int NumCovars, double coolness);
  void SampleLinearRegressionParams(double* beta, const double* Y, const double* X, int NumIndivs, int NumCovars);
  void SampleLinearRegressionParametersWithAnnealing(const double* Y, const double* X, double* beta, double *lambda, 
							       double coolness);

  LinearRegression();
  void Initialise(double priorPrecision, const DataMatrix& Covars, const DataMatrix& Outcome, 
		  LogWriter &Log);

};
#endif
