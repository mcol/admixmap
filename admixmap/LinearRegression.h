// *-*-C++-*-*
#ifndef LINEARREGRESSION_H
#define LINEARREGRESSION_H 1

#include "Regression.h"
#include "IndividualCollection.h"
#include "Gaussian.h"

///Class to sample the parameters of a linear regression
class LinearRegression : public Regression{
public:
  LinearRegression();
  ~LinearRegression();
  void Initialise(unsigned RegNumber, double priorPrecision, const IndividualCollection* const, LogWriter &);
  void InitializeOutputFile(const std::vector<std::string>& CovariateLabels, unsigned NumOutcomes);
  //void Initialise(unsigned Number, const IndividualCollection* const individuals);
  double getDispersion()const;
  void OutputParams(ostream* out);
  void Update(bool sumbeta, IndividualCollection* individuals, double coolness
#ifdef PARALLEL
	      , MPI::Intracomm &Comm
#endif
	      );
  double getLogLikelihood(const IndividualCollection* const IC)const;
  double getLogLikelihoodAtPosteriorMeans(IndividualCollection *IC, int iterations);
private:
    // ** Linear Regression Objects
  double lambda0; //parameters of
  double lambda1; //Gamma prior for lambda
  double *R, *QY, *QX, *V, *betahat;
  Gaussian DrawBeta;//sampler

  void QRSolve(int dim1, int dim2, const double* a, const double* b, double* x);
  void SamplePrecision(double* lambda, const double* Y, const double* X, int NumIndivs, int NumCovars, double coolness);
  void SampleLinearRegressionParams(double* beta, const double* Y, const double* X, int NumIndivs, int NumCovars);
  void SampleLinearRegressionParametersWithAnnealing(const double* Y, const double* X, double* beta, double *lambda, 
							       double coolness);

};
#endif
