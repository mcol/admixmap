// *-*-C++-*-*
#ifndef COXREGRESSION_H
#define COXREGRESSION_H 1

#include "Regression.h"
#include "GaussianProposalMH.h"

///Struct to hold arguments for sampling logistic regression parameters
typedef struct{
  int n;  ///< number of individuals
  int NumIntervals; ///< number of intervals
  int d;  ///< numberof covariates
  int index;///< index of current parameter
  double beta0;  ///< prior mean
  double priorprecision;  ///< prior precision 
  const double* Covariates;
  const double* beta;  ///< regression parameters
  std::vector<double> IntervalLengths;
  std::vector<bool> atRisk;
  std::vector<double> HazardRates;
  std::vector<unsigned> events;
  double coolness;  ///<for annealing

}CoxBetaArgs;

///class to sample the parameters of a Cox regression
class CoxRegression : public Regression{
public:
  CoxRegression();
  ~CoxRegression();
  void Initialise(unsigned Number, double priorPrecision, const DataMatrix& Covars, const DataMatrix& Outcome, 
		  LogWriter &Log);
  void InitializeOutputFile(const std::vector<std::string>& CovariateLabels, unsigned NumOutcomes);
  //void Initialise(unsigned Number, const IndividualCollection* const individuals);
  void ReadData(const DataMatrix& CoxData);
  double DerivativeInverseLinkFunction(unsigned i)const;
  double getDispersion()const;
  void OutputParams(std::ostream* out);
  void Update(bool sumbeta, const std::vector<double>& Outcome, double coolness);
  double getLogLikelihood(const std::vector<double>& Outcome)const;
  double getLogLikelihood(const double* const _beta, const std::vector<double>& _HazardRates, 
			  const std::vector<double>& Outcome)const;
  double getLogLikelihoodAtPosteriorMeans(int iterations, const std::vector<double>& Outcome);
private:
  GaussianProposalMH* BetaSampler;//to sample regression parameters
  CoxBetaArgs BetaParameters;
  int acceptbeta;
  double c;///< precision of mu 
  double mu;///<initial guess at average hazard rate per unit time
  std::vector<double> sum_nr;///< sumof n_it * r_it over individuals
  static const double* EY;

  static void getExpectedOutcome(const double* const beta, const double* const X, double* EY, int n, int d);
  static void getExpectedOutcome(const double* const beta, const double* const X, double* Y, int n, int dim, int index, double betaj);
  static void SetExpectedY(const double* const Covariates, const double* const beta, double* Ey);

  static double lr( const double beta, const void* const vargs );
  
  static double dlr( const double beta, const void* const vargs );
  
  static double ddlr( const double beta, const void* const vargs );

  void plotloglikelihood(int j, const double* Covariates);
};
#endif
