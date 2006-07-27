// *-*-C++-*-*
#ifndef COXREGRESSION_H
#define COXREGRESSION_H 1

#include "Regression.h"
#include "IndividualCollection.h"
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
  std::vector<int>::const_iterator startTimes;
  std::vector<int>::const_iterator endTimes;
  std::vector<int>::const_iterator endpoints;
  std::vector<double>::const_iterator HazardRates;
  std::vector<unsigned>::const_iterator events;
  double c;///< precision of mu 
  double mu;///<initial guess at average hazard rate per unit time
  double coolness;  ///<for annealing

}CoxBetaArgs;

///class to sample the parameters of a Cox regression
class CoxRegression : public Regression{
public:
  CoxRegression(const DataMatrix& CoxData);
  ~CoxRegression();
  void Initialise(unsigned RegNumber, double priorPrecision, const IndividualCollection* const, LogWriter &);
  //void Initialise(unsigned Number, const IndividualCollection* const individuals);
  void ReadData(const DataMatrix& CoxData);
  double DerivativeInverseLinkFunction(unsigned i)const;
  double getDispersion()const;
  void OutputParams(ostream* out);
  void OutputExpectedY();
  void Update(bool sumbeta, const std::vector<double>& Outcome, const double* const Covariates, double coolness
#ifdef PARALLEL
	      , MPI::Intracomm &Comm
#endif
	      );
  double getLogLikelihood(const std::vector<double>& Outcome)const;
  double getLogLikelihood(const double* const _beta, const std::vector<double>& _HazardRates, 
			  const std::vector<double>& Outcome)const;
  double getLogLikelihoodAtPosteriorMeans(int iterations, const std::vector<double>& Outcome);
private:
  GaussianProposalMH* BetaSampler;//to sample regression parameters
  CoxBetaArgs BetaParameters;
  int acceptbeta;
  std::vector<int> startTimes;
  std::vector<int> endTimes;
  std::vector<int> endpoints;//endpoints of subintervals
  std::vector<unsigned> events;//counts of events
  std::vector<double> HazardRates;
  static const double* EY;

  bool atRisk(unsigned ind, unsigned interval)const;

  static bool atRisk(const std::vector<int>::const_iterator start, 
		     const std::vector<int>::const_iterator finish, const std::vector<int>::const_iterator endpts);
  unsigned intervalLength(unsigned t)const;
  unsigned numFailures(unsigned ind, unsigned interval)const;

  static void getExpectedOutcome(const double* const beta, const double* const X, double* EY, int n, int d);
  static void getExpectedOutcome(const double* const beta, const double* const X, double* Y, int n, int dim, int index, double betaj);
  static void SetExpectedY(const double* const Covariates, const double* const beta, double* Ey);

  static double lr( const double beta, const void* const vargs );
  
  static double dlr( const double beta, const void* const vargs );
  
  static double ddlr( const double beta, const void* const vargs );

  CoxRegression();
};
#endif
