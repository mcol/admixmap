// *-*-C++-*-*
#ifndef LOGISTICREGRESSION_H
#define LOGISTICREGRESSION_H 1

#include "Regression.h"
#include "IndividualCollection.h"
#include "GaussianProposalMH.h"

///Struct to hold arguments for sampling logistic regression parameters
typedef struct{
  ///number of individuals
  int n;
  ///numberof covariates
  int d;
  ///index of current parameter
  int index;
  ///prior mean
  double beta0;
  ///prior precision 
  double priorprecision;
  double XtY;
  const double* Covariates;
  ///regression parameters
  const double* beta;
  ///for annealing
  double coolness;

}BetaArgs;

///class to sample the parameters of a logistic regression
class LogisticRegression : public Regression{
public:
  LogisticRegression();
  ~LogisticRegression();
  void Initialise(unsigned RegNumber, double priorPrecision, const IndividualCollection* const, LogWriter &);
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
  // ** Logistic Regression Objects
  GaussianProposalMH** BetaDrawArray;
  BetaArgs BetaParameters;
  int acceptbeta;
  int *dims;
 
  static double lr( const double beta, const void* const vargs );
  
  static double dlr( const double beta, const void* const vargs );
  
  static double ddlr( const double beta, const void* const vargs );

};
#endif
