// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   AdmixedIndividual.h 
 *   header file for Individual class
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef ADMIXED_INDIVIDUAL_H
#define ADMIXED_INDIVIDUAL_H 1
#include "Individual.h"
#include "samplers/StepSizeTuner.h"
#include "chib.h"
#include <gsl/gsl_cdf.h>
#include "AffectedsOnlyTest.h"

///Class to represent an individual in an admixture model
class AdmixedIndividual : public Individual 
{
public:
  AdmixedIndividual();
  AdmixedIndividual(int number, const AdmixOptions* const options, const InputData* const Data, bool undertest);
  ~AdmixedIndividual();

  static void SetStaticMembers(Genome* const pLoci, const AdmixOptions* const options);
  static void DeleteStaticMembers();
  void drawInitialAdmixtureProps(const vector<vector<double> > &alpha); 

  double getSumrho()const;
  const std::vector<double> getRho()const;
  double getLogLikelihood(const AdmixOptions* const , const bool forceUpdate, const bool store);
  double getLogLikelihoodAtPosteriorMeans(const AdmixOptions* const options);
  double getLogLikelihoodOnePop();

  void ResetSufficientStats();
  void UpdateScores(const AdmixOptions* const options, DataMatrix *Outcome, DataMatrix *Covariates, 
		    const vector<Regression*> R, AffectedsOnlyTest& affectedsOnlyTest);
  void SampleJumpIndicators(bool sampleArrivals);
  void SampleRho(const AdmixOptions* const options, double rhoalpha, double rhobeta,  
		 bool updateSumLogRho);
  void SampleTheta( const int iteration, double *SumLogTheta, const DataMatrix* const Outcome, 
		    const DataType* const OutcomeType, const std::vector<double> lambda, const int NumCovariates,
		    DataMatrix *Covariates, const std::vector<const double*> beta, const double* const poptheta,
		    const AdmixOptions* const options, const vector<vector<double> > &alpha, 
		    double DInvLink, const double dispersion, const bool RW, const bool anneal);

  void FindPosteriorModes(const AdmixOptions* const options, const vector<vector<double> > &alpha,  
			  double rhoalpha, double rhobeta, AlleleFreqs* A, ofstream &modefile);  
  void resetStepSizeApproximator(int k);
  void setChibNumerator(const AdmixOptions* const options, const vector<vector<double> > &alpha, double rhoalpha, 
	    double rhobeta, chib *MargLikelihood, AlleleFreqs *A);
  void updateChib(const AdmixOptions* const options, const vector<vector<double> > &alpha, double rhoalpha, 
	    double rhobeta, chib *MargLikelihood, AlleleFreqs *A);

  static void ResetScores(const AdmixOptions* const options);
  static void SumScoresForLinkageAffectedsOnly(int j, double *SumAffectedsScore, double *SumAffectedsVarScore, 
					       double *SumAffectedsScore2, double *SumAffectedsInfo);
  static void SumScoresForAncestry(int j, double *SumAncestryScore, double *SumAncestryInfo, 
				   double *SumAncestryScore2, double *SumAncestryVarScore);
  static void OutputLikRatios(const char* const filename, int iterations, const Vector_s& PopLabels);
  double getLogPosteriorTheta()const;
  double getLogPosteriorRho()const;
  double getLogPosteriorAlleleFreqs()const;

private:
  bool IAmUnderTest;//true if not in Individual array
  double *dirparams; // dirichlet parameters of full conditional for conjugate updates
  double *thetahat;
  double loglikhat;///< loglikelihood at posterior mode
  double *SumSoftmaxTheta;
  double *ThetaProposal;// proposal admixture proportions
  int *SumLocusAncestry, *SumLocusAncestry_X;
  std::vector<unsigned> SumNumArrivals;
  std::vector< double > rhohat;
  std::vector<double> sumlogrho;
  
  std::vector<double> logPosterior[3]; // elements 0, 1, 2 are for theta, rho, freqs
  
  //RWM sampler for individual admixture
  StepSizeTuner ThetaTuner;
  int w, NumberOfUpdates;
  double step, step0;
  
  //score test objects, static so they can accumulate sums over individuals
  //  static double *AffectedsScore;
  //  static double *AffectedsVarScore;
  //  static double *AffectedsInfo;
  static double **AncestryScore;
  static double **AncestryInfo;
  static double **AncestryVarScore;
  static double **AncestryInfoCorrection;
  static double *B;//used for ancestry score test
  static double *PrevB;//holds B for previous iteration while B accumulates for this iteration
  static double *Xcov; //column matrix of covariates used to calculate B and for score test, 
                       //static only for convenience since it is reused each time
  
  //  static double *LikRatio1;
  //  static double *LikRatio2;
 
  void InitialiseSumIntensities(const AdmixOptions* const options); 
  void setAdmixtureProps(const double* const, size_t);
  //void setAdmixturePropsX(const double* const, size_t);
  const int *getSumLocusAncestry()const;
  const int *getSumLocusAncestryX()const;
  const std::vector<unsigned> getSumNumArrivals()const;
  const std::vector<unsigned> getSumNumArrivals_X()const;
  void getSumNumArrivals(std::vector<unsigned> *sum)const;
  void UpdateAdmixtureForRegression( int Populations, int NumCovariates, const double* const poptheta, 
				     bool ModelIndicator, DataMatrix *Covariates);
  void Accept_Reject_Theta( double p, int Populations, bool ModelIndicator, bool RW );
  double LogAcceptanceRatioForRegressionModel( RegressionType RegType, bool RandomMatingModel, 
					       int Populations, int NumCovariates, 
					       const DataMatrix* const Covariates, const double* beta, 
					       const double Outcome, const double* const poptheta, const double lambda);
  void UpdateHMMInputs(unsigned int j, const AdmixOptions* const options, 
			     const double* const theta, const vector<double> rho);
  void ProposeTheta(const AdmixOptions* const options, const vector<vector<double> > &alpha,
		    int *SumLocusAncestry, int* SumLocusAncestry_X);
  double ProposeThetaWithRandomWalk(const AdmixOptions* const options, const vector<vector<double> > &alpha);
  double LogPriorTheta_Softmax(const double* const theta, 
			       const AdmixOptions* const options, const vector<vector<double> > &alpha) const ;
  double LogPriorRho_LogBasis(const vector<double> rho, const AdmixOptions* const options, 
			      double rhoalpha, double rhobeta) const;
  double LogPosteriorTheta_Softmax(const AdmixOptions* const options, const double* const theta, 
				    const vector<vector<double> > &alpha)const;
  double LogPosteriorRho_LogBasis(const AdmixOptions* const options, const vector<double> rho, 
				  double rhoalpha, double rhobeta)const;
  
  //  static void UpdateScoreForLinkageAffectedsOnly(unsigned int locus, int Pops, int k0, const double* const _Theta_, 
  //					 bool RandomMatingModel, bool diploid, const vector<vector<double> > AProbs);
  void UpdateScoreForAncestry(int locus, const double* admixtureCovars, double phi, double EY, double DInvLink, 
			      const vector<vector<double> > AProbs);
  void UpdateB(double DInvLink, double dispersion, const double* admixtureCovars);
  void UpdateScoreTests(const AdmixOptions* const options, const double* admixtureCovars, DataMatrix *Outcome, 
				  Chromosome* chrm, const vector<Regression*> R, AffectedsOnlyTest& affectedsOnlyTest);
  double getLogLikelihood(const AdmixOptions* const options, 
			  const double* const theta, const vector<double > rho, bool updateHMM);
};

#endif /* ADMIXED_INDIVIDUAL_H */

