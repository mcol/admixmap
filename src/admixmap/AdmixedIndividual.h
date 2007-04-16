// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   AdmixedIndividual.h 
 *   header file for Individual class
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
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
#include "bcppcl/StepSizeTuner.h"
#include "chib.h"
#include <gsl/gsl_cdf.h>
#include "AffectedsOnlyTest.h"
#include "CopyNumberAssocTest.h"

class AdmixOptions;

///Class to represent an individual in an admixture model
class AdmixedIndividual : public Individual 
{
public:
  AdmixedIndividual();
  AdmixedIndividual(int number, const AdmixOptions* const options, const InputData* const Data, bool undertest);
  ~AdmixedIndividual();

  static void SetStaticMembers(Genome* const pLoci, const Options* const options);
  static void DeleteStaticMembers();
  void drawInitialAdmixtureProps(const vector<vector<double> > &alpha); 
  void SetMissingGenotypes();

  double getSumrho()const;
  const std::vector<double> getRho()const;
  double getLogLikelihood(const Options* const , const bool forceUpdate, const bool store);
  double getLogLikelihoodAtPosteriorMeans(const Options* const options);
  double getLogLikelihoodOnePop();

  void ResetSufficientStats();
  void UpdateScores(const AdmixOptions* const options, DataMatrix *Outcome, DataMatrix *Covariates, 
		    const vector<Regression*> R, AffectedsOnlyTest& affectedsOnlyTest, CopyNumberAssocTest& ancestryAssocTest);
  void SampleJumpIndicators(bool sampleArrivals);
  void SampleRho(const AdmixOptions* const options, double rhoalpha, double rhobeta,  
		 bool updateSumLogRho);
  void SampleTheta( const int iteration, double *SumLogTheta, const DataMatrix* const Outcome, 
		    const DataType* const OutcomeType, const std::vector<double> lambda, const int NumCovariates,
		    DataMatrix *Covariates, const std::vector<const double*> beta, const double* const poptheta,
		    const AdmixOptions* const options, const vector<vector<double> > &alpha, 
		    double DInvLink, const double dispersion, CopyNumberAssocTest& ancestryAssocTest,const bool RW, const bool anneal);

  void FindPosteriorModes(const AdmixOptions* const options, const vector<vector<double> > &alpha,  
			  double rhoalpha, double rhobeta, AlleleFreqs* A, ofstream &modefile);  
  void resetStepSizeApproximator(int k);
  void setChibNumerator(const AdmixOptions* const options, const vector<vector<double> > &alpha, double rhoalpha, 
	    double rhobeta, chib *MargLikelihood, AlleleFreqs *A);
  void updateChib(const AdmixOptions* const options, const vector<vector<double> > &alpha, double rhoalpha, 
	    double rhobeta, chib *MargLikelihood, AlleleFreqs *A);

  static void SumScoresForLinkageAffectedsOnly(int j, double *SumAffectedsScore, double *SumAffectedsVarScore, 
					       double *SumAffectedsScore2, double *SumAffectedsInfo);
  static void SumScoresForAncestry(int j, double *SumAncestryScore, double *SumAncestryInfo, 
				   double *SumAncestryScore2, double *SumAncestryVarScore);
  static void OutputLikRatios(const char* const filename, int iterations, const Vector_s& PopLabels);
  double getLogPosteriorTheta()const;
  double getLogPosteriorRho()const;
  double getLogPosteriorAlleleFreqs()const;

  void SetGenotypeProbs(int j, int jj, unsigned locus, const double* const AlleleProbs);
  void SetGenotypeProbs(int j, int jj, unsigned locus, bool chibindicator);
  void AnnealGenotypeProbs(int j, const double coolness);
  
  // Override functions from Individual, so AdmixedIndividual can
  // be instantiated.
  virtual const std::vector<std::vector<double> >& getUnorderedProbs(const unsigned int) const {
    throw string("AdmixedIndividual::getUnorderedProbs(const unsigned) is not implemented."); }
//   virtual void calculateUnorderedGenotypeProbs() {
//     throw string("AdmixedIndividual::calculateUnorderedGenotypeProbs() is not implemented."); }
//   virtual void calculateUnorderedGenotypeProbs(unsigned) {
//     throw string("AdmixedIndividual::calculateUnorderedGenotypeProbs(unsigned) is not implemented."); }

private:
  bool IAmUnderTest;//true if not in Individual array
  double *dirparams; // dirichlet parameters of full conditional for conjugate updates
  double *thetahat;
  double loglikhat;///< loglikelihood at posterior mode
  double *SumSoftmaxTheta;
  double *ThetaProposal;// proposal admixture proportions
  static double* ThetaThetaPrime;///< for diploid HMM updates
  static double* ThetaThetaInv;///< for diploid HMM updates
  int *SumLocusAncestry, *SumLocusAncestry_X;
  std::vector<unsigned> SumNumArrivals;
  std::vector< double > rhohat;
  std::vector<double> sumlogrho;
  double** GenotypeProbs;///<array to hold GenotypeProbs
  
  std::vector<double> logPosterior[3]; // elements 0, 1, 2 are for theta, rho, freqs
  
  //RWM sampler for individual admixture
  StepSizeTuner ThetaTuner;
  int w, NumberOfUpdates;
  double step, step0;
  
  void InitialiseSumIntensities(const AdmixOptions* const options); 
  void setAdmixtureProps(const double* const, size_t);
  ///set possible happairs, SNPs only. Required for parallel code
  void SetPossibleHaplotypePairs(const vector<vector<unsigned short> > Genotype, vector<hapPair> &PossibleHapPairs);
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
  void UpdateHMMInputs(unsigned int j, const Options* const options, 
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
  
  void UpdateScoreTests(const AdmixOptions* const options, const double* admixtureCovars, DataMatrix *Outcome, 
			Chromosome* chrm, const vector<Regression*> R, AffectedsOnlyTest& affectedsOnlyTest, CopyNumberAssocTest& ancestryAssocTest);
  double getLogLikelihood(const Options* const options, 
			  const double* const theta, const vector<double > rho, bool updateHMM);

};

#endif /* ADMIXED_INDIVIDUAL_H */

