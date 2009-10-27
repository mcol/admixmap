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
#include "bclib/StepSizeTuner.h"
#include "chib.h"
#include <gsl/gsl_cdf.h>
#include "AffectedsOnlyTest.h"
#include "CopyNumberAssocTest.h"
#include "PopAdmix.h"

class AdmixOptions;
class InputAdmixData;

///Class to represent an individual in an admixture model
class AdmixedIndividual : public Individual
{
public:
  AdmixedIndividual( int number, const AdmixOptions * const options, const InputAdmixData* const Data, bool undertest);
  ~AdmixedIndividual();

  static void SetStaticMembers( Genome & pLoci, const Options & options );
  static void DeleteStaticMembers();
  void drawInitialAdmixtureProps( const AlphaType & alpha );
  void SetMissingGenotypes();

  double getSumrho()const;
  const genepi::RhoType & getRho() const;
  double getLogLikelihood(const Options&, const bool forceUpdate, const bool store);
  double getLogLikelihoodAtPosteriorMeans(const Options& options);
  double getLogLikelihoodOnePop();

  void ResetSufficientStats();
  void UpdateScores(const AdmixOptions& options, bclib::DataMatrix *Outcome, bclib::DataMatrix *Covariates,
		    const vector<bclib::Regression*> & R, AffectedsOnlyTest& affectedsOnlyTest, CopyNumberAssocTest& ancestryAssocTest);
  void SampleJumpIndicators(bool sampleArrivals);
  void SampleRho(const AdmixOptions& options, double rhoalpha, double rhobeta,
		 bool updateSumLogRho);
  void SampleTheta( int iteration, double * SumLogTheta, const bclib::DataMatrix * Outcome,
		    const DataType * OutcomeType, const std::vector<double> & lambda, int NumCovariates,
		    bclib::DataMatrix * Covariates, const std::vector<const double*> & beta, const PopAdmix::PopThetaType & poptheta,
		    const AdmixOptions & options, const AlphaType & alpha,
		    double DInvLink, double dispersion, CopyNumberAssocTest & ancestryAssocTest, bool RW, bool anneal);

  void FindPosteriorModes(const AdmixOptions& options, const AlphaType &alpha,
			  double rhoalpha, double rhobeta, AlleleFreqs* A, ofstream &modefile);
  void resetStepSizeApproximator(int k);
  void setChibNumerator(const AdmixOptions& options, const AlphaType &alpha, double rhoalpha,
			double rhobeta, chib *MargLikelihood, AlleleFreqs *A);
  void updateChib(const AdmixOptions& options, const AlphaType &alpha, double rhoalpha,
		  double rhobeta, chib *MargLikelihood, AlleleFreqs *A);

  double getLogPosteriorTheta()const;
  double getLogPosteriorRho()const;
  double getLogPosteriorAlleleFreqs()const;

  void SetGenotypeProbs(int j, int jj, unsigned locus, const double* const AlleleProbs);
  void SetGenotypeProbs(int j, int jj, unsigned locus, bool chibindicator);
  void AnnealGenotypeProbs(int j, const double coolness);

  void WritePosteriorMeans(ostream& os, unsigned samples, bool globalrho)const;
  void WritePosteriorMeansLoci(ostream& os)const;
private:
  const bool IAmUnderTest; //true if not in Individual array
  bool AncestryProbs; // option LocusAncestryProbsIndicator
  double *dirparams; // dirichlet parameters of full conditional for conjugate updates
  AdmixtureProportions thetahat;
  double loglikhat;///< loglikelihood at posterior mode
  double *SumSoftmaxTheta;
  AdmixtureProportions ThetaProposal; ///< proposal admixture proportions
  int *SumLocusAncestry, *SumLocusAncestry_X;
  std::vector<unsigned> SumNumArrivals;
  genepi::RhoType rhohat;
  genepi::RhoType sumlogrho;
  double** GenotypeProbs;///<array to hold GenotypeProbs
  double* SumProbs; // array to accumulate sums of unordered hidden state probs
  std::vector<double> logPosterior[3]; // elements 0, 1, 2 are for theta, rho, freqs

  //RWM sampler for individual admixture
  bclib::StepSizeTuner ThetaTuner;
  int w, NumberOfUpdates;
  double step, step0;

  AdmixedIndividual();
  void InitialiseSumIntensities(const AdmixOptions* const options);
  void setAdmixtureProps(const AdmixtureProportions& rhs);
  ///set possible happairs, SNPs only.
  void SetPossibleHaplotypePairs(const vector<vector<unsigned short> > Genotype, vector<hapPair> &PossibleHapPairs);
  //void setAdmixturePropsX(const double* const, size_t);
  const int *getSumLocusAncestry()const;
  const int *getSumLocusAncestryX()const;
  const std::vector<unsigned> getSumNumArrivals()const;
  const std::vector<unsigned> getSumNumArrivals_X()const;
  void getSumNumArrivals(std::vector<unsigned> *sum)const;
  void UpdateAdmixtureForRegression( int Populations, int NumCovariates, const PopAdmix::PopThetaType & poptheta,
				     bool ModelIndicator, bclib::DataMatrix *Covariates);
  void Accept_Reject_Theta( double p, int Populations, bool ModelIndicator, bool RW );
  double LogAcceptanceRatioForRegressionModel( RegressionType RegType, bool RandomMatingModel,
					       int Populations, int NumCovariates,
					       const bclib::DataMatrix * Covariates, const double * beta,
					       const double Outcome, const PopAdmix::PopThetaType & poptheta, double lambda);
  void UpdateHMMInputs(unsigned int j, const Options& options,
                       const AdmixtureProportions& theta,
                       const genepi::RhoType& rho);
  void ProposeTheta(const AdmixOptions& options, const AlphaType &alpha,
		    int *SumLocusAncestry, int* SumLocusAncestry_X);
  double ProposeThetaWithRandomWalk( const AdmixOptions & options, const AlphaType & alpha );
  double LogPriorTheta_Softmax(const AdmixtureProportions& theta,
			       const AdmixOptions& options,
                               const AlphaType &alpha) const;
  double LogPriorRho_LogBasis( const genepi::RhoType & rho, const AdmixOptions& options,
			      double rhoalpha, double rhobeta) const;
  double LogPosteriorTheta_Softmax(const AdmixOptions& options,
                                   const AdmixtureProportions& theta,
                                   const AlphaType& alpha) const;
  double LogPosteriorRho_LogBasis(const AdmixOptions& options, const genepi::RhoType & rho,
				  double rhoalpha, double rhobeta)const;

  void UpdateScoreTests(const AdmixOptions& options, const double* admixtureCovars, bclib::DataMatrix *Outcome,
			Chromosome* chrm, const vector<bclib::Regression*> R, AffectedsOnlyTest& affectedsOnlyTest,
			CopyNumberAssocTest& ancestryAssocTest);
  double getLogLikelihood(const Options& options,
                          const AdmixtureProportions& theta,
                          const genepi::RhoType& rho, bool updateHMM);
  void getPosteriorMeans( double * ThetaMean, genepi::RhoType & rhoMean /* output parameter? */, unsigned samples ) const;


    // ====== DEBUGGING METHODS (overridden from PedBase) ======
    #if PEDBASE_DEBUG_METHODS
	virtual void dumpTheta( const char * prefix ) const;
    #endif
};

#endif /* ADMIXED_INDIVIDUAL_H */
