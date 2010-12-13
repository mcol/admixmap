//=============================================================================
//
// Copyright (C) 2002-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file AdmixedIndividual.h
/// Definition of the AdmixedIndividual class.
//=============================================================================

#ifndef ADMIXED_INDIVIDUAL_H
#define ADMIXED_INDIVIDUAL_H 1

#include "Individual.h"
#include "AffectedsOnlyTest.h"
#include "PopAdmix.h"
#include "bclib/StepSizeTuner.h"
#include <gsl/gsl_cdf.h>

class AdmixOptions;
class CopyNumberAssocTest;
class InputAdmixData;
class chib;


/** \addtogroup admixmap
 * @{ */


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

  /// Sample individual sumintensities (globalrho=0)
  void SampleRho(const AdmixOptions& options, double rhoalpha, double rhobeta,
		 bool updateSumLogRho);

  /// Sample individual odds ratios (globalpsi=0)
  void SamplePsi(const AdmixOptions& options,
                 const genepi::cvector<double>& priormean,
                 const genepi::cvector<double>& priorprec,
                 bool updateSumLogPsi);

  /// Sample individual admixture proportions
  void SampleTheta(int iteration, double *SumLogTheta,
                   const bclib::DataMatrix *Outcome,
                   const DataType *OutcomeType,
                   const std::vector<double>& lambda, int NumCovariates,
                   bclib::DataMatrix *Covariates,
                   const std::vector<const double*>& beta,
                   const PopAdmix::PopThetaType& poptheta,
                   const AdmixOptions& options, const AlphaType& alpha,
                   double DInvLink, double dispersion,
                   CopyNumberAssocTest& ancestryAssocTest,
                   bool RW, bool updateSumLogTheta);

  /// Use an EM algorithm to search for posterior modes of individual
  /// parameters theta and rho
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

  void SetGenotypeProbs(int chr, unsigned sizeOfChromosome, bool chibindicator);
  void AnnealGenotypeProbs(int j, const double coolness);

  void WritePosteriorMeans(ostream& os, unsigned samples, bool globalrho)const;
  void WritePosteriorMeansXChr(ostream& os, unsigned samples) const;
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

  /// Number of arrivals between each pair of adjacent loci (used only for
  /// models with individual-specific sum-intensities, that is globalrho=0)
  std::vector<unsigned> SumNumArrivals;
  genepi::RhoType rhohat;
  genepi::RhoType sumlogrho;
  double** GenotypeProbs;///<array to hold GenotypeProbs

  /// Vector to accumulate the sums of unordered hidden state probabilities
  /// (used only for LocusAncestryProbsIndicator). It's mutable so that it
  /// can be updated in WritePosteriorMeans().
  mutable genepi::cvector<double> SumProbs;

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

  /// Return the number of arrivals across the genome
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

  void UpdateScoreTests(const AdmixOptions& options,
                        const double *admixtureCovs, bclib::DataMatrix *Outcome,
                        Chromosome *chrm, const vector<bclib::Regression*>& R,
                        AffectedsOnlyTest& affectedsOnlyTest,
			CopyNumberAssocTest& ancestryAssocTest);
  double getLogLikelihood(const Options& options,
                          const AdmixtureProportions& theta,
                          const genepi::RhoType& rho, bool updateHMM);
  void getPosteriorMeansTheta(genepi::cvector<double>& ThetaMean,
                              double dSamples) const;
  void getPosteriorMeansRho(genepi::RhoType& rhoMean, double dSamples) const;


    // ====== DEBUGGING METHODS (overridden from PedBase) ======
    #if PEDBASE_DEBUG_METHODS
	virtual void dumpTheta( const char * prefix ) const;
    #endif
};


/** @} */


#endif /* ADMIXED_INDIVIDUAL_H */
