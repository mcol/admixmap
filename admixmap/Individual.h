// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   Individual.h 
 *   header file for Individual class
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */
#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H 1
#include "common.h"
#include "Genome.h"
#include "Chromosome.h"
#include "AlleleFreqs.h"
#include "LogWriter.h"
#include "StepSizeTuner.h"
#include "chib.h"
#include <gsl/gsl_cdf.h>

using namespace::std;

class AlleleFreqs;

class Individual
{
public:
  Individual();

  Individual(int number, const AdmixOptions* const options, const InputData* const Data, const Genome& Loci, 
		       const Chromosome* const * chrm);
 
  ~Individual();

  void HMMIsBad(bool loglikisbad);

  static void SetStaticMembers(const Genome* const pLoci, const AdmixOptions* const options);

  static void DeleteStaticMembers();

  Sex getSex()const ;

  const double* getAdmixtureProps()const;

  void setAdmixtureProps(const double* const, size_t);
  
  void setAdmixturePropsX(const double* const, size_t);
  
  bool GetAlleleCountsAtLocus(int complocus, int locus, int a, int* count)const;

  const std::vector<std::vector<unsigned short> > getGenotype(unsigned int locus)const;

  void setGenotypesToMissing();

  const std::vector<hapPair > &getPossibleHapPairs(unsigned int locus)const;

  const int* getSampledHapPair(int locus)const;

  bool IsMissing(unsigned int locus)const;

  double getSumrho()const;

  const std::vector<double> getRho()const;

  const int *getSumLocusAncestry()const;

  double getLogLikelihood(const AdmixOptions* const options, Chromosome **chrm, 
			  const double* const theta, const double* const thetaX,
			  const vector<double > rho, const vector<double> rho_X, bool updateHMM);
  double getLogLikelihood(const AdmixOptions* const , Chromosome**);
  double getLogLikelihoodAtPosteriorMeans(const AdmixOptions* const options, Chromosome **chrm);

  double getLogLikelihoodOnePop();

  double getLogPosteriorProb();

  void GetLocusAncestry(int chrm, int locus, int Ancestry[2])const;

  int GetLocusAncestry(int, int, int)const;
   
  double IntegratingConst( double alpha, double beta, double a, double b )const;

  void SampleParameters( double *SumLogTheta, AlleleFreqs *A, int iteration , DataMatrix *Outcome,
			 const DataType* const OutcomeType, const double* const * ExpectedY, 
			 const std::vector<double> lambda, int NumCovariates,
			 DataMatrix *Covariates, const std::vector<const double*> beta, const double *poptheta, 
			 const AdmixOptions* const options,
			 Chromosome **chrm, const vector<vector<double> > &alpha,  
			 double rhoalpha, double rhobeta, const vector<double> sigma, 
			 double DInvLink, double dispersion, bool anneal);

  void SampleTheta( int iteration, double *SumLogTheta, const DataMatrix* const Outcome, Chromosome ** C,
		    const DataType* const OutcomeType, const double* const* ExpectedY, 
		    const std::vector<double> lambda, int NumCovariates,
		    DataMatrix *Covariates, const std::vector<const double*> beta, const double* const poptheta,
		    const AdmixOptions* const options, const vector<vector<double> > &alpha, const vector<double> sigma,
		    double DInvLink, double dispersion, bool RW, bool anneal);

  void Chib(int iteration, double *SumLogLikelihood, double *MaxLogLikelihood,
	    const AdmixOptions* const options, Chromosome **chrm, const vector<vector<double> > &alpha, double globalrho,
	    double rhoalpha, double rhobeta, double *thetahat, double *thetahatX,
	    vector<double> &rhohat, vector<double> &rhohatX,
	    LogWriter& Log, chib *MargLikelihood, AlleleFreqs *A);

  static void ResetScores(const AdmixOptions* const options);
 
  static void SumScoresForLinkageAffectedsOnly(int j, double *SumAffectedsScore, 
					       double *SumAffectedsVarScore, double *SumAffectedsScore2, double *SumAffectedsInfo);
  static void SumScoresForAncestry(int j, double *SumAncestryScore, double *SumAncestryInfo, double *SumAncestryScore2,
				   double *SumAncestryVarScore);

  static void OutputLikRatios(const char* const filename, int iterations, const std::string* const PopLabels);

  double getLogPosteriorTheta()const;
  double getLogPosteriorRho()const;
  double getLogPosteriorAlleleFreqs()const;
  void SetGenotypeProbs(int j, const Chromosome* C, bool chibindicator);

private:
  unsigned myNumber;//number of this individual, counting from 1
  std::vector<genotype> genotypes;
  std::vector<hapPair > *PossibleHapPairs;//possible haplotype pairs compatible with genotype
  double **GenotypeProbs;
  bool **GenotypesMissing;
  std::vector<hapPair> sampledHapPairs;

  static unsigned int numChromosomes;
  static int Populations;
  static const Genome *Loci;
  double *Theta, *ThetaX;//admixture proportions
  double *SumSoftmaxTheta;
  double *ThetaProposal, *ThetaXProposal;// proposal admixture proportions

  int **LocusAncestry, *SumLocusAncestry, *SumLocusAncestry_X;
  unsigned SumNumArrivals[2], SumNumArrivals_X[2];

  std::vector< double > _rho; //sum of intensities
  std::vector< double > _rho_X;//sum of intensities for X chromosome
  std::vector<double> sumlogrho;
  double TruncationPt; // upper truncation point for sum intensities parameter rho

  Sex sex; 
  std::vector< unsigned int > gametes;// number of gametes on each chromosome
  unsigned int X_posn;  //number of X chromosome

  struct{
    double value;//loglikelihood at current parameter values, provided 'ready' is true
    bool ready;//true iff value is the loglikelihood at the current parameter values
    bool HMMisOK;//true iff values in HMM objects correspond to current parameter values for this individual
  }logLikelihood;

  std::vector<double> logPosterior[3];

  //RWM sampler for individual admixture
  StepSizeTuner ThetaTuner;
  int w, NumberOfUpdates;
  double step, step0;

  //score test objects, static so they can accumulate sums over individuals
  static double *AffectedsScore;
  static double *AffectedsVarScore;
  static double *AffectedsInfo;
  static double **AncestryScore;
  static double **AncestryInfo;
  static double **AncestryVarScore;
  static double **AncestryInfoCorrection;
  static double *B;//used for ancestry score test
  static double *PrevB;//holds B for previous iteration while B accumulates for this iteration
  static double *Xcov; //column matrix of covariates used to calculate B and for score test, 
                       //static only for convenience since it is reused each time

  static double *LikRatio1;
  static double *LikRatio2;

  void UpdateAdmixtureForRegression( int Populations, int NumCovariates, const double* const poptheta, 
				     bool ModelIndicator, DataMatrix *Covariates);
  void Accept_Reject_Theta( double p, bool xdata, int Populations, bool ModelIndicator, bool RW );
  double LogAcceptanceRatioForTheta_XChrm(const std::vector<double> &sigma, int Populations );
  double LogAcceptanceRatioForRegressionModel( RegressionType RegType, int TI,  bool RandomMatingModel, 
					       int Populations, int NumCovariates, 
					       const DataMatrix* const Covariates, const std::vector<const double*> beta, 
					       const double* const* ExpectedY, const DataMatrix* const Outcome, 
					       const double* const poptheta, const std::vector<double> lambda);

  void UpdateHMMForwardProbs(unsigned int j, Chromosome* const chrm, const AdmixOptions* const options, 
			  const double* const theta, const double* const thetaX,
			     const vector<double> rho, const vector<double> rhoX);

  void SampleRho(const AdmixOptions* const options, bool X_data, double rhoalpha, double rhobeta,  
		 unsigned int SumN[], unsigned int SumN_X[]);

  void ProposeTheta(const AdmixOptions* const options, const vector<double> sigma, const vector<vector<double> > &alpha);
  double ProposeThetaWithRandomWalk(const AdmixOptions* const options, Chromosome **C, const vector<vector<double> > &alpha);

  double LogPrior(const double* const theta, const double* const thetaX, const vector<double> rho, const vector<double> rhoX, 
		  const AdmixOptions* const options, const AlleleFreqs* const A, double rhoalpha, double rhobeta, 
		  const vector<vector<double> > &alpha)const;
  double CalculateLogPosteriorTheta(const AdmixOptions* const options, const double* const theta, const double* const thetaX, 
					      const vector<vector<double> > &alpha)const;
  double CalculateLogPosteriorRho(const AdmixOptions* const options,  
				  const vector<double> rho, const vector<double> rhoX,
				  double rhoalpha, double rhobeta)const;


  void UpdateScoreForLinkageAffectedsOnly(int locus, int Pops, int k0, bool RandomMatingModel, 
						    const vector<vector<double> > AProbs);
  void UpdateScoreForAncestry(int locus, double phi, double EY, double DInvLink, const vector<vector<double> > AProbs);
  void UpdateB(double DInvLink, double dispersion);

  void SampleMissingOutcomes(DataMatrix *Outcome, const DataType* const OutcomeType, 
			     const double* const* ExpectedY, const vector<double> lambda);
  void UpdateScoreTests(const AdmixOptions* const options, DataMatrix *Outcome, const DataType* const OutcomeType, Chromosome* chrm, 
			double DInvLink, double dispersion, const double* const* ExpectedY);
};

#endif /* INDIVIDUAL_H */

