// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   Individual.h 
 *   header file for Individual class
 *   Copyright (c) 2002-2006 LSHTM
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
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

///Class to represent an individual and update individual-level parameters
class Individual
{
public:
  Individual();

  Individual(int number, const AdmixOptions* const options, const InputData* const Data, bool undertest);
 
  ~Individual();

  void DeleteGenotypes();

  void HMMIsBad(bool loglikisbad);

  static void SetStaticMembers(Genome* const pLoci, const AdmixOptions* const options);

  static void DeleteStaticMembers();

  void setOutcome(double*);
  void setCovariates(double*);

  Sex getSex()const ;

  const double* getAdmixtureProps()const;

  void drawInitialAdmixtureProps(const vector<vector<double> > &alpha); 
  void setAdmixtureProps(const double* const, size_t);
  
  void setAdmixturePropsX(const double* const, size_t);
  
  void setGenotypesToMissing();
  void SetMissingGenotypes();

  const std::vector<hapPair > &getPossibleHapPairs(unsigned int locus)const;

  const int* getSampledHapPair(int locus)const;

  bool GenotypeIsMissing(unsigned int locus)const;//locus is a comp locus
  bool simpleGenotypeIsMissing(unsigned locus)const;//locus is a simple locus

  double getSumrho()const;

  const std::vector<double> getRho()const;

  const int *getSumLocusAncestry()const;
  const int *getSumLocusAncestryX()const;
  const std::vector<unsigned> getSumNumArrivals()const;
  const std::vector<unsigned> getSumNumArrivals_X()const;
  void getSumNumArrivals(std::vector<unsigned> *sum)const;

  double getLogLikelihood(const AdmixOptions* const options, 
			  const double* const theta, const double* const thetaX,
			  const vector<double > rho, //const vector<double> rho_X, 
			  bool updateHMM);

  double getLogLikelihood(const AdmixOptions* const , const bool forceUpdate, const bool store);
  void storeLogLikelihood(const bool setHMMAsOK); // to call if a Metropolis proposal is accepted

  double getLogLikelihoodAtPosteriorMeans(const AdmixOptions* const options);

  double getLogLikelihoodOnePop();

  double getLogPosteriorProb();

  void GetLocusAncestry(int locus, int Ancestry[2])const;
  void GetLocusAncestry(int chrm, int locus, int Ancestry[2])const;

  int GetLocusAncestry(int, int, int)const;
   
  void ResetSufficientStats();
  void SampleLocusAncestry(const AdmixOptions* const options);
  void AccumulateAncestry(int* SumAncestry);
  void UpdateScores(const AdmixOptions* const options, DataMatrix *Outcome, const DataType* const OutcomeType, 
		    DataMatrix *Covariates, double DInvLink, double dispersion,const double* const * ExpectedY);
  void SampleHapPair(unsigned chr, unsigned jj, unsigned locus, AlleleFreqs *A, bool hapmixmodel, bool anneal);
  void SampleHapPair(unsigned j, unsigned jj, unsigned locus, AlleleFreqs *A, bool hapmixmodel, bool anneal, 
		     const double* const AlleleProbs);
  void SampleJumpIndicators(bool sampleArrivals);
  void SampleRho(const AdmixOptions* const options, double rhoalpha, double rhobeta,  
		 bool updateSumLogRho);

  void SampleTheta( int iteration, double *SumLogTheta, 
		    const DataMatrix* const Outcome, 
		    const DataType* const OutcomeType, //const double* const* ExpectedY, 
		    const std::vector<double> lambda, int NumCovariates,
		    DataMatrix *Covariates, const std::vector<const double*> beta, const double* const poptheta,
		    const AdmixOptions* const options, const vector<vector<double> > &alpha, 
		    double DInvLink, double dispersion, bool RW, bool anneal);

  void SampleMissingOutcomes(DataMatrix *Outcome, const DataType* const OutcomeType, 
			     const double* const* ExpectedY, const vector<double> lambda);
 
  void FindPosteriorModes(const AdmixOptions* const options, const vector<vector<double> > &alpha,  
			  double rhoalpha, double rhobeta, //const vector<double> sigma, 
			  ofstream &modefile,
			  double *thetahat, double *thetahatX, vector<double> &rhohat); 
			  //, vector<double> &rhohatX);

  void resetStepSizeApproximator(int k);

  void Chib(int iteration, // double *SumLogLikelihood, double *MaxLogLikelihood,
	    const AdmixOptions* const options, const vector<vector<double> > &alpha, // double globalrho,
	    double rhoalpha, double rhobeta, double *thetahat, double *thetahatX,
	    vector<double> &rhohat, //vector<double> &rhohatX, 
	    chib *MargLikelihood, AlleleFreqs *A);

  static void ResetScores(const AdmixOptions* const options);
 
  static void SumScoresForLinkageAffectedsOnly(int j, double *SumAffectedsScore, 
					       double *SumAffectedsVarScore, double *SumAffectedsScore2, double *SumAffectedsInfo);
  static void SumScoresForAncestry(int j, double *SumAncestryScore, double *SumAncestryInfo, double *SumAncestryScore2,
				   double *SumAncestryVarScore);

  static void OutputLikRatios(const char* const filename, int iterations, const std::string* const PopLabels);

  double getLogPosteriorTheta()const;
  double getLogPosteriorRho()const;
  double getLogPosteriorAlleleFreqs()const;
  void SetGenotypeProbs(int j, int jj, unsigned locus, const double* const AlleleProbs);
  void SetGenotypeProbs(int j, int jj, unsigned locus, bool chibindicator);
  void AnnealGenotypeProbs(int j, const double coolness);

private:
  unsigned myNumber;//number of this individual, counting from 1
  bool IAmUnderTest;//true if not in Individual array
  bool SexIsFemale;
  bool Xdata;
  Sex sex; 
  std::vector< unsigned int > gametes;// number of gametes on each chromosome
  unsigned int X_posn;  //number of X chromosome
  unsigned NumIndGametes; // 1 if assortative mating, 2 if random mating

  std::vector<genotype> genotypes;
  std::vector<hapPair > *PossibleHapPairs;//possible haplotype pairs compatible with genotype
  double **GenotypeProbs;
  bool **GenotypesMissing;//indicators for missing genotypes at comp loci
  bool *missingGenotypes;//indicators for missing genotypes at simple loci
  std::vector<hapPair> sampledHapPairs;

  static unsigned int numChromosomes;
  static int Populations;
  static Genome *Loci;
  double *dirparams; // dirichlet parameters of full conditional for conjugate updates
  double *Theta, *ThetaX;//admixture proportions
  double* ThetaMode;
  double *SumSoftmaxTheta;
  double *ThetaProposal, *ThetaXProposal;// proposal admixture proportions

  int **LocusAncestry, *SumLocusAncestry, *SumLocusAncestry_X;
  std::vector<unsigned> SumNumArrivals;

  std::vector< double > _rho; //sum of intensities
  std::vector< double > rhoMode;
			  //std::vector< double > _rho_X;//sum of intensities for X chromosome
  std::vector<double> sumlogrho;
  double TruncationPt; // upper truncation point for sum intensities parameter rho

  double* Outcome;
  double* Covariates;

  struct {
    double value; //loglikelihood at current parameter values, annealed if coolness < 1.  Valid iff 'ready' is true
    double tempvalue; // to store values temporarily: holds unnanealed value (-energy), or value at proposed update   
    bool ready;//true iff value is the loglikelihood at the current parameter values
    bool HMMisOK;//true iff values in HMM objects correspond to current parameter values for this individual
  } logLikelihood;
  
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
  
  void SetDefaultAdmixtureProps();
  double IntegratingConst( double alpha, double beta, double a, double b )const;
  void UpdateAdmixtureForRegression( int Populations, int NumCovariates, const double* const poptheta, 
				     bool ModelIndicator, DataMatrix *Covariates);
  void Accept_Reject_Theta( double p, /*bool xdata,*/ int Populations, bool ModelIndicator, bool RW );
  //double LogAcceptanceRatioForTheta_XChrm(const std::vector<double> &sigma, int Populations );
  double LogAcceptanceRatioForRegressionModel( RegressionType RegType, bool RandomMatingModel, 
					       int Populations, int NumCovariates, 
					       const DataMatrix* const Covariates, const double* beta, 
					       //const double ExpectedY, 
					       const double Outcome, const double* const poptheta, const double lambda);
  
  void UpdateHMMInputs(unsigned int j, const AdmixOptions* const options, 
			     const double* const theta, const double* const thetaX,
			     const vector<double> rho);
  
  void ProposeTheta(const AdmixOptions* const options, /*const vector<double> sigma,*/ const vector<vector<double> > &alpha,
		    int *SumLocusAncestry, int* SumLocusAncestry_X);

  double ProposeThetaWithRandomWalk(const AdmixOptions* const options, const vector<vector<double> > &alpha);
  
  double LogPriorTheta(const double* const theta, const double* const thetaX,  
		       const AdmixOptions* const options, const vector<vector<double> > &alpha) const ;
  double LogPriorTheta_Softmax(const double* const theta, const double* const thetaX,  
			       const AdmixOptions* const options, const vector<vector<double> > &alpha) const ;
  //Computes LogPrior density in softmax basis at supplied parameter values 

  double LogPriorRho(const vector<double> rho, //const vector<double> rhoX, 
		     const AdmixOptions* const options, double rhoalpha, double rhobeta)const; 
  double LogPriorRho_LogBasis(const vector<double> rho, //const vector<double> rhoX, 
			    const AdmixOptions* const options, double rhoalpha, double rhobeta) const;

  double CalculateLogPosteriorTheta(const AdmixOptions* const options, const double* const theta, const double* const thetaX, 
				    const vector<vector<double> > &alpha)const;
  double CalculateLogPosteriorRho(const AdmixOptions* const options,  
				  const vector<double> rho, //const vector<double> rhoX,
				  double rhoalpha, double rhobeta)const;
  
  void UpdateScoreForLinkageAffectedsOnly(unsigned int locus, int Pops, int k0, bool RandomMatingModel, 
					  const vector<vector<double> > AProbs);
  void UpdateScoreForAncestry(int locus, const double* admixtureCovars, double phi, double EY, double DInvLink, const vector<vector<double> > AProbs);
  void UpdateB(double DInvLink, double dispersion, const double* admixtureCovars);
  
  void UpdateScoreTests(const AdmixOptions* const options, const double* admixtureCovars, DataMatrix *Outcome, 
			const DataType* const OutcomeType, 
			Chromosome* chrm, double DInvLink, double dispersion, const double* const* ExpectedY);

  static void SetPossibleHaplotypePairs(const vector<vector<unsigned short> > Genotype, vector<hapPair> &PossibleHapPairs);
};

#endif /* INDIVIDUAL_H */

