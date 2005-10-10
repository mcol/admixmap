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

  Individual(int i,AdmixOptions*,InputData *Data,Genome&,Chromosome **);
 
  ~Individual();

  static void SetStaticMembers(Genome *pLoci, AdmixOptions *options);

  static void DeleteStaticMembers();

  Sex getSex();

  double *getAdmixtureProps();

  void setAdmixtureProps(double *, size_t);
  
  void setAdmixturePropsX(double *, size_t);
  
  unsigned short **getGenotype(unsigned int locus);

  std::vector<hapPair > &getPossibleHapPairs(unsigned int locus);

  bool IsMissing(unsigned int locus);

  double getSumrho();

  std::vector<double> getRho();

  int *getSumLocusAncestry();

  double getLogLikelihood(AdmixOptions*, Chromosome**);
  double getLogLikelihoodAtPosteriorMeans(AdmixOptions* options, Chromosome **chrm);
  double getLogLikelihood( AdmixOptions* options, Chromosome **chrm, double *theta, double* theta, 
			   vector<double > rho, vector<double> rho_X, bool chibindicator);

  double getLogLikelihoodOnePop(bool chibindicator);

  double getLogPosteriorProb();

  void GetLocusAncestry(int chrm, int locus, int Ancestry[2]);

  int GetLocusAncestry(int, int, int);
   
  double IntegratingConst( double alpha, double beta, double a, double b );

  void SampleParameters( int i, double *SumLogTheta, double *LogLikelihood, AlleleFreqs *A, int iteration , DataMatrix *Outcome,
			 int NumOutcomes, DataType* OutcomeType, double **ExpectedY, double *lambda, int NoCovariates,
			 DataMatrix *Covariates, double **beta, const double *poptheta, AdmixOptions* options,
			 Chromosome **chrm, vector<vector<double> > &alpha,  
			 double rhoalpha, double rhobeta, vector<double> sigma, 
			 double DInvLink, double dispersion);

  void SampleTheta( int i, int iteration, double *SumLogTheta, DataMatrix *Outcome, Chromosome **C,
		    int NumOutcomes, DataType* OutcomeType, double **ExpectedY, double *lambda, int NoCovariates,
		    DataMatrix *Covariates, double **beta, const double *poptheta,
		    AdmixOptions* options, vector<vector<double> > &alpha, vector<double> sigma, double, double, bool);

 void OnePopulationUpdate( int i, DataMatrix *Outcome, int NumOutcomes, DataType* OutcomeType, double **ExpectedY, double *lambda,
			   Chromosome **chrm, AlleleFreqs *A );

  void Chib(int iteration, double *LogLikelihood, double *SumLogLikelihood, double *MaxLogLikelihood,
		      AdmixOptions *options, Chromosome **chrm, vector<vector<double> > &alpha, double globalrho,
		      double rhoalpha, double rhobeta, double *thetahat, double *thetahatX,
		      vector<double> &rhohat, vector<double> &rhohatX,
		      LogWriter *Log, chib *MargLikelihood, AlleleFreqs *A);

  static void ResetScores(AdmixOptions *options);
 
  static void SumScoresForLinkageAffectedsOnly(int j, double *SumAffectedsScore, 
					       double *SumAffectedsVarScore, double *SumAffectedsScore2, double *SumAffectedsInfo);
  static void SumScoresForAncestry(int j, double *SumAncestryScore, double *SumAncestryInfo, double *SumAncestryScore2,
				   double *SumAncestryVarScore);

  static void OutputLikRatios(const char *filename, int iterations, std::string *PopLabels);

private:
  unsigned short ***genotypes;
  std::vector<hapPair > *PossibleHapPairs;//possible haplotype pairs compatible with genotype

  static unsigned int numChromosomes;
  static int Populations;
  static Genome *Loci;
  double *Theta, *ThetaX;//admixture proportions
  double *SumSoftmaxTheta;
  double *ThetaProposal, *ThetaXProposal;// proposal admixture proportions

  int **LocusAncestry, *SumLocusAncestry, *SumLocusAncestry_X;
  unsigned SumN[2], SumN_X[2];

  std::vector< double > _rho; //sum of intensities
  std::vector< double > _rho_X;//sum of intensities for X chromosome
  std::vector<double> sumlogrho;

  //double LogPosterior;
  Sex sex; // 0 = missing, 1 = male, 2 = female 
  std::vector< unsigned int > gametes;// number of gametes on each chromosome
  unsigned int X_posn;  //number of X chromosome
  double TruncationPt; // upper truncation point for sum intensities parameter rho

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

  void UpdateAdmixtureForRegression( int i,int Populations, int NoCovariates, const double *poptheta, bool ModelIndicator,
				     DataMatrix *Covariates);
  void Accept_Reject_Theta( double p, bool xdata, int Populations, bool ModelIndicator, bool RW );
  double AcceptanceProbForTheta_XChrm(std::vector<double> &sigma, int Populations );
  double LogAcceptanceRatioForRegressionModel( int i, RegressionType RegType, int TI,  bool RandomMatingModel, int Populations,
					       int NoCovariates, DataMatrix *Covariates, double **beta, double **ExpectedY,
					       DataMatrix *Outcome, const double *poptheta, double *lambda);

  bool UpdateForBackProbs(unsigned int j, Chromosome *chrm, AdmixOptions *options, 
			  double* theta, double *thetaX,
			  vector<double> rho, vector<double> rhoX, bool chibindicator, bool calcbackprobs);

  void SumAncestry(unsigned int j, Chromosome *chrm);

  void SampleRho(bool XOnly, bool RandomMatingModel, bool X_data, double rhoalpha, double rhobeta,  
		 unsigned int SumN[], unsigned int SumN_X[]);

  void ProposeTheta(AdmixOptions *options, vector<double> sigma, vector<vector<double> > &alpha);
  double ProposeThetaWithRandomWalk(AdmixOptions *options, Chromosome **C, vector<vector<double> > &alpha);

  double LogPrior(double *theta, double *thetaX, vector<double> rho, vector<double> rhoX, 
				AdmixOptions *options, AlleleFreqs *A, double rhoalpha, double rhobeta, 
				vector<vector<double> > &alpha);
  double CalculateLogPosterior(AdmixOptions *options, double *theta, double *thetaX, vector<double> rho, vector<double> rhoX,
				       vector<vector<double> > &alpha, double rhoalpha, double rhobeta);


  void setIsMissing(vector<unsigned int >& decoded);

  void UpdateScoreForLinkageAffectedsOnly(int j, bool ModelIndicator, Chromosome **);
  void UpdateScoreForAncestry(int j,double phi, double EY,double DInvLink, Chromosome **);
  void UpdateB(double DInvLink, double dispersion);
};


#endif /* INDIVIDUAL_H */

