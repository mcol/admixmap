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

  static void SetStaticMembers(int nchr, Genome *pLoci);

  static void ResetStaticSums();

  static void DeleteStaticMembers();

  int getSex();

  Matrix_d& getAdmixtureProps();

  void setAdmixtureProps(Matrix_d);
  
  Matrix_d& getAdmixturePropsX();

  void setAdmixturePropsX(Matrix_d);
  
  unsigned short **getGenotype(unsigned int locus);

  std::vector<hapPair > &getPossibleHapPairs(unsigned int locus);

  bool IsMissing(unsigned int locus);

  static int *getSumXi();

  static int getSumXi(int j);

  static double getSumrho0();

  double getSumrho();

  std::vector<double> getRho();

  double getLogLikelihood(AdmixOptions*,AlleleFreqs*,Chromosome**);

  double getLogLikelihoodOnePop(AlleleFreqs *);

  double getLogPosteriorProb();

  Vector_i GetLocusAncestry( int, int );
  //int *GetLocusAncestry(int, int);

  int GetLocusAncestry(int, int, int);
   
  double getLogLikelihood(AdmixOptions*,AlleleFreqs*,Chromosome **, Matrix_d, std::vector<double>,Matrix_d, std::vector<double>);
  double getLogLikelihoodXOnly(AdmixOptions*,AlleleFreqs*,Chromosome**, Matrix_d, std::vector<double>);
  double IntegratingConst( double alpha, double beta, double a, double b );

  void SampleParameters( int i, Vector_d *SumLogTheta, AlleleFreqs *A, int iteration , Matrix_d *Outcome,
			 int NumOutcomes, Vector_i &OutcomeType, double **ExpectedY, Vector_d &lambda, int NoCovariates,
			 Matrix_d &Covariates0, double **beta, Vector_d &poptheta, AdmixOptions* options,
			 Chromosome **chrm, vector<Vector_d> alpha, bool _symmetric, vector<bool> _admixed, 
			 double rhoalpha, double rhobeta, vector<double> sigma, 
			 double DInvLink, double dispersion);

 void OnePopulationUpdate( int i, Matrix_d *Outcome, int NumOutcomes, Vector_i &OutcomeType, double **ExpectedY, Vector_d &lambda,
			   int AnalysisTypeIndicator);

  void ChibLikelihood(int iteration, double *LogLikelihood, double *SumLogLikelihood, double *MaxLogLikelihood,
		      AdmixOptions *options, Chromosome **chrm, vector<Vector_d> alpha,  
		      vector<bool> _admixed, double rhoalpha, double rhobeta, Matrix_d &thetahat, Matrix_d &thetahatX,
		      vector<double> &rhohat, vector<double> &rhohatX,
		      std::ofstream *LogFileStreamPtr, chib *MargLikelihood, AlleleFreqs *A);

  static void InitialiseAffectedsOnlyScores(int L, int K);
  static void InitialiseAncestryScores(int L, int K);
  static void ResetScores(AdmixOptions *options);
 
  static void SumScoresForLinkageAffectedsOnly(int j,int Populations, Matrix_d *SumAffectedsScore, 
					       Matrix_d *SumAffectedsVarScore,Matrix_d *SumAffectedsScore2, Matrix_d *SumAffectedsInfo);
  static void SumScoresForAncestry(int j, int Populations,  
				      Matrix_d *SumAncestryScore, Matrix_d *SumAncestryInfo, Matrix_d *SumAncestryScore2,
				      Matrix_d *SumAncestryVarScore);


private:
  unsigned short ***genotypes;
  std::vector<hapPair > *PossibleHapPairs;//possible haplotype pairs compatible with genotype

  static int *sumxi;//sum of jump indicators over individuals, gametes
  static double Sumrho0;//? sum of distances between loci where there are no arrivals, summed over individuals

  static unsigned int numChromosomes;
  static Genome *Loci;
  Matrix_d Theta, ThetaX;//admixture proportions
  Matrix_d ThetaProposal, ThetaXProposal;// proposal admixture proportions

  Matrix_d AdmixtureHat;
  Matrix_d XAdmixtureHat;

  Matrix_i *LocusAncestry;
  Matrix_i SumLocusAncestry, SumLocusAncestry_X;

  std::vector< double > _rho;
  std::vector< double > _rho_X;
  std::vector< double > _rhoHat;
  std::vector< double > _rhoHat_X;
  // f0 and f1 are arrays of scalars of the form exp - rho*x, where x is distance between loci
  // required to calculate transition matrices 
  double *f[2]; 

  double LogPosterior;
  short unsigned int sex; // 0 = missing, 1 = male, 2 = female 
  std::vector< unsigned int > gametes;
  unsigned int X_posn;
  double TruncationPt; // upper truncation point for sum intensities parameter rho

  //score test objects, static so they can accumulate sums over individuals
  static Matrix_d AffectedsScore;
  static Matrix_d AffectedsVarScore;
  static Matrix_d AffectedsInfo;
  static Matrix_d *AncestryScore;
  static Matrix_d *AncestryInfo;
  static Matrix_d AncestryVarScore;
  static Matrix_d AncestryInfoCorrection;
  static Matrix_d B;//used for ancestry score test
  static Matrix_d PrevB;//holds B for previous iteration while B accumulates for this iteration
  static Matrix_d Xcov; //column matrix of covariates used to calculate B and for score test, 
                       //static only for convenience since it is reused each time

  void UpdateAdmixtureForRegression( int i,int Populations, int NoCovariates, Vector_d &poptheta, bool ModelIndicator,
				     Matrix_d *Covariates0);
  void Accept_Reject_Theta( double p, bool xdata, int Populations, bool ModelIndicator );
  double AcceptanceProbForTheta_XChrm(std::vector<double> &sigma, int Populations );
  double AcceptanceProbForTheta_LogReg( int i, int TI, bool ModelIndicator,int Populations, 
					int NoCovariates, Matrix_d &Covariates0, double **beta, double **ExpectedY,
					Matrix_d *Outcome, Vector_d &poptheta);
  double AcceptanceProbForTheta_LinearReg( int i, int TI, bool ModelIndicator,int Populations,
					   int NoCovariates, Matrix_d &Covariates0, double **beta, double **ExpectedY,
					   Matrix_d *Outcome, Vector_d &poptheta, Vector_d &lambda);

  bool UpdateForBackProbs(unsigned int j, Chromosome *chrm, AdmixOptions *options, bool randomAlleleFreqs);

  void SumAncestry(unsigned int j, Chromosome *chrm);

  void SampleRho(bool XOnly, bool RandomMatingModel, bool X_data, double rhoalpha, double rhobeta, double L, double L_X, 
		 unsigned int SumN[], unsigned int SumN_X[]);
  void SampleTheta( int i, Vector_d *SumLogTheta, Matrix_d *Outcome,
		    int NumOutcomes,  Vector_i &OutcomeType, double **ExpectedY, Vector_d &lambda, int NoCovariates,
		    Matrix_d &Covariates0, double **beta, Vector_d &poptheta,
		    AdmixOptions* options, vector<Vector_d> alpha, vector<double> sigma);

  void ProposeTheta(AdmixOptions *options, vector<double> sigma, vector<Vector_d> alpha);
  void CalculateLogPosterior(AdmixOptions *options, bool isX_data, vector<Vector_d> alpha, 
						 bool _symmetric, vector<bool> _admixed, double rhoalpha, double rhobeta, double L, 
			     double L_X, unsigned int SumN[], unsigned int SumN_X[]);

  void InitializeChib(Matrix_d theta, Matrix_d thetaX, vector<double> rho, vector<double> rhoX, 
		 AdmixOptions *options, AlleleFreqs *A, Chromosome **chrm, double rhoalpha, double rhobeta, 
		 vector<Vector_d> alpha, vector<bool> _admixed, chib *MargLikelihood, std::ofstream *LogFileStreamPtr);

  void setIsMissing(vector<unsigned int >& decoded);

  void UpdateScoreForLinkageAffectedsOnly(int j,int Populations, bool ModelIndicator, Chromosome **);//could be private
  void UpdateScoreForAncestry(int j,double phi, double EY,double DInvLink, Chromosome **, int Populations);
};


#endif /* INDIVIDUAL_H */

