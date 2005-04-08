// *-*-C++-*-*
#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H 1
#include "common.h"
#include "IndividualVisitor.h"
#include "Genome.h"
#include "Chromosome.h"
//#include "AlleleFreqs.h"
#include "chib.h"
#include <gsl/gsl_cdf.h>

using namespace::std;

class AlleleFreqs;
class Individual
{
public:
  Individual();

  Individual(AdmixOptions*,const Vector_s& data,Genome&,Chromosome **);

  ~Individual();

  int getSex();

  Matrix_d& getAdmixtureProps();

  void setAdmixtureProps(Matrix_d);
  
  Matrix_d& getAdmixturePropsX();

  void setAdmixturePropsX(Matrix_d);
  
  //   std::vector< std::vector<unsigned int> >&
  //   getGenotype();

  std::vector<unsigned int>&  getGenotype(unsigned int locus);

  Vector_i getPossibleHaplotypes(unsigned int locus);

  std::vector< std::vector<unsigned int> >&  IsMissing();

  std::vector<unsigned int>&
  IsMissing(unsigned int locus);

  void setGenotype(unsigned int locus,std::vector<unsigned int> genotype);

  std::vector<bool>&  getXi(unsigned int locus);

  const std::vector< std::vector<bool> >&  getXi();

  Vector_i getSumXi();

  double getSumrho0();

  double getSumrho();

  std::vector<double> getRho();

  Matrix_d getAncestryProbs( int , int);
  Matrix_d getAncestryProbs( int );

  double getLogLikelihood(AdmixOptions*,AlleleFreqs*,Chromosome**);

  double getLogLikelihoodOnePop(AlleleFreqs *);

  double getLogPosteriorProb();

  Vector_i GetLocusAncestry( int, int );
   
  Matrix_d SampleParameters(int, AdmixOptions*, AlleleFreqs*, Chromosome**, std::vector<Vector_d>,
			    bool, std::vector<bool>, double, double, int, std::vector<double>, Matrix_d& );
  double getLogLikelihood(AdmixOptions*,AlleleFreqs*,Chromosome **, Matrix_d, std::vector<double>,Matrix_d, std::vector<double>);
  double getLogLikelihoodXOnly(AdmixOptions*,AlleleFreqs*,Chromosome**, Matrix_d, std::vector<double>);
  double IntegratingConst( double alpha, double beta, double a, double b );

  void SampleIndividualParameters( int i, Vector_d *SumLogTheta, AlleleFreqs *A, int iteration , MatrixArray_d *Target, 
				   Vector_i &OutcomeType, MatrixArray_d &ExpectedY, Vector_d &lambda, int NoCovariates, 
				   Matrix_d &Covariates0,MatrixArray_d &beta, Vector_d &poptheta, AdmixOptions* options, 
				   Chromosome **chrm, vector<Vector_d> alpha, bool _symmetric, vector<bool> _admixed, 
				   double rhoalpha, double rhobeta, vector<double> sigma);

 void OnePopulationUpdate( int i, MatrixArray_d *Target, Vector_i &OutcomeType, MatrixArray_d &ExpectedY, Vector_d &lambda, 
			   AlleleFreqs *Loci, int AnalysisTypeIndicator);

  void ChibLikelihood(int i,int iteration, double *LogLikelihood, double *SumLogLikelihood, vector<double> MaxLogLikelihood, 
		      AdmixOptions *options, Chromosome **chrm, vector<Vector_d> alpha,  
		      vector<bool> _admixed, double rhoalpha, double rhobeta, MatrixArray_d &thetahat, MatrixArray_d &thetahatX, 
		      vector<vector<double> > &rhohat,  vector<vector<double> > &rhohatX,
		      std::ofstream *LogFileStreamPtr, chib *MargLikelihood, AlleleFreqs *A);

private:
  static std::vector<unsigned int>
  encodeGenotype(std::vector<unsigned int>&); // doesn't appear to do anything - should add 1 to decoded genotypes 
  static void
  s2c(char *c, std::string s);
   
  std::vector< std::vector<unsigned int> > _genotype; // stores "encoded" genotypes: missing=0, alleles numbered from 1
  Vector_i *PossibleHaplotypes;
  std::vector< std::vector<unsigned int> > new_genotype; // also stores encoded genotypes    
  std::vector< std::vector<bool> > _xi;
  std::vector< unsigned int > numCompLoci;
  unsigned int numChromosomes;
  Matrix_d AdmixtureProps;
  Matrix_d XAdmixtureProps;
  Matrix_d AdmixtureHat;
  Matrix_d XAdmixtureHat;
  MatrixArray_i LocusAncestry;
  double Sumrho0;
  Vector_i sumxi;
  std::vector< double > _rho;
  std::vector< double > _rho_X;
  std::vector< double > _rhoHat;
  std::vector< double > _rhoHat_X;
  std::vector< Matrix_d > AncestryProbs; //Conditional probabilities of locus ancestry
  double LogPosterior;
  unsigned int sex;
  std::vector< unsigned int > gametes;
  unsigned int X_posn;
  double TruncationPt; // upper truncation point for sum intensities parameter rho

  void UpdateAdmixtureForRegression( int i,int Populations, int NoCovariates, Vector_d &poptheta, bool ModelIndicator,
				     Matrix_d *Covariates0);
  void Accept_Reject_Theta( double p, Matrix_d &theta, Matrix_d &thetaX, bool xdata, int Populations, bool ModelIndicator );
  double AcceptanceProbForTheta_XChrm(Matrix_d &Theta, Matrix_d &ThetaX,std::vector<double> &sigma, int Populations );
  double AcceptanceProbForTheta_LogReg( int i, int TI, Matrix_d &theta ,bool ModelIndicator,int Populations, 
					int NoCovariates, Matrix_d &Covariates0, MatrixArray_d &beta, MatrixArray_d &ExpectedY, 
					MatrixArray_d &Target, Vector_d &poptheta);
  double AcceptanceProbForTheta_LinearReg( int i, int TI,  Matrix_d &theta ,bool ModelIndicator,int Populations,
					   int NoCovariates, Matrix_d &Covariates0, MatrixArray_d &beta, MatrixArray_d &ExpectedY,
					   MatrixArray_d &Target, Vector_d &poptheta, Vector_d &lambda);

  void SampleLocusAncestry(Chromosome **chrm, AdmixOptions *options, AlleleFreqs *A, vector< Vector_d > f, 
			   Matrix_i *, Matrix_i *);

  void SampleNumberOfArrivals(AlleleFreqs *A, AdmixOptions *options, vector< unsigned int > *SumN, 
			      vector< unsigned int > *SumN_X);

  void SampleRho(bool XOnly, bool RandomMatingModel, bool X_data, double rhoalpha, double rhobeta, double L, double L_X, 
			   vector< unsigned int > SumN, vector< unsigned int > SumN_X);

  void CalculateLogPosterior(AdmixOptions *options, bool isX_data, vector<Vector_d> alpha, 
						 bool _symmetric, vector<bool> _admixed, double rhoalpha, double rhobeta, double L, 
						 double L_X, vector< unsigned int > SumN, vector< unsigned int > SumN_X, 
						 Matrix_i &SumLocusAncestry, Matrix_i &SumLocusAncestry_X);

  void InitializeChib(Matrix_d theta, Matrix_d thetaX, vector<double> rho, vector<double> rhoX, 
		 AdmixOptions *options, AlleleFreqs *A, Chromosome **chrm, double rhoalpha, double rhobeta, 
		 vector<Vector_d> alpha, vector<bool> _admixed, chib *MargLikelihood, std::ofstream *LogFileStreamPtr);

  void setIsMissing(vector<unsigned int >& decoded);
};

#endif /* INDIVIDUAL_H */
