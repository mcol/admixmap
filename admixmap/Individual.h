// *-*-C++-*-*
#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H 1
#include "common.h"
#include "IndividualVisitor.h"
#include "Genome.h"
#include "chib.h"
#include <gsl/gsl_cdf.h>

using namespace::std;

class Individual
{
public:
  Individual();

  Individual(AdmixOptions*,const Vector_s& data,Genome&,Genome&);

  ~Individual();

  void
  accept(IndividualVisitor&, double, std::vector<int>, double);

  int getSex();

  Matrix_d&
  getAncestry();

  void
  setAncestry(Matrix_d);
  
  Matrix_d&
  getAncestryX();

  void
  setAncestryX(Matrix_d);
  
  //   std::vector< std::vector<unsigned int> >&
  //   getGenotype();

  std::vector<unsigned int>&
  getGenotype(unsigned int locus);

  std::vector< std::vector<unsigned int> >&
  IsMissing();

  std::vector<unsigned int>&
  IsMissing(unsigned int locus);

  void
  setGenotype(unsigned int locus,std::vector<unsigned int> genotype);

  std::vector<bool>&
  getXi(unsigned int locus);

  const std::vector< std::vector<bool> >&
  getXi();

  Vector_i getSumXi();

  double getSumrho0();

  double getSumrho();

  std::vector<double> getRho();

  Matrix_d getExpectedAncestry( int );

  double getLogLikelihood(AdmixOptions*,Genome&,Genome&);

  double getLogLikelihoodOnePop(Genome&);

  double getLogPosteriorProb();

  Vector_i GetLocusAncestry( int, int );
   
  Matrix_d SampleParameters(int, AdmixOptions*, Vector_d&, Genome&,
			    Genome&, std::vector<Vector_d>,
			    bool, std::vector<bool>, double,
			    double, int, std::vector<double>, Matrix_d& );
  double getLogLikelihood(AdmixOptions*,Genome&,Genome&,
			  Matrix_d, std::vector<double>,Matrix_d, std::vector<double>);
  double getLogLikelihoodXOnly(AdmixOptions*,Genome&,Genome&,
			       Matrix_d, std::vector<double>);
  double IntegratingConst( double alpha, double beta, double a, double b );

  void IndivUpdate(int i,int iteration,  
		   Vector_d *SumLogTheta, MatrixArray_d *Target, Vector_i &TargetType, MatrixArray_d &ExpectedY, 
		   Vector_d &lambda,
		   int NoCovariates, Matrix_d &Covariates0, MatrixArray_d &beta, Vector_d &poptheta, AdmixOptions *options,
		   Vector_d &f, Genome *Loci, Genome *chrm, vector<Vector_d> alpha, bool _symmetric, 
		   vector<bool> _admixed, double rhoalpha, double rhobeta, vector<double> sigma);

  void ChibLikelihood(int i,int iteration, double *LogLikelihood, double *SumLogLikelihood, vector<double> MaxLogLikelihood, 
		      AdmixOptions *options, Genome *Loci, Genome *chrm, vector<Vector_d> alpha,  
		      vector<bool> _admixed, double rhoalpha, double rhobeta, MatrixArray_d &thetahat, MatrixArray_d &thetahatX, 
		      vector<vector<double> > &rhohat,  vector<vector<double> > &rhohatX,
		      std::ofstream *LogFileStreamPtr, chib *MargLikelihood);

private:
  static std::vector<unsigned int>
  encodeGenotype(std::vector<unsigned int>&);
  static void
  s2c(char *c, std::string s);
   
  std::vector< std::vector<unsigned int> > _genotype;
  std::vector< std::vector<unsigned int> > new_genotype;
  std::vector< std::vector<bool> > _xi;
  std::vector< unsigned int > numCompLoci;
  unsigned int numChromosomes;
  Matrix_d _ancestry;
  Matrix_d _ancestryX;
  Matrix_d _ancestryHat;
  Matrix_d _ancestryHat_X;
  MatrixArray_i LocusAncestry;
  double Sumrho0;
  Vector_i sumxi;
  std::vector< double > _rho;
  std::vector< double > _rho_X;
  std::vector< double > _rhoHat;
  std::vector< double > _rhoHat_X;
  std::vector< Matrix_d > ExpectedAncestry;
  double LogPosterior;
  unsigned int sex;
  std::vector< unsigned int > gametes;
  unsigned int X_posn;
  double TruncationPt;


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
  void SampleIndividualParameters( int i, Vector_d *SumLogTheta, int iteration , MatrixArray_d *Target, Vector_i &TargetType, 
				   MatrixArray_d &ExpectedY, 
				   Vector_d &lambda, int NoCovariates, Matrix_d &Covariates0,MatrixArray_d &beta, 
				   Vector_d &poptheta, AdmixOptions* options, Vector_d &f, 
				   Genome &Loci, Genome &chrm, vector<Vector_d> alpha, bool _symmetric, 
				   vector<bool> _admixed, double rhoalpha, double rhobeta, vector<double> sigma);

  void OnePopulationUpdate( int i, MatrixArray_d *Target, Vector_i &TargetType, MatrixArray_d &ExpectedY, Vector_d &lambda, 
			    Genome *Loci, int AnalysisTypeIndicator);

  void
  InitializeChib(Matrix_d theta, Matrix_d thetaX, vector<double> rho, vector<double> rhoX, 
		 AdmixOptions *options, Genome *Loci, Genome *chrm, double rhoalpha, double rhobeta, 
		 vector<Vector_d> alpha, vector<bool> _admixed, chib *MargLikelihood, std::ofstream *LogFileStreamPtr );

};

#endif /* INDIVIDUAL_H */
