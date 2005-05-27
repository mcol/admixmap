// *-*-C++-*-*
#ifndef INDIVIDUAL_COLLECTION_H
#define INDIVIDUAL_COLLECTION_H 1

#include "Genome.h"
#include "Chromosome.h"
#include "AlleleFreqs.h"
#include "chib.h"
#include "Individual.h"
#include "IndAdmixOutputter.h"
#include "matrix.h"
#include "MatrixArray_d.h"
#include "matrix_d.h"
#include "matrix_i.h"
#include "vector.h"
#include "vector_i.h"
#include <vector>
#include <string.h>
#include <string>

class Regression;
class Individual;
class IndAdmixOutputter;

class IndividualCollection
{
private:
  Individual **_child;
  void getLabels( const string, Vector_i temporary, string *labels );
  void getLabels(const Vector_s& data, Vector_i temporary, string *labels);

  void LoadCovariates(AdmixOptions *options, InputData *, LogWriter *Log);
  void LoadOutcomeVar(AdmixOptions *options, InputData *, LogWriter *Log);
  void LoadRepAncestry(AdmixOptions *options, InputData *, LogWriter *Log);
  void CheckGenotypes(Genome *Loci, LogWriter *Log);
  void InitialiseMLEs(double, double, AdmixOptions *, const Matrix_d&);

  unsigned int NumInd, NumCompLoci;
  //MLEs of Individual admixture and sumintensities
  //used to calculate marginal likelihood
  vector< vector<double> > rhohat, rhohatX;
  MatrixArray_d thetahat;
  MatrixArray_d thetahatX;
  vector<double> MaxLogLikelihood;

  //Regression Objects
  MatrixArray_d ExpectedY;
  MatrixArray_d Target;
  Matrix_d Covariates;
  Matrix_d Input;
  std::string *CovariateLabels;
  std::string *TargetLabels;
  Vector_i OutcomeType;

  MatrixArray_d ReportedAncestry;
  std::vector<double> sigma;
  IndAdmixOutputter* indadmixoutput;
  double LogLikelihood, SumLogLikelihood;
  std::vector< int > _locusfortest;
  Vector_d SumLogTheta;

public:
  IndividualCollection();
  ~IndividualCollection();
  IndividualCollection(AdmixOptions*,InputData *Data,Genome&,Chromosome **);

  void Initialise(AdmixOptions *, MatrixArray_d *,Genome *,std::string *PopulationLabels, double rhoalpha,double rhobeta, LogWriter *Log,
		  const Matrix_d &MLEMatrix);

  void LoadGenotypes(AdmixOptions *options, InputData *, LogWriter *Log, Genome *Loci);
  
  void getOnePopOneIndLogLikelihood(LogWriter *Log, AlleleFreqs *A, std::string *PopulationLabels);

  void Update(int iteration, AlleleFreqs *A, Regression *R, 
	      Vector_d &poptheta, AdmixOptions *options, Chromosome **chrm, 
	      vector<Vector_d> alpha, bool _symmetric, vector<bool> _admixed, double rhoalpha, double rhobeta,
	      std::ofstream *LogFileStreamPtr, chib *MargLikelihood);
  
  void OutputIndAdmixture();

  void Output(std::ofstream *);

  void OutputErgodicAvg(int samples, chib *MargLikelihood, std::ofstream *avgstream);

  int getSize();

  void add(Individual*);

  Individual* getIndividual(int);

  MatrixArray_d getAncestries();
  void setAdmixtureProps(Matrix_d);
  void setAdmixturePropsX(Matrix_d);

  int *GetSumXi();
  int GetSumXi(int j);
  double GetSumrho0();
  double GetSumrho();
  double getSumLogTheta(int);
  MatrixArray_d getOutcome();
  Matrix_d getOutcome(int);
  Vector_d getTargetCol(int,int);
  int getTargetSize();
  int GetNumberOfInputRows();
  int GetNumberOfInputCols();

  Matrix_d getCovariates();
  int getOutcomeType(int);
  Vector_i *getOutcomeType();

  void SetExpectedY(int,Matrix_d);
  void calculateExpectedY(int);
  double getExpectedY(int);

  std::string getTargetLabels(int);
  std::string getCovariateLabels(int);
  std::string *getCovariateLabels();

  double getLL();
  double DerivativeInverseLinkFunction(int AnalysisType,int i);
};

#endif /* !defined INDIVIDUAL_COLLECTION_H */
