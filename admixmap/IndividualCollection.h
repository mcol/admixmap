// *-*-C++-*-*
#ifndef INDIVIDUAL_COLLECTION_H
#define INDIVIDUAL_COLLECTION_H 1

#include "Genome.h"
#include "Chromosome.h"
#include "AlleleFreqs.h"
#include "chib.h"
#include "Individual.h"
#include "IndividualVisitor.h"
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

class IndividualCollection
{
private:
  std::vector<Individual*> _child;
  void getLabels( const string, Vector_i temporary, string *labels );
  void getLabels(const Vector_s& data, Vector_i temporary, string *labels);

  void LoadCovariates(AdmixOptions *options, InputData *, LogWriter *Log);
  void LoadOutcomeVar(AdmixOptions *options, InputData *, LogWriter *Log);
  void LoadRepAncestry(AdmixOptions *options, InputData *, LogWriter *Log);
  void CheckGenotypes(Genome *Loci, LogWriter *Log);

  vector< vector<double> > rhohat, rhohatX;
  MatrixArray_d thetahat;
  MatrixArray_d thetahatX;
  vector<double> MaxLogLikelihood;
  MatrixArray_d ExpectedY;
  MatrixArray_d Target;
  MatrixArray_d Covariates;
  MatrixArray_d ReportedAncestry;
  Matrix_d Input;
  std::string *CovariateLabels;
  std::string *TargetLabels;
  Vector_i OutcomeType;
  std::vector<double> sigma;
  IndividualVisitor* indadmixoutput;
  double LogLikelihood, SumLogLikelihood;
  std::vector< int > _locusfortest;

public:
  IndividualCollection();
  ~IndividualCollection();
  IndividualCollection(AdmixOptions*,const Matrix_s& data,Genome&,Chromosome **);

  void Initialise(AdmixOptions *, MatrixArray_d *,Genome *,std::string *PopulationLabels);

  void PreUpdate(double, double, AdmixOptions *);

  void LoadGenotypes(AdmixOptions *options, InputData *, LogWriter *Log, Genome *Loci);
  
  void getOnePopOneIndLogLikelihood(LogWriter *Log, AlleleFreqs *A, std::string *PopulationLabels);

  void Update(int iteration, Vector_d *SumLogTheta, AlleleFreqs *A, Vector_d *lambda, int NoCovariates, MatrixArray_d *beta, 
	      Vector_d &poptheta, AdmixOptions *options, Chromosome **chrm, 
	      vector<Vector_d> alpha, bool _symmetric, vector<bool> _admixed, double rhoalpha, double rhobeta,
	      std::ofstream *LogFileStreamPtr, chib *MargLikelihood);
  
  void accept(Matrix_d);

  void Output(std::ofstream *);

  void OutputErgodicAvg(int samples, chib *MargLikelihood, std::ofstream *avgstream);

  int getSize();

  void add(Individual*);

  Individual* getIndividual(int);

  MatrixArray_d getAncestries();
  void setAncestry(Matrix_d);
  void setAncestryX(Matrix_d);

  Vector_i GetSumXi();
  double GetSumrho0();
  double GetSumrho();
  MatrixArray_d getOutcome();
  Matrix_d getOutcome(int);
  Vector_d getTargetCol(int,int);
  int getTargetSize();
  int GetNumberOfInputRows();
  int GetNumberOfInputCols();

  MatrixArray_d getCovariates();
  Matrix_d getCovariates(int);
  int getOutcomeType(int);
  Vector_i *getOutcomeType();

  void SetExpectedY(int,Matrix_d);
  void calculateExpectedY(int);
  double getExpectedY(int);

  std::string getTargetLabels(int);
  std::string getCovariateLabels(int);
  std::string *getCovariateLabels();

  double getLL();
  double DerivativeInverseLinkFunction(int AnalysisType, int i);
};

#endif /* !defined INDIVIDUAL_COLLECTION_H */
