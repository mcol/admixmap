// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   IndividualCollection.h 
 *   header file for IndividualCollection class
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
#ifndef INDIVIDUAL_COLLECTION_H
#define INDIVIDUAL_COLLECTION_H 1

#include "Genome.h"
#include "Chromosome.h"
#include "AlleleFreqs.h"
#include "chib.h"
#include "Individual.h"
#include "IndAdmixOutputter.h"
#include "matrix.h"
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
class LogWriter;

class IndividualCollection
{
public:
  IndividualCollection();
  ~IndividualCollection();
  IndividualCollection(AdmixOptions*,InputData *Data,Genome&,Chromosome **);

  void Initialise(AdmixOptions *, double **,Genome *,std::string *PopulationLabels, double rhoalpha,double rhobeta, LogWriter *Log,
		  const Matrix_d &MLEMatrix);

  void LoadData(AdmixOptions *options, InputData *, LogWriter *Log);
  
  void getOnePopOneIndLogLikelihood(LogWriter *Log, std::string *PopulationLabels);

  void Update(int iteration, AlleleFreqs *A, Regression *R, 
	      const double *poptheta, AdmixOptions *options, Chromosome **chrm, 
	      std::vector<std::vector<double> > &alpha, double rhoalpha, double rhobeta,
	      LogWriter *Log, chib *MargLikelihood);
  
  void OutputIndAdmixture();

  void OutputChibEstimates(LogWriter *, int);

  void OutputErgodicAvg(int samples, chib *MargLikelihood, std::ofstream *avgstream);

  int getSize();

  void add(Individual*);

  Individual* getIndividual(int);

  void setAdmixtureProps(double *, size_t);
  void setAdmixturePropsX(double *, size_t);

  int *GetSumXi();
  int GetSumXi(int j);
  double GetSumrho0();
  double GetSumrho();
  double getSumLogTheta(int);
  double *getSumLogTheta();
  Matrix_d getOutcome(int);
  double getOutcome(int, int);
  Vector_d getTargetCol(int,int);
  int getNumberOfOutcomeVars();
  int GetNumCovariates() const;
  int GetNumberOfInputCovariates();

  Matrix_d getCovariates();
  int getOutcomeType(int);
  Vector_i *getOutcomeType();

  void SetExpectedY(int,Matrix_d);
  void SetExpectedY(int,double *);
  void calculateExpectedY(int);
  double getExpectedY(int);

  std::string getTargetLabels(int);
  std::string getCovariateLabels(int);
  std::string *getCovariateLabels();

  double getLL();
  double DerivativeInverseLinkFunction(int AnalysisType,int i);
private:
  Individual **_child; //array of pointers to Individual
  void getLabels(const Vector_s& data,  string *labels);

  void LoadCovariates(InputData *);
  void LoadOutcomeVar(AdmixOptions *options, InputData *, LogWriter *Log);
  void LoadRepAncestry(InputData *);
  void InitialiseMLEs(double, double, AdmixOptions *, const Matrix_d&);

  unsigned int NumInd, NumCompLoci;
  //MLEs of Individual admixture and sumintensities
  //used to calculate marginal likelihood
  vector< vector<double> > rhohat, rhohatX;
  double **thetahat;
  double **thetahatX;
  double *SumLogTheta;//sums of log individual admixture proportions
  vector<double> MaxLogLikelihood;

  //Regression Objects
  double **ExpectedY;
  Matrix_d *Outcome;
  int NumOutcomes;
  int NumCovariates;
  Matrix_d Covariates;//all covariates, including admixture props
  Matrix_d InputCovariates;//covariates from file
  std::string *CovariateLabels;
  std::string *OutcomeVarLabels;
  int *OutcomeType;

  Matrix_d *ReportedAncestry;
  std::vector<double> sigma;
  IndAdmixOutputter* indadmixoutput;
  double LogLikelihood, SumLogLikelihood;
  std::vector< int > _locusfortest;
 


};

#endif /* !defined INDIVIDUAL_COLLECTION_H */
