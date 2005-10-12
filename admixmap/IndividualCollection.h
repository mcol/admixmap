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
#include "DataMatrix.h"

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

  void Initialise(AdmixOptions *, Genome *, std::string *PopulationLabels, 
		  double rhoalpha, double rhobeta, LogWriter *Log, const DataMatrix &MLEMatrix);

  void LoadData(AdmixOptions *options, InputData *);
  
  void getOnePopOneIndLogLikelihood(LogWriter *Log, std::string *PopulationLabels);

  void Update(int iteration, AdmixOptions *options, Chromosome **chrm, AlleleFreqs *A,
          Regression *R0, Regression *R1, const double *poptheta,
	      std::vector<std::vector<double> > &alpha, double globalrho, double rhoalpha, double rhobeta,
	      LogWriter *Log);

  void ConjugateUpdateIndAdmixture(int iteration, Regression *R0, Regression *R1, const double *poptheta, 
				   AdmixOptions *options, Chromosome **chrm, vector<vector<double> > &alpha);
  
  void OutputIndAdmixture();

  void OutputDeviance(AdmixOptions *options, Chromosome** C, LogWriter *Log, double SumRho, unsigned numChromosomes);

  void OutputChibEstimates(LogWriter *, int);

  void OutputErgodicAvg(int samples, std::ofstream *avgstream);

  int getSize();

  void add(Individual*);

  Individual* getIndividual(int);

  void setAdmixtureProps(double *, size_t);
  void setAdmixturePropsX(double *, size_t);

  double GetSumrho();
  double getSumLogTheta(int);
  double *getSumLogTheta();
  std::vector<double> getOutcome(int);
  double getOutcome(int, int);
  int getNumberOfOutcomeVars();
  int GetNumCovariates() const;
  int GetNumberOfInputCovariates();

  const double* getCovariates()const;
  DataType getOutcomeType(int);

  void SetExpectedY(int, const double*const );
  void calculateExpectedY(int);
  double getExpectedY(int);

  std::string getTargetLabels(int);
  std::string getCovariateLabels(int);
  std::string *getCovariateLabels();

  double getLL();
  double DerivativeInverseLinkFunction(int i);
private:
  Individual **_child; //array of pointers to Individual
  void getLabels(const Vector_s& data,  string *labels);

  void LoadCovariates(InputData *, AdmixOptions* options);
  void LoadOutcomeVar(InputData *);
  void LoadRepAncestry(InputData *);
  void InitialiseMLEs(double, double, AdmixOptions *, const DataMatrix&);

  unsigned int NumInd, NumCompLoci;
  //MLEs of Individual admixture and sumintensities
  //used to calculate marginal likelihood
  vector<double> rhohat, rhohatX;
  double *thetahat;
  double *thetahatX;

  double *SumLogTheta;//sums of log individual admixture proportions
  vector<double> MaxLogLikelihood;

  DataMatrix *ReportedAncestry;
  std::vector<double> sigma;
  IndAdmixOutputter* indadmixoutput;
  double LogLikelihood, SumLogLikelihood;
  double SumDeviance, SumDevianceSq;
  std::vector< int > _locusfortest;
 
  //Regression Objects
  double **ExpectedY;
  DataMatrix Outcome;
  int NumOutcomes;
  int NumCovariates;//covariates including admixture
  int NumberOfInputCovariates;//covariates in file
  DataMatrix Covariates;//all covariates, including admixture props
  std::string *CovariateLabels;
  DataType *OutcomeType;

  chib MargLikelihood;
};

#endif /* !defined INDIVIDUAL_COLLECTION_H */
