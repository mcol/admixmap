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
  IndividualCollection(const AdmixOptions* const options, const InputData* const Data, const Genome& Loci, 
					   const Chromosome* const* chrm);

  void Initialise(const AdmixOptions* const options, const Genome* const Loci, const string* const PopulationLabels,
				      double rhoalpha, double rhobeta, LogWriter &Log, const DataMatrix &MLEMatrix);

  void LoadData(const AdmixOptions* const options, const InputData* const);
  
  void getOnePopOneIndLogLikelihood(LogWriter &Log, const std::string* const PopulationLabels);

  void Update(int iteration, const AdmixOptions* const options, Chromosome **chrm, AlleleFreqs *A,
	      const Regression* const R, const double* const poptheta,
	      const std::vector<std::vector<double> > &alpha, double globalrho, double rhoalpha, double rhobeta,
	      LogWriter &Log, double coolness, bool anneal);

  void ConjugateUpdateIndAdmixture(int iteration, const Regression* const R, const double* const poptheta, 
				   const AdmixOptions* const options, Chromosome **chrm, 
				   const vector<vector<double> > &alpha);
  void setGenotypeProbs(Chromosome** C, unsigned nchr);
  
  void OutputIndAdmixture();

  void OutputDeviance(const AdmixOptions* const options, Chromosome** C, Regression *R, LogWriter &Log, 
		      double SumRho, unsigned numChromosomes);

  void OutputChibEstimates(LogWriter &, int)const;
  void OutputChibResults(LogWriter&)const;

  void OutputErgodicAvg(int samples, bool ML, std::ofstream *avgstream);
  void OutputResiduals(const char* ResidualFilename, const Vector_s Labels, int iterations);

  int getSize()const;

  void add(Individual*);

  Individual* getIndividual(int)const;

  void setAdmixtureProps(const double* const, size_t);
  void setAdmixturePropsX(const double* const, size_t);

  double GetSumrho()const;
  double getSumLogTheta(int)const;
  const double *getSumLogTheta()const;
  std::vector<double> getOutcome(int)const;
  double getOutcome(int, int)const;
  int getNumberOfOutcomeVars()const;
  int GetNumCovariates() const;
  int GetNumberOfInputCovariates()const;

  const double* getCovariates()const;
  DataType getOutcomeType(int)const;

  void SetExpectedY(int, const double*const );
  //void calculateExpectedY(int);
  void UpdateSumResiduals();
  double getExpectedY(int)const;

  const std::string getCovariateLabels(int)const;
  const std::string *getCovariateLabels()const;

  double getLogLikelihood(const AdmixOptions* const options, Chromosome **C, const Regression* R, bool sumdeviance);
  double getModifiedLogLikelihood(const AdmixOptions* const options, Chromosome **C, double coolness);
  double DerivativeInverseLinkFunction(int i)const;
  void ResetChib();
  const chib* getChib()const;
private:
  Individual **_child; //array of pointers to Individual
  void getLabels(const Vector_s& data,  string *labels);

  void LoadCovariates(const InputData*, const AdmixOptions* const options);
  void LoadOutcomeVar(const InputData* const);
  void LoadRepAncestry(const InputData* const);
  void InitialiseMLEs(double, double, const AdmixOptions* const, const DataMatrix&);

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
  double SumLogLikelihood;
  double SumDeviance, SumDevianceSq;
  std::vector< int > _locusfortest;
 
  //Regression Objects
  double **ExpectedY;
  double **SumResiduals;
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
