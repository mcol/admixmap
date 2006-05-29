// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   IndividualCollection.h 
 *   header file for IndividualCollection class
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef INDIVIDUAL_COLLECTION_H
#define INDIVIDUAL_COLLECTION_H 1

#include "Genome.h"
//#include "Chromosome.h"
//#include "AlleleFreqs.h"
#include "chib.h"
#include "Individual.h"
#include "IndAdmixOutputter.h"
#include "DataMatrix.h"
#include "admixmap.h"//for PARALLEL

#include <vector>
#include <string.h>
#include <string>

class Regression;
class Individual;
class IndAdmixOutputter;
class LogWriter;
class AlleleFreqs;
class AlleleFreqSampler;

class IndividualCollection
{
public:
  IndividualCollection();
  ~IndividualCollection();
  IndividualCollection(const AdmixOptions* const options, const InputData* const Data, Genome* Loci);
  void SetNullValues();
  void DeleteGenotypes(bool);

  void Initialise(const AdmixOptions* const options, const Genome* const Loci, const string* const PopulationLabels,
		  const std::vector<std::vector<double> > &alpha, double rhoalpha, double rhobeta, 
		  LogWriter &Log);

  void LoadData(const AdmixOptions* const options, const InputData* const);
  
  void getOnePopOneIndLogLikelihood(LogWriter &Log, const std::string* const PopulationLabels);

  void SampleLocusAncestry(int iteration, const AdmixOptions* const options,
			   const Regression* const R, const double* const poptheta,
			   const vector<vector<double> > &alpha, 
			   bool anneal);
  void SampleHapPairs(const AdmixOptions* const options, AlleleFreqs *A, const Genome* const Loci,
		      bool anneal);
  void SampleParameters(int iteration, const AdmixOptions* const options,
			const Regression* const R, const double* const poptheta,
			const vector<vector<double> > &alpha, double rhoalpha, double rhobeta,
			bool anneal);
  void UpdateChib(int iteration, const AdmixOptions* const options,const vector<vector<double> > &alpha, 
		  double rhoalpha, double rhobeta, AlleleFreqs *A);

  void FindPosteriorModes(const AdmixOptions* const options, 
			  const Regression* const R, 
			  const vector<vector<double> > &alpha, double rhoalpha, double rhobeta, 
			  const std::string* const PopulationLabels);

  void setGenotypeProbs(const Genome * const G, const AlleleFreqs* const A);
  void annealGenotypeProbs(unsigned nchr, const double coolness, const double* Coolnesses);
  
  void OutputIndAdmixture();

  double getDevianceAtPosteriorMean(const AdmixOptions* const options, Regression *R, Genome* Loci, LogWriter &Log, 
				    const vector<double>& SumRho, unsigned numChromosomes, AlleleFreqs* A);

  void OutputChibEstimates(bool, LogWriter &, int)const;
  void OutputChibResults(LogWriter&)const;

  //void OutputResiduals(const char* ResidualFilename, const Vector_s Labels, int iterations);
  void OpenExpectedYFile(const char* Filename, LogWriter & Log);
  void OutputExpectedY(int k);
  void FinishWritingEYAsRObject(unsigned NumIterations, const Vector_s Labels);

  int getSize()const;
  void add(Individual*);

  Individual* getIndividual(int)const;

  void setAdmixtureProps(const double* const, size_t);
  //void setAdmixturePropsX(const double* const, size_t);

  double GetSumrho()const;
  double getSumLogTheta(int)const;
  const double *getSumLogTheta()const;
  const int* getSumAncestry()const;
  const vector<int> getSumLocusAncestry(int k)const;
  const vector<int> getSumLocusAncestryX(int k)const;
  const vector<unsigned> getSumNumArrivals();
  unsigned GetSNPAlleleCounts(unsigned locus, int allele)const;
  int getNumberOfMissingGenotypes(unsigned locus)const;
  const vector<int> getAlleleCounts(unsigned locus, int pop, unsigned NumStates)const;

  std::vector<double> getOutcome(int)const;
  double getOutcome(int, int)const;
    bool isMissingOutcome(int, int)const;
  int getNumberOfOutcomeVars()const;
  int GetNumCovariates() const;
  int GetNumberOfInputCovariates()const;
  const double* getCovariates()const;
  DataType getOutcomeType(int)const;
  void SetExpectedY(int, const double*const );
  //void UpdateSumResiduals();
  double getExpectedY(int)const;
  double getExpectedY(int, int)const;
  const std::string getCovariateLabels(int)const;
  const Vector_s getCovariateLabels()const;
  double getSampleVarianceOfOutcome(int j)const;
  double getSampleVarianceOfCovariate(int j)const;

  double getEnergy(const AdmixOptions* const options, const Regression* R, 
			  const bool & annealed);
  void accumulateEnergyArrays(const AdmixOptions* const options);
  double* getSumEnergy();
  double* getSumEnergySq();

  double DerivativeInverseLinkFunction(int i)const;
  void ResetChib();
  void OutputErgodicChib(std::ofstream *avgstream);
  const chib* getChib()const;

  void HMMIsBad(bool b);
  void resetStepSizeApproximators(int k); 

private:
  Individual **_child; //array of pointers to Individual objects
  Individual** TestInd;// pointer to individual for whom to estimate marginal likelihood
  int sizeTestInd;
  int worker_rank, NumWorkers;
  double* SumEnergy, *SumEnergySq;//to store sum over iters of energy of test ind at each coolness
  void getLabels(const Vector_s& data, Vector_s& labels);
  
  void LoadCovariates(const InputData*, const AdmixOptions* const options);
  void LoadOutcomeVar(const InputData* const);
  void LoadRepAncestry(const InputData* const);
  void InitialiseMLEs(double, double, const AdmixOptions* const);

  unsigned int NumInd, size, NumCompLoci;
  //MLEs of Individual admixture and sumintensities
  //used to calculate marginal likelihood
  vector<double> rhohat, rhohatX;
  double *thetahat;
  double *thetahatX;
  std::vector< std::vector<double> > admixtureprior;

  double *SumLogTheta;//sums of log individual admixture proportions
  //vector<double> MaxLogLikelihood;

  DataMatrix *ReportedAncestry;
  //std::vector<double> sigma;
  IndAdmixOutputter* indadmixoutput;
  double SumLogLikelihood;
  double SumDeviance, SumDevianceSq;
  std::vector< int > _locusfortest;
  int* SumAncestry;//for update of locus-specific sumintensities
#ifdef PARALLEL
  int* GlobalSumAncestry;//SumAncestry summed over processes, kept on master processes
  int rank_with_freqs;
  //communicators for workers (Individual updaters),  workers+freqsampler and workers+master
  MPI::Intracomm workers, workers_and_freqs, workers_and_master;
  unsigned Populations;
#endif
 
  //Regression Objects
  double **ExpectedY;
  //double **SumResiduals;
  DataMatrix Outcome;
  int NumOutcomes;
  int NumCovariates;//covariates including admixture
  int NumberOfInputCovariates;//covariates in file
  DataMatrix Covariates;//all covariates, including admixture props
  Vector_s CovariateLabels;
  DataType *OutcomeType;
  std::ofstream EYStream;//output file for expected outcome

  chib MargLikelihood;
};

#endif /* !defined INDIVIDUAL_COLLECTION_H */


