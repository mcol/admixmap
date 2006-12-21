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
#include "AdmixedIndividual.h"
#include "IndAdmixOutputter.h"
#include "utils/DataMatrix.h"

#include <vector>
#include <string.h>
#include <string>

class Regression;
class Individual;
class IndAdmixOutputter;
class LogWriter;
class AlleleFreqs;
class AlleleFreqSampler;

///Class to hold an array of Individuals
class IndividualCollection
{
public:
  IndividualCollection();
  virtual ~IndividualCollection();
  IndividualCollection(const AdmixOptions* const options, const InputData* const Data, Genome* Loci);

  void DeleteGenotypes(bool);
  virtual void Initialise(const AdmixOptions* const options, const Genome* const Loci, 
		  const Vector_s& PopulationLabels, LogWriter &Log);
  void DrawInitialAdmixture(const std::vector<std::vector<double> > &alpha);
  void LoadData(const AdmixOptions* const options, const InputData* const, bool admixtureAsCovariate);
  void getOnePopOneIndLogLikelihood(LogWriter &Log, const Vector_s& PopulationLabels);

  void SampleLocusAncestry(int iteration, const AdmixOptions* const options,
			   const vector<Regression*> &R,  
			   AffectedsOnlyTest& affectedsOnlyTest, AncestryAssocTest& ancestryAssocTest, bool anneal);
  void SampleAdmixtureWithRandomWalk(int iteration, const AdmixOptions* const options,
				     const vector<Regression*> &R, const double* const poptheta,
				     const vector<vector<double> > &alpha, AncestryAssocTest& ancestryAssocTest, bool anneal);

//   void SampleLocusAncestry(int iteration, const AdmixOptions* const options,
// 			   const vector<Regression*> &R);
  
  void SampleHapPairs(const AdmixOptions* const options, AlleleFreqs *A, const Genome* const Loci,
		      bool skipMissingGenotypes, bool anneal);
  void SampleParameters(int iteration, const AdmixOptions* const options,
			const vector<Regression*> &R, const double* const poptheta,
			const vector<vector<double> > &alpha, double rhoalpha, double rhobeta,
			AncestryAssocTest& ancestryAssocTest, bool anneal);
  void setChibNumerator(const AdmixOptions* const options,const vector<vector<double> > &alpha, 
		  double rhoalpha, double rhobeta, AlleleFreqs *A);
  void updateChib(const AdmixOptions* const options,const vector<vector<double> > &alpha, 
		  double rhoalpha, double rhobeta, AlleleFreqs *A);

  void FindPosteriorModes(const AdmixOptions* const options, 
			  const vector<Regression*> &R, 
			  const vector<vector<double> > &alpha, double rhoalpha, double rhobeta, AlleleFreqs* A, 
			  const Vector_s& PopulationLabels);

  void setGenotypeProbs(const Genome * const G, const AlleleFreqs* const 
#ifdef PARALLEL
			A
#endif
			);

  void annealGenotypeProbs(unsigned nchr, const double coolness, const double* Coolnesses);
  void OutputIndAdmixture();
  double getDevianceAtPosteriorMean(const AdmixOptions* const options, vector<Regression *>&R, Genome* Loci, LogWriter &Log, 
				    const vector<double>& SumRho, unsigned numChromosomes, AlleleFreqs* A);
  void OutputChibResults(LogWriter&)const;
  int getSize()const;
  int getNumDiploidIndividuals();
  virtual int getNumberOfIndividualsForScoreTests()const{return getSize();}
  virtual int getFirstScoreTestIndividualNumber()const{return 0;};

  Individual* getIndividual(int)const;
  //void setAdmixtureProps(const double* const, size_t);
  double GetSumrho()const;
  double getSumLogTheta(int)const;
  const double *getSumLogTheta()const;

  const vector<int> getSumLocusAncestry(int k)const;
  const vector<int> getSumLocusAncestryX(int k)const;
  const vector<unsigned> getSumNumArrivals();
  unsigned GetSNPAlleleCounts(unsigned locus, int allele)const;
  int getNumberOfMissingGenotypes(unsigned locus)const;
  const vector<int> getAlleleCounts(unsigned locus, int pop, unsigned NumStates)const;
  const std::vector<double> getOutcome(int)const;
  double getOutcome(int, int)const;
  bool isMissingOutcome(int, int)const;
  int getNumberOfOutcomeVars()const;
  int GetNumCovariates() const;
  int GetNumberOfInputCovariates()const;
  DataType getOutcomeType(int)const;
  const DataMatrix& getCovariatesMatrix()const;
  const DataMatrix& getOutcomeMatrix()const;

  double getLogLikelihood(const AdmixOptions* const options, bool forceupdate);
  double getEnergy(const AdmixOptions* const options, const vector<Regression*> &R, 
			  const bool & annealed);
  void accumulateEnergyArrays(const AdmixOptions* const options);
  double* getSumEnergy();
  double* getSumEnergySq();
  void ResetChib();
  void OutputErgodicChib(std::ofstream *avgstream, bool fixedfreqs);
  const chib* getChib()const;
  void HMMIsBad(bool b);
  void resetStepSizeApproximators(int k); 

  //functions specific to hapmixmodel
  virtual void SampleLocusAncestry(const AdmixOptions* const ){};
  virtual const int* getSumAncestry()const{return 0;};
  virtual void AccumulateConditionalGenotypeProbs(const AdmixOptions* const, const Genome& ){};
  virtual void OutputCGProbs(const char* ){};

protected:
  unsigned int NumInd, size;
  int Populations, NumCompLoci;
  int worker_rank, NumWorkers;
  void SetNullValues();
  Individual** _child;//pointer to _child array

private:
  AdmixedIndividual** pchild;
  AdmixedIndividual** TestInd;// pointer to individual for whom to estimate marginal likelihood
  int sizeTestInd;

  double* SumEnergy, *SumEnergySq;//to store sum over iters of energy of test ind at each coolness
  void getLabels(const Vector_s& data, Vector_s& labels);
  
  void LoadCovariates(const InputData*, const AdmixOptions* const options, bool admixtureAsCovariate);
  void LoadOutcomeVar(const InputData* const);
  void LoadRepAncestry(const InputData* const);

  std::vector< std::vector<double> > admixtureprior;
  double *SumLogTheta;//sums of log individual admixture proportions
  DataMatrix *ReportedAncestry;
  IndAdmixOutputter* indadmixoutput;
  double SumLogLikelihood;
  double SumDeviance, SumDevianceSq;
  std::vector< int > _locusfortest;

#ifdef PARALLEL
  // int rank_with_freqs;
  //communicators for workers (Individual updaters)
  MPI::Intracomm workers;
#endif
 
  //Regression Objects
  DataMatrix Outcome;
  int NumOutcomes;
  int NumCovariates;//covariates including admixture
  int NumberOfInputCovariates;//covariates in file
  DataMatrix Covariates;//all covariates, including admixture props
  DataType *OutcomeType;
  std::ofstream EYStream;//output file for expected outcome
  chib MargLikelihood;
};

#endif /* !defined INDIVIDUAL_COLLECTION_H */


