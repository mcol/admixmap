// *-*-C++-*-*
/** 
 *   IndividualCollection.h
 *   header file for IndividualCollection class
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
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
#include "bclib/DataMatrix.h"
#include "Individual.h"
#include <vector>
#include <string.h>
#include <string>

class Regression;
class Individual;
class AlleleFreqs;
class AlleleFreqSampler;

namespace bclib{
  class LogWriter;
}
///Class to hold an array of Individuals
class IndividualCollection
{
public:
  virtual ~IndividualCollection();
  IndividualCollection(unsigned numIndividuals, unsigned numPopulations, unsigned numCompositeLoci);

  void DeleteGenotypes(bool);
  //virtual void Initialise(const AdmixOptions* const options, const Genome* const Loci,
		  //const Vector_s& PopulationLabels, bclib::LogWriter &Log) = 0;
  virtual void LoadData(const Options* const options, const InputData* const, bool admixtureAsCovariate);
  virtual void getOnePopOneIndLogLikelihood(bclib::LogWriter &, const Vector_s& ){};

  void SampleHapPairs(const Options& options, AlleleFreqs *A, const Genome* const Loci,
		      bool skipMissingGenotypes, bool anneal, bool UpdateCounts);

  void AccumulateAlleleCounts(const Options& options, AlleleFreqs *A, const Genome* const Loci,
                              bool anneal);

  int getSize()const;
  int getNumDiploidIndividuals()const;
  virtual int getNumberOfIndividualsForScoreTests()const{return getSize();}
  virtual unsigned int getFirstScoreTestIndividualNumber()const{return 0;};

  Individual* getIndividual(int)const;

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
  const bclib::DataMatrix& getCovariatesMatrix()const;
  const bclib::DataMatrix& getOutcomeMatrix()const;

  double getLogLikelihood(const Options& options, bool forceupdate);
  double getEnergy(const Options& options, const vector<bclib::Regression*> &R, 
			  const bool & annealed);

  virtual void HMMIsBad(bool b);
  virtual void resetStepSizeApproximators(int ){};

protected:
  const unsigned int NumInd;
  unsigned int size;
  unsigned NumDiploidIndividuals;
  const int Populations, NumCompLoci;
  virtual void SetNullValues();
  Individual** _child;//pointer to _child array

  //Regression Objects
  bclib::DataMatrix Outcome;
  int NumOutcomes;
  int NumCovariates;//covariates including admixture
  int NumberOfInputCovariates;//covariates in file
  bclib::DataMatrix Covariates;//all covariates, including admixture props
  DataType *OutcomeType;

  bclib::DataMatrix *ReportedAncestry;
  
  double SumLogLikelihood;
  double SumDeviance, SumDevianceSq;
  //std::vector< int > _locusfortest;

private:
  IndividualCollection();
  //double* SumEnergy, *SumEnergySq;//to store sum over iters of energy of test ind at each coolness
  void getLabels(const Vector_s& data, Vector_s& labels);
  
  void LoadCovariates(const InputData*, const Options* const options, bool admixtureAsCovariate);
  void LoadOutcomeVar(const InputData* const);

  std::ofstream EYStream;//output file for expected outcome

};

#endif /* !defined INDIVIDUAL_COLLECTION_H */


