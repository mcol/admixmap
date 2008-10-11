// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   AdmixIndividualCollection.h
 *   header file for Admixed IndividualCollection class
 *   Copyright (c) 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef ADMIX_INDIVIDUAL_COLLECTION_H
#define ADMIX_INDIVIDUAL_COLLECTION_H 1

#include "IndividualCollection.h"
#include "chib.h"
#include "AdmixedIndividual.h"
#include "IndAdmixOutputter.h"
class InputAdmixData;
class IndAdmixOutputter;

///Class to hold an array of AdmixedIndividuals
class AdmixIndividualCollection : public IndividualCollection
{
public:
  ~AdmixIndividualCollection();
  AdmixIndividualCollection(const AdmixOptions* const options, const InputAdmixData* const Data, Genome* Loci);

  //void DeleteGenotypes(bool);
  void Initialise(const AdmixOptions& options, const Genome& Loci,
		  const Vector_s& PopulationLabels, bclib::LogWriter &Log);
  void DrawInitialAdmixture(const std::vector<std::vector<double> > &alpha);
  void LoadData(const AdmixOptions* const options, const InputAdmixData* const);
  void getOnePopOneIndLogLikelihood(bclib::LogWriter &Log, const Vector_s& PopulationLabels);

  void HMMUpdates(int iteration, const AdmixOptions& options,
                  const vector<bclib::Regression*> &R, const double* const poptheta,
                  const vector<vector<double> > &alpha, 
                  AffectedsOnlyTest& affectedsOnlyTest, CopyNumberAssocTest& ancestryAssocTest, bool anneal);
  void SampleAdmixtureWithRandomWalk(int iteration, const AdmixOptions* const options,
				     const vector<bclib::Regression*> &R, const double* const poptheta,
				     const vector<vector<double> > &alpha, CopyNumberAssocTest& ancestryAssocTest, bool anneal);

  void SampleParameters(int iteration, const AdmixOptions& options,
			const vector<bclib::Regression*> &R, const double* const poptheta,
			const vector<vector<double> > &alpha, double rhoalpha, double rhobeta,
			CopyNumberAssocTest& ancestryAssocTest, bool anneal);
  void setChibNumerator(const AdmixOptions& options, const vector<vector<double> > &alpha, 
		  double rhoalpha, double rhobeta, AlleleFreqs *A);
  void updateChib(const AdmixOptions& options,const vector<vector<double> > &alpha, 
		  double rhoalpha, double rhobeta, AlleleFreqs *A);

  void FindPosteriorModes(const AdmixOptions& options, 
			  const vector<bclib::Regression*> &R, 
			  const vector<vector<double> > &alpha, double rhoalpha, double rhobeta, AlleleFreqs* A, 
			  const Vector_s& PopulationLabels);

  void OutputIndAdmixture();
  double getDevianceAtPosteriorMean(const Options& options, vector<bclib::Regression *> &R, Genome* Loci,
				    bclib::LogWriter &Log, const vector<double>& SumLogRho, unsigned numChromosomes,  AlleleFreqs* A);

  void WritePosteriorMeans(const AdmixOptions& options, const vector<string>& PopLabels, 
			   Genome* Loci)const;
  void OutputChibResults(bclib::LogWriter&)const;

  AdmixedIndividual* getIndividual(int)const;
  //void setAdmixtureProps(const double* const, size_t);
  double GetSumrho()const;
  double getSumLogTheta(int)const;
  const double *getSumLogTheta()const;

//  double getLogLikelihood(const AdmixOptions* const options, bool forceupdate);
  //double getEnergy(const AdmixOptions* const options, const vector<bclib::Regression*> &R,
			 // const bool & annealed);
  void setGenotypeProbs(const Genome * const G);

  void annealGenotypeProbs(unsigned nchr, const double coolness, const double* Coolnesses);

  void HMMIsBad(bool b);
  void resetStepSizeApproximators(int k);
  void accumulateEnergyArrays(const Options& options);
  double* getSumEnergy()const;
  double* getSumEnergySq()const;
  void ResetChib();
  void OutputErgodicChib(std::ofstream *avgstream, bool fixedfreqs);
  const chib* getChib()const;

private:
  AdmixedIndividual** AdmixedChild;
  AdmixedIndividual** TestInd;// pointer to individual for whom to estimate marginal likelihood
  int sizeTestInd;

  double* SumEnergy, *SumEnergySq;//to store sum over iters of energy of test ind at each coolness

  std::vector< std::vector<double> > admixtureprior;
  double *SumLogTheta;//sums of log individual admixture proportions

  IndAdmixOutputter* indadmixoutput;
  double SumLogLikelihood;
  double SumDeviance, SumDevianceSq;
  std::vector< int > _locusfortest;

  chib MargLikelihood;

  AdmixIndividualCollection();
  void SetNullValues();
  void LoadRepAncestry(const InputAdmixData* const);
};

#endif /* !defined ADMIX_INDIVIDUAL_COLLECTION_H */


