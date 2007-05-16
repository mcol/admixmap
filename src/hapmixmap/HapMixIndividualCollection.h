// *-*-C++-*-*
/** 
 *   HAPMIXMAP
 *   HapMixIndividualCollection.h 
 *   header file for HapMixIndividualCollection class
 *   Copyright (c) 2006, 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef HAPMIX_INDIVIDUAL_COLLECTION_H
#define HAPMIX_INDIVIDUAL_COLLECTION_H 1

#include "IndividualCollection.h"
#include "HapMixIndividual.h"
class HapMixOptions;
class HapMixFreqs;

///Class to hold an array of Individuals for hapmixmodel
class HapMixIndividualCollection : public IndividualCollection
{
public:
  ///constructor
  HapMixIndividualCollection(const HapMixOptions* const options, InputData* const Data, Genome* Loci, const double* theta);
  //destructor
  ~HapMixIndividualCollection();
  ///samples hidden states for each individual and calculates sufficient statistics for population-level parameters
  void SampleHiddenStates(const HapMixOptions* const options, unsigned iteration);
  ///returns sufficient statistics for update of arrival rates
  const int* getConcordanceCounts()const;
  ///returns sufficient statistics for update of mixture proportions
  const int* getSumArrivalCounts()const;
  const HapMixIndividual* getHapMixIndividual(int num)const;
  ///determines if individual i is a case/control ie its genotype came from ccgenotypesfile
  bool isCaseControl(unsigned i)const;
  int getNumberOfIndividualsForScoreTests()const;
  unsigned int getFirstScoreTestIndividualNumber()const;
  void AccumulateConditionalGenotypeProbs(const HapMixOptions* const options, const Genome& Loci);
  void OutputCGProbs(const char* filename, const Vector_s& LocusLabels);
  double getDevianceAtPosteriorMean(const Options* const options, vector<Regression *>&R, Genome* Loci, LogWriter &Log,
				    const double* const MixtureProps, const vector<double>& SumRho, unsigned numChromosomes, AlleleFreqs* A );
private:
  unsigned NumCaseControls;
  int* ConcordanceCounts;//for update of locus-specific sumintensities
  int* GlobalConcordanceCounts;//ConcordanceCounts summed over processes, kept on master processes
  int* SumArrivalCounts;
  int* GlobalSumArrivalCounts;
  GenotypeProbOutputter GPO;

  std::vector<HapMixIndividual*> HapMixChild; 

  bool isMaskedIndividual(unsigned i, const vector<unsigned>& maskedIndividuals)const;
  //default c'tor not implemented
  HapMixIndividualCollection();

};

#endif /* !defined HAPMIX_INDIVIDUAL_COLLECTION_H */


