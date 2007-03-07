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
class HapMixOptions;
class HapMixFreqs;

///Class to hold an array of Individuals for hapmixmodel
class HapMixIndividualCollection : public IndividualCollection
{
public:
  ~HapMixIndividualCollection();
  HapMixIndividualCollection(const HapMixOptions* const options, const InputData* const Data, IGenome* Loci, const HapMixFreqs* A);
  void SampleLocusAncestry(const Options* const options);
  const int* getSumAncestry()const;
  //Individual* getIndividual(int)const;
  int getNumberOfIndividualsForScoreTests()const;
  unsigned int getFirstScoreTestIndividualNumber()const;
  void AccumulateConditionalGenotypeProbs(const HapMixOptions* const options, const IGenome& Loci);
  void OutputCGProbs(const char* filename);
  double getDevianceAtPosteriorMean(const Options* const options, vector<Regression *>&R, IGenome* Loci, LogWriter &Log,
				    const vector<double>& SumRho, unsigned numChromosomes, AlleleFreqs* A );
private:
  unsigned NumCaseControls;
  int* SumAncestry;//for update of locus-specific sumintensities
  int* GlobalSumAncestry;//SumAncestry summed over processes, kept on master processes
  GenotypeProbOutputter GPO;
};

#endif /* !defined HAPMIX_INDIVIDUAL_COLLECTION_H */


