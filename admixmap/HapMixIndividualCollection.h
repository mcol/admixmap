// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   HapMixIndividualCollection.h 
 *   header file for HapMixIndividualCollection class
 *   Copyright (c) 2006 David O'Donnell and Paul McKeigue
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

///Class to hold an array of Individuals for hapmixmodel
class HapMixIndividualCollection : public IndividualCollection
{
public:
  ~HapMixIndividualCollection();
  HapMixIndividualCollection(const AdmixOptions* const options, const InputData* const Data, Genome* Loci);
  void SampleLocusAncestry(const AdmixOptions* const options);
  const int* getSumAncestry()const;
  //Individual* getIndividual(int)const;
  int getNumberOfIndividualsForScoreTests()const;
  int getFirstScoreTestIndividualNumber()const;
  void AccumulateConditionalGenotypeProbs(const AdmixOptions* const options, const Genome& Loci);
  void OutputCGProbs(const char* filename);
private:
  unsigned NumCaseControls;
  int* SumAncestry;//for update of locus-specific sumintensities
  int* GlobalSumAncestry;//SumAncestry summed over processes, kept on master processes
  GenotypeProbOutputter GPO;
};

#endif /* !defined HAPMIX_INDIVIDUAL_COLLECTION_H */


