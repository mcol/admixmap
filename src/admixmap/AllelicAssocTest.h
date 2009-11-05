// *-*-C++-*-*
/** 
 *   AllelicAssocTest.h 
 *   header file for AllelicAssocTest class
 *   Copyright (c) 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef ALLELICASSOCTEST_H
#define ALLELICASSOCTEST_H 1

#include "ScoreTestBase.h"
#include "AdmixOptions.h"
#include "IndividualCollection.h"
#include "Genome.h"
#include "AllelicAssocSubTest.h"
#include "bclib/RObjectWriter.h"


/** \addtogroup admixmap
 * @{ */


/**
 *  Class to implement 
 *   (1) Score test for allelic association
 *   (2) Score test for within-haplotype association
 */
class AllelicAssocTest : public ScoreTestBase{

public:
  AllelicAssocTest();
  ~AllelicAssocTest();

  void Initialise( AdmixOptions* , const IndividualCollection* const, const Genome* const ,
		  bclib::LogWriter &);

  void Reset();

  void Output();
  void WriteFinalTables(bclib::LogWriter& Log);

  void MergeRareHaplotypes(const std::vector<double> &alpha0);

  void Update( const PedBase * , double, double, double, bool);

  void Accumulate();


private:
  std::vector<AllelicAssocSubTest*> SubTests;

  std::vector<HaplotypeTest*> HaplotypeAssocTests; 

  int *locusObsIndicator;

  bclib::RObjectWriter AllelicAssocRObject;
  bclib::RObjectWriter HaplotypeAssocRObject;
  //outputfile = allelicAssocStream
  //std::ofstream allelicAssocScoreStream;

  const AdmixOptions *options;
  const IndividualCollection *individuals;
  const Genome* Lociptr;
  unsigned NumCompositeLoci;
  std::vector<unsigned> NumLoci;//number of loci within each comp locus
  std::vector<unsigned> NumMergedHaplotypes;
  const Chromosome* const* chrm;//Copy of pointer to array of chromosomes
  unsigned NumOutputs;//counts calls to output function for dimensions of R objects
  void ROutput();
};


/** @} */


#endif /* !defined ALLELICASSOCTEST_H */
