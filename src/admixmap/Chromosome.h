// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   Chromosome.h 
 *   header file for Chromosome class
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef CHROMOSOME_H
#define CHROMOSOME_H 1

#include <vector>
#include <iostream>
#include <string>
#include "HMM.h"

using std::vector;
using std::string;

class Individual;
class AdmixOptions;
class GenotypeProbIterator;

/// Represents a chromosome and holds HMM object
class Chromosome
{
public:
  Chromosome();
  Chromosome(int, int size,int start, int, bool isx);
   ~Chromosome();
// ******** Chromosome information **********************************
  void SetLabel(std::string );
  void SetDistance(int,double);
  const std::string GetLabel( )const;
  int GetLocus(int)const;
  unsigned int GetSize()const;
  unsigned int GetNumberOfCompositeLoci()const;
  const double *GetDistances()const;
  double GetDistance(int)const;

  bool isXChromosome()const;
// ****************** Setting of locus correlation, f *************************
  void SetLocusCorrelation(const double rho);
  void SetLocusCorrelation(const std::vector<double>::const_iterator rho);
  void SetLocusCorrelation(const std::vector<double> rho_, bool global, bool RandomMating);

// ********** Interface to HMM ****************************************
  void SetGenotypeProbs(const GenotypeProbIterator& GenotypeProbs, const bool* const GenotypesMissing);
  void SetHMMTheta(const double* const Admixture, bool RandomMating, bool diploid);
  void SetStateArrivalProbs(bool RandomMating, bool isdiploid);

  void SampleLocusAncestry(int *OrderedStates, bool diploid);
  std::vector<std::vector<double> > getAncestryProbs(const bool isDiploid, int);
  double getLogLikelihood(const bool isDiploid);
  void SampleJumpIndicators(const int* const LocusAncestry, const unsigned int gametes, 
			    int *SumLocusAncestry, std::vector<unsigned> &SumN, 
			    bool SampleArrivals)const;
  const vector<double> getHiddenStateProbs(const bool, int);
private:
  double *Distances;
  unsigned int NumberOfCompositeLoci;
 
  int Number;//number of chromosome
  int _startLocus;
  int populations;
  std::string _Label;
  HMM SampleStates;
  bool isX;
  
  // f0 and f1 are arrays of scalars of the form exp(- rho*x), where x is distance between loci
  // With a global rho model, this array is same for all individuals and calculated only once.
  // required to calculate transition matrices 
  double *f; 
  int *CodedStates;//used to sample hidden states from HMM

  // UNIMPLEMENTED
  // to avoid use
 // Private default constructor
  Chromosome(const Chromosome&);
  Chromosome& operator=(const Chromosome&);
};

#endif /* !defined CHROMOSOME_H */
