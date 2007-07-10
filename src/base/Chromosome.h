// *-*-C++-*-*
/** 
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
#include <string>
#include "HiddenMarkovModel.h"

/// Represents a chromosome and holds HMM object
class Chromosome
{
public:
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
  void SetGlobalLocusCorrelation(const double rho);
  void SetLocusCorrelation(const std::vector<double>& vrho, bool RandomMating);

  // ********** Interface to HMM ****************************************

  ///
  std::vector<std::vector<double> > getHiddenStateCopyNumberProbs(const bool isDiploid, int);
  ///
  void SampleJumpIndicators(const int* const LocusAncestry, const unsigned int gametes, 
			    int *SumLocusAncestry, std::vector<unsigned> &SumN, 
			    bool SampleArrivals)const;
  
  ///HMM Object, public for convenient access
  HiddenMarkovModel* HMM;

protected:
  // f0 and f1 are arrays of scalars of the form exp(- rho*x), where x is distance between loci
  // With a global rho model, this array is same for all individuals and calculated only once.
  // required to calculate transition matrices 
  double *f; 
  const unsigned int NumberOfCompositeLoci;
  const bool isX;

  double LocusCorrelation(unsigned locus, double drho);

private:
  double *Distances;
 
  const int Number;//number of chromosome
  const int _startLocus;
  const int NumHiddenStates;
  std::string _Label;
  
  int *CodedStates;//used to sample hidden states from HMM

  // UNIMPLEMENTED
  // to avoid use
 // Private default constructor
  Chromosome();
  Chromosome(const Chromosome&);
  Chromosome& operator=(const Chromosome&);
};

#endif /* !defined CHROMOSOME_H */
