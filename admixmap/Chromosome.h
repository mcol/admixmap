// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   Chromosome.h 
 *   header file for Chromosome class
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef CHROMOSOME_H
#define CHROMOSOME_H 1

#include "HMM.h"

class Individual;
class AdmixOptions;

class Chromosome
{
public:
  Chromosome(int size,int start, int, bool isx);
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

  void isDiploid(bool b){Diploid = b;};
  bool isDiploid()const{return Diploid;};
  bool isXChromosome()const;
// ****************** Setting of locus correlation, f *************************
  void InitialiseLocusCorrelation(const double rho);
  void SetLocusCorrelation(const double rho);
  void InitialiseLocusCorrelation(const std::vector<double> rho);
  void SetLocusCorrelation(const std::vector<double> rho);
  void SetLocusCorrelation(const std::vector<double> rho_, bool global, bool RandomMating);

// ********** Interface to HMM ****************************************
  void SetGenotypeProbs(double* const GenotypeProbs, bool* const GenotypesMissing);
  void SetStateArrivalProbs(const double* const Admixture, bool RandomMating, bool diploid);
  void UpdateHMMInputs(const double* const Admixture, double* const GenotypeProbs, bool* const GenotypesMissing, 
  		     const AdmixOptions* const options, const std::vector< double > _rho, bool diploid);

  void SampleLocusAncestry(int *OrderedStates);
  std::vector<std::vector<double> > getAncestryProbs(const bool isDiploid, int);
  double getLogLikelihood(const bool isDiploid);
  void SampleJumpIndicators(const int* const LocusAncestry, const unsigned int gametes, 
			    int *SumLocusAncestry, std::vector<unsigned> &SumN, 
			    bool SampleArrivals)const;
private:
  double *Distances;
  unsigned int NumberOfCompositeLoci;
 
  int _startLocus;
  int populations;
  std::string _Label;
  HMM SampleStates;
  bool isX;
  bool Diploid;
  
  // f0 and f1 are arrays of scalars of the form exp(- rho*x), where x is distance between loci
  // With a global rho model, this array is same for all individuals and calculated only once.
  // required to calculate transition matrices 
  double *f; 
  int *CodedStates;//used to sample hidden states from HMM

  // UNIMPLEMENTED
  // to avoid use
 // Private default constructor
  Chromosome();
  Chromosome(const Chromosome&);
  Chromosome& operator=(const Chromosome&);
};

#endif /* !defined CHROMOSOME_H */
