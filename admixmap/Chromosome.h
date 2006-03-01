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

  void InitialiseLociCorr(const double rho);
  void SetLociCorr(const double rho);
  void InitialiseLociCorr(const std::vector<double> rho);
  void SetLociCorr(const std::vector<double> rho);

  void UpdateHMMForwardProbs(const double* const Admixture, double* const GenotypeProbs, bool* const GenotypesMissing, 
			     const AdmixOptions* const options, const std::vector< double > _rho, bool diploid);
  void UpdateHMMBackwardProbs(const double* const hapAdmixture, const double* const GenotypeProbs);  
  //call only after a call to UpdateHMMForwardProbs
  //this is ok as whenever we need backward probs we also need forward probs but not vice versa


  void SampleLocusAncestry(int *OrderedStates, const double* const Admixture)const;
  std::vector<std::vector<double> > getAncestryProbs(const bool isDiploid, int)const;
  double getLogLikelihood(const bool isDiploid)const;
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
