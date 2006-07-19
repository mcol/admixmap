// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   Genome.h
 *   header file for Genome class (formerly known as GeneticArray)
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#ifndef GENETIC_ARRAY_H
#define GENETIC_ARRAY_H 1

#include "CompositeLocus.h"
#include "Chromosome.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "AdmixOptions.h"
#include "InputData.h"

///Container class for Chromosome and CompositeLocus objects.
class Genome 
{
public:

  Genome();
  Genome(int);
  virtual ~Genome();

  void Initialise(const InputData* const data_, int populations, LogWriter &Log, int rank);

  const std::vector< int > GetChrmAndLocus( int )const;

  const std::vector<  std::vector< int > >GetChrmAndLocus( )const;
  void GetChrmAndLocus(unsigned locus, unsigned* c, unsigned* l);

  bool isX_data()const;

  CompositeLocus* operator()(int) const;

  void SetLabels(const std::vector<std::string> &labels, const std::vector<double> &distances);

  const double *GetDistances()const;

  unsigned int GetNumberOfCompositeLoci()const;

  unsigned int GetNumberOfChromosomes()const;

  unsigned int GetTotalNumberOfLoci()const;

  int getNumberOfLoci(int)const;

  const unsigned int *GetSizesOfChromosomes()const;
  unsigned GetSizeOfChromosome(unsigned)const;

  unsigned getFirstXLocus()const;
  unsigned isXChromosome(unsigned);

  double GetDistance(int)const;

  void GetChromosomes(int);
  const Chromosome* const* getChromosomes()const;
  Chromosome* getChromosome(unsigned);

  void PrintSizes(LogWriter &Log, const string unit)const;
  void PrintLocusTable(const char* filename, const std::vector<double>& Distances)const;

  unsigned GetChrNumOfLocus(unsigned locus); 
  int GetNumberOfStates()const;
  int GetNumberOfStates(int locus)const;
  
  double GetLengthOfGenome()const;
  double GetLengthOfXchrm()const;
  void SetLocusCorrelation(const std::vector<double> rho);
  void SetLocusCorrelation(double rho);

private:
  Chromosome **C;
  double *Distances;
  unsigned int NumberOfCompositeLoci;
  CompositeLocus *LocusArray; 
  double LengthOfGenome;
  double LengthOfXchrm;

  unsigned int NumberOfChromosomes;
  unsigned int TotalLoci;//number of simple loci;
  unsigned int *SizesOfChromosomes;
  bool X_data;
  unsigned XChromosomeIndex;
  std::vector< std::vector< int > > LocusTable;
  std::vector<std::string> ChrmLabels;

  void InitialiseChromosomes(const vector<unsigned> cstart, int populations);

  // UNIMPLEMENTED
  // to avoid use
  Genome(const Genome&);
  Genome& operator=(const Genome&);
 
};

#endif /* !GENETIC_ARRAY_H */
