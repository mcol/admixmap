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
#include "GeneticDistanceUnit.h"
#include "InputData.h"

///Container class for Chromosome and CompositeLocus objects.
class Genome{
public:

  Genome();
  Genome(int);
  virtual ~Genome();

  void Initialise(const InputData* const data_, int populations, LogWriter &Log);

  const vector<int> GetChrmAndLocus(int) const;
  const int getChromosomeNumber(int) const;
  const int getRelativeLocusNumber(int) const;

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
  unsigned isXChromosome(unsigned)const;
  bool isXLocus(unsigned j)const;

  double GetDistance(int)const;

  void GetChromosomes(int);
  const Chromosome* const* getChromosomes()const;
  Chromosome *getChromosome(unsigned);

  void PrintLocusTable(const char* filename, const std::vector<double>& Distances)const;

  unsigned GetChrNumOfLocus(unsigned locus); 
  int GetNumberOfStates()const;
  int GetNumberOfStates(int locus)const;
  
  double GetLengthOfGenome()const;
  double GetLengthOfXchrm()const;
  virtual void SetLocusCorrelation(const std::vector<double>& rho);
  void SetLocusCorrelation(double rho);

protected:
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

  virtual void CreateChromosome(unsigned i, unsigned size, bool isX, unsigned cstart, int NumHiddenStates );
  void SetupChromosome(Chromosome* C, bool isX, unsigned cstart, const string& label);

private:
  /** Index of chromosomes and loci.
   *
   * LocusTable[i][0] is the chromosome
   * LucusTable[i][1] is the locus index
   */
  std::vector<std::vector<int> > LocusTable;
  std::vector<std::string> ChrmLabels;

  void InitialiseChromosomes(const std::vector<unsigned> cstart, int populations);
  void PrintSizes(LogWriter &Log, GeneticDistanceUnit u)const;


  // UNIMPLEMENTED
  // to avoid use
  Genome(const Genome&);
  Genome& operator=(const Genome&);
 
};

#endif /* !GENETIC_ARRAY_H */
