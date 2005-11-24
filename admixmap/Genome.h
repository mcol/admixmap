// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   Genome.h
 *   header file for Genome class (formerly known as GeneticArray)
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#ifndef GENETIC_ARRAY_H
#define GENETIC_ARRAY_H 1

#include "CompositeLocus.h"
//#include "Chromosome.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "AdmixOptions.h"
#include "InputData.h"

class Chromosome;//declared here to avoid circular includes
                 //not necessary if getChromosomes function moved out
class Genome 
// Object is an array of objects of class CompositeLocus
// Used to loop over composite loci on a single chromosome, or an entire genome  
{

public:

  Genome();
  Genome(int);
  virtual ~Genome();

  const std::vector< int > GetChrmAndLocus( int )const;

  const std::vector<  std::vector< int > >GetChrmAndLocus( )const;

  bool isX_data()const;

  CompositeLocus *&
  operator()(int) const;

  void SetLabels(const std::vector<std::string> &labels, std::vector<double> temp);

  void loadAlleleStatesAndDistances(const AdmixOptions* const options, const InputData* const data_);

  const double *GetDistances()const;

  unsigned int GetNumberOfCompositeLoci()const;

  unsigned int GetNumberOfChromosomes()const;

  unsigned int GetTotalNumberOfLoci()const;

  int getNumberOfLoci(int)const;

  const unsigned int *GetSizesOfChromosomes()const;

  double GetDistance(int)const;

  void SetDistance(int,double);

  Chromosome **GetChromosomes(int);

  void SetSizes(LogWriter &Log);
 
  int GetNumberOfStates()const;
  int GetNumberOfStates(int locus)const;
  
  double GetLengthOfGenome()const;
  double GetLengthOfXchrm()const;

protected:// to make available to Chromosome
  double *Distances;
  unsigned int NumberOfCompositeLoci;
  CompositeLocus **TheArray;
  bool isChromosome;

private: 
  double LengthOfGenome;
  double LengthOfXchrm;

  unsigned int NumberOfChromosomes;
  unsigned int TotalLoci;//number of simple loci;
  unsigned int *SizesOfChromosomes;
  bool X_data;
  std::vector< std::vector< int > > _chrmandlocus;
  std::vector<std::string> ChrmLabels;

  void InitialiseCompositeLoci();

  // UNIMPLEMENTED
  // to avoid use
  Genome(const Genome&);
  Genome& operator=(const Genome&);
 
};

#endif /* !GENETIC_ARRAY_H */
