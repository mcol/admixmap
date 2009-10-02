// *-*-C++-*-*
/** 
 *   Genome.h
 *   header file for Genome class (formerly known as GeneticArray)
 */

/*
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
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
#include "GeneticDistanceUnit.h"
#include "InputData.h"
#include "config.h"	// AGGRESSIVE_RANGE_CHECK

namespace bclib{
  class LogWriter;
}



/** \addtogroup base
 * @{ */


///Container class for Chromosome and CompositeLocus objects.
class Genome{

protected:
  void throwErr( const std::string & msg ) const;
  void checkCIdx( unsigned int cIdx ) const
    {
    #if AGGRESSIVE_RANGE_CHECK
	if ( cIdx >= NumberOfChromosomes )
	    throwErr( "chromosome range" );
    #else
	if ( cIdx ) {;} // Suppress compiler warning
    #endif
    }

public:

    Genome();
//  Genome(int); ///< Unimplemented?  What is the int?
  virtual ~Genome();

  void Initialise(const InputData* const data_, int populations, bool hapmixmodelindicator, 
		  bclib::LogWriter &Log);

  const vector<int> GetChrmAndLocus(int) const;
  int getRelativeLocusNumber(int) const;

  const std::vector<  std::vector< int > >GetChrmAndLocus( )const;
  void GetChrmAndLocus(unsigned locus, unsigned* c, unsigned* l);

  bool isX_data()const;

  CompositeLocus* operator()(int) const;
  const CompositeLocus* GetLocus(int)const;

  CompositeLocus &	 operator[]( unsigned int locusIdx );
  const CompositeLocus & operator[]( unsigned int locusIdx ) const
	    { return const_cast<Genome*>( this )->operator[]( locusIdx ); }

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
  const Chromosome & getChromosomeRef( unsigned cIdx ) const { checkCIdx(cIdx); return *C[cIdx]; }
  Chromosome &	     getChromosomeRef( unsigned cIdx )	     { checkCIdx(cIdx); return *C[cIdx]; }

  #if USE_GENOTYPE_PARSER
    void PrintLocusTable( const char * filename, const SimpleLocusArray & simpleLoci ) const;
  #else
    void PrintLocusTable(const char* filename, const std::vector<double>& Distances, const std::string& unitString)const;
  #endif

  unsigned GetChrNumOfLocus(unsigned locus); 
  int GetNumberOfStates()const;
  int GetNumberOfStates(int locus)const;
  
  double GetLengthOfGenome()const;
  double GetLengthOfXchrm()const;

  /// @parm rho Vector has length 1 if globalrho option, otherwise indexed on
  ///	    << something, perhaps composite-locus-index >>
  ///	    (note however that is ignored if the length is not 1).
  virtual void SetLocusCorrelation(const genepi::cvector<double>& rho);

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
   * LocusTable[i][1] is the locus index
   */
  std::vector<std::vector<int> > LocusTable;

  #if ! USE_GENOTYPE_PARSER
    std::vector<std::string> ChrmLabels;
  #endif

  #if USE_GENOTYPE_PARSER
  void InitialiseChromosomes(const std::vector<unsigned> & cstart, const std::vector<size_t> & sstart,
			     int populations, const SimpleLocusArray & sLoci );

  #else
    void InitialiseChromosomes(const std::vector<unsigned> cstart, int populations );
  #endif
  void PrintSizes(bclib::LogWriter &Log, const std::string& distanceUnit)const;


  // UNIMPLEMENTED
  // to avoid use
  Genome(const Genome&);
  Genome& operator=(const Genome&);
 
};


/** @} */


#endif /* !GENETIC_ARRAY_H */
