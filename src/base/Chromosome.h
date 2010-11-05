// *-*-C++-*-*
/*
 *   Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

//=============================================================================
/// \file Chromosome.h
/// Definition of the Chromosome class.
//=============================================================================

#ifndef CHROMOSOME_H
#define CHROMOSOME_H 1


#include <vector>
#include <string>
#include "HiddenMarkovModel.h"
#include "RhoType.h"


/** \addtogroup base
 * @{ */


/// Represents a chromosome and holds HMM object
class Chromosome
{
public:
  ///HMM Object (should be private for safety, use getHMM() instead):
  HiddenMarkovModel * HMM;

public:
  Chromosome(int, int size,int start, int, bool isx);
   ~Chromosome();
// ******** Chromosome information **********************************
  void SetLabel(std::string );

  /// Set the distance from locus @a locus to the _previous_ locus.  This implies
  /// that @a locus is >=1.
  void SetDistance( unsigned int locus, double dist );

  const std::string GetLabel() const;


  /// Returns the number of the num'th compositelocus on this chromosome
  /// eg if chromosome 2 consists of compositeloci 5, 6, 7 and 8,
  /// GetLocus(i) returns 5 + i.
  int GetLocus(int) const;


  unsigned int GetSize			() const { return NumberOfCompositeLoci; }
  unsigned int GetNumberOfCompositeLoci () const { return NumberOfCompositeLoci; }


  /// Get the distance from locus @a locus to the _previous_ locus.  This implies
  /// that @a locus is >=1.
  double GetDistance( unsigned int locus ) const
    {
    gp_assert_gt( locus, 0			    );
    gp_assert_lt( locus, GetNumberOfCompositeLoci() );

    return Distances[ locus ];
    }

  bool isXChromosome()const;
  // ****************** Setting of locus correlation, f *************************
  void SetGlobalLocusCorrelation(const double rho);
  void SetLocusCorrelation( const genepi::RhoType & vrho, bool RandomMating );

  // ********** Interface to HMM ****************************************

  /// Get the conditional ancestry probabilities at the specified locus
  void getHiddenStateCopyNumberProbs(std::vector<std::vector<double> >& AProbs,
                                     const bool isDiploid, int locus);
  ///
  void SampleJumpIndicators(const int* const LocusAncestry, const unsigned int gametes, 
			    int *SumLocusAncestry, std::vector<unsigned> &SumN, 
			    bool SampleArrivals)const;

  HiddenMarkovModel &	    getHMM()	   { return *HMM; }
  const HiddenMarkovModel & getHMM() const { return *HMM; }


private:

  /// Compute the locus correlation at a given locus given the value of the
  /// sum-intensities parameter rho
  double LocusCorrelation(unsigned locus, double drho);

  // f0 and f1 are arrays of scalars of the form exp(- rho*x), where x is distance between loci
  // With a global rho model, this array is same for all individuals and calculated only once.
  // required to calculate transition matrices 
  double *f; 
  const unsigned int NumberOfCompositeLoci;
  const bool isX;

  double *Distances;
 
  const int Number;//number of chromosome
  const int _startLocus;
  const int NumHiddenStates; ///< Not really the number of hidden states -- we're just funning you!
  std::string _Label;
  
  int *CodedStates;//used to sample hidden states from HMM

  // UNIMPLEMENTED
  // to avoid use
 // Private default constructor
  Chromosome();
  Chromosome(const Chromosome&);
  Chromosome& operator=(const Chromosome&);
};


/** @} */


#endif /* !defined CHROMOSOME_H */

