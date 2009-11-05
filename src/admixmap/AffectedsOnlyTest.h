// *-*-C++-*-*
/** 
 *   ADMIXMAP
 *   AffectedsOnlyTest.h 
 *   header file for AffectedsOnly class
 *   Copyright (c) 2006, 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef AFFECTEDSONLYTEST_H
#define AFFECTEDSONLYTEST_H 1

#include "ScoreTestBase.h"
#include "bclib/RObjectWriter.h"

#include "config.h"		// AGGRESSIVE_RANGE_CHECK
#include "SimpleLocusArray.h"	// SLocIdxType
#include "Organism.h"		// PopIdx
#include "bclib/estr.h"
using genepi::SLocIdxType;
using genepi::PopIdx;

class Genome;
namespace bclib{
  class LogWriter;
}


/** \addtogroup admixmap
 * @{ */


/**
   Class to implement affecteds-only score test for linkage with locus ancestry
 */
class AffectedsOnlyTest : public ScoreTestBase{

public:
  AffectedsOnlyTest();
  ~AffectedsOnlyTest();

  void Initialise(const char* filename, const int NumPopulations, const int NumLoci);

  void Reset();
  void Update(unsigned int locus, int k0, const double* const Theta, 
	      bool RandomMatingModel, bool diploid, const std::vector<std::vector<double> > AProbs);

  void Accumulate();

  void Output(const Vector_s& PopLabels, const Genome& Loci);
  void WriteFinalTable(const char* filename, const Vector_s& PopLabels, 
		       const Genome& Loci, bclib::LogWriter& Log);
  void OutputLikRatios(const char* const filename, const Vector_s& PopLabels, const Genome& Loci);



  //-------------------------------------------------------------------------
  // Expose internal data structures so pedigrees can directly update them:
  //-------------------------------------------------------------------------

  /// Range-check indexes.
  void rc( SLocIdxType t, PopIdx k ) const
    {
    #if AGGRESSIVE_RANGE_CHECK
	if ( t >= L )
	    throw std::runtime_error( genepi::estr("Locus-index (") + t + ") out-of-range (" + L + ')' );
	if ( k >= K )
	    throw std::runtime_error( genepi::estr("Population-index (") + k + ") out-of-range (" + K + ')' );
    #else
	if ( t|k ) {;} // Suppress compiler warning
    #endif
    }
  double & getAffectedsScore	( SLocIdxType t, PopIdx k ) { rc(t,k); return AffectedsScore   [ t * K + k ]; }
  double & getAffectedsVarScore ( SLocIdxType t, PopIdx k ) { rc(t,k); return AffectedsVarScore[ t * K + k ]; }
  double & getAffectedsInfo	( SLocIdxType t, PopIdx k ) { rc(t,k); return AffectedsInfo    [ t * K + k ]; }



  //-------------------------------------------------------------------------
  // Private members:
  //-------------------------------------------------------------------------

private:
  unsigned K, L;
  unsigned firstpoplabel;
  double* SumAffectedsScore2;
  double* SumAffectedsVarScore;
  double* SumAffectedsScore;
  double* SumAffectedsInfo;

  //score test objects, static so they can accumulate sums over individuals
  double *AffectedsScore;
  double *AffectedsVarScore;
  double *AffectedsInfo;

  double *LikRatio1;
  double *LikRatio2;

  bclib::RObjectWriter R;

  void OutputAffectedsOnlyTest(bclib::DelimitedFileWriter& outfile, const Vector_s& PopLabels, 
			       const Genome& Loci, const std::string& sep, bool final);
};


/** @} */


#endif /* !defined SCORETESTS_H */
