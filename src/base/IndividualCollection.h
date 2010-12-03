//=============================================================================
//
// Copyright (C) 2002-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
// Portions copyright (C) 2009  David D. Favro
//
// This program is free software distributed WITHOUT ANY WARRANTY.  You can
// redistribute it and/or modify it under the terms of the GNU General Public
// License, version 2 or later, as published by the Free Software Foundation.
// See the file COPYING for details.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.	 If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file IndividualCollection.h
/// Header file for IndividualCollection class
//=============================================================================


#ifndef __base_IndividualCollection_h
#define __base_IndividualCollection_h


#include "Genome.h"
#include "bclib/DataMatrix.h"
#include "PedBase.h"
#include "bclib/exceptions.h" // gp_assert()
#include <vector>


class AlleleFreqs;
class AlleleFreqSampler;

namespace bclib{
  class LogWriter;
}


using genepi::PedBase;


/** \addtogroup base
 * @{ */



//=============================================================================
/// Class to hold an array of PedBase.
//=============================================================================

class IndividualCollection
{
private:
  IndividualCollection();
  //double* SumEnergy, *SumEnergySq;//to store sum over iters of energy of test ind at each coolness
  void getLabels(const Vector_s& data, Vector_s& labels);

  void LoadCovariates(const InputData*, const Options* const options, bool admixtureAsCovariate);
  void LoadOutcomeVar(const InputData* const);

  std::ofstream EYStream;//output file for expected outcome

protected:
  const unsigned int	NumInd;		///< Some number of PedBases in this container
  unsigned int		size;		///< Some other number of PedBases in this container
  const int		Populations;
  const int		NumCompLoci;
  unsigned		NumDiploidIndividuals;

  PedBase * *	_child; ///< Pointer to _child array: should be private, use getElement() or operator[].

  // Regression Objects
  bclib::DataMatrix	Outcome			;
  int			NumOutcomes		;
  int			NumCovariates		; ///< Covariates including admixture
  int			NumberOfInputCovariates ; ///< Covariates in file
  bclib::DataMatrix	Covariates		; ///< All covariates, including admixture props
  DataType *		OutcomeType		;

  bclib::DataMatrix *	ReportedAncestry	;

  double		SumLogLikelihood	;
  double		SumDeviance		;
  double		SumDevianceSq		;

  virtual void SetNullValues(); ///< This will zero out @a size, but leave NumInd as initialized.


public:
  virtual ~IndividualCollection();
  IndividualCollection( unsigned numElements, unsigned numPopulations, unsigned numCompositeLoci );

  void DeleteGenotypes(bool);
  //virtual void Initialise( const AdmixOptions * options, const Genome * Loci,
		  //const Vector_s& PopulationLabels, bclib::LogWriter &Log) = 0;

  /// Be aware: this is virtual as if it expects to be overridden, yet in
  /// AdmixIndividualCollection, the signature changes:
  virtual void LoadData( const Options & options, const InputData &, bool admixtureAsCovariate );

  virtual void getOnePopOneIndLogLikelihood(bclib::LogWriter &, const Vector_s& );

  void SampleHapPairs(const Options& options, AlleleFreqs *A, const Genome* const Loci,
		      bool skipMissingGenotypes, bool anneal, bool UpdateCounts);

  void AccumulateAlleleCounts(const Options& options, AlleleFreqs *A, const Genome* const Loci,
			      bool anneal);

  int getSize			() const { return size			; }
  int getNumDiploidIndividuals	() const { return NumDiploidIndividuals ; }

  virtual int		getNumberOfIndividualsForScoreTests () const;
  virtual unsigned int	getFirstScoreTestIndividualNumber   () const;


  int GetNumberOfInputCovariates() const { return NumberOfInputCovariates ; }
  int GetNumCovariates		() const { return NumCovariates		  ; }
  const bclib::DataMatrix & getCovariatesMatrix () const { return Covariates; }
  const bclib::DataMatrix & getOutcomeMatrix	() const { return Outcome   ; }

  PedBase & getElement( size_t idx ) const
	{
	gp_assert_lt( idx, size );
	return *(_child[ idx ]);
	}

  #if 0 // Declared here non-virtual, yet only implemented in AdmixIndCol
      const vector<int>	     getSumLocusAncestry (int k) const;
      const vector<int>	     getSumLocusAncestryX(int k) const;
      const vector<unsigned> getSumNumArrivals();
  #endif

  unsigned	 GetSNPAlleleCounts(unsigned locus, int allele)const;
  int		 getNumberOfMissingGenotypes(unsigned locus)const;
  void		 getAlleleCounts(std::vector<int>& counts,
                                 unsigned locus, int pop) const;
  vector<double> getOutcome( size_t j )const;

  double    getOutcome	     ( int j, int ind ) const { return Outcome.get	 ( ind, j ); }
  bool	    isMissingOutcome ( int j, int i   ) const { return Outcome.isMissing ( i  , j ); }

  int	    getNumberOfOutcomeVars() const { return NumOutcomes	   ; }
  DataType  getOutcomeType(int i)    const { return OutcomeType[i] ; }

  double    getLogLikelihood(const Options& options, bool forceupdate);
  double    getEnergy(const Options& options, const vector<bclib::Regression*> &R,
			    const bool & annealed);

  virtual void HMMIsBad( bool b );
  virtual void resetStepSizeApproximators( int );
  };



/** @} */



#endif // ! __base_IndividualCollection_h
