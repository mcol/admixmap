//=============================================================================
//
// Copyright (C) 2009  David D. Favro
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 3 as published by the Free
// Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
///
/// \file PedBase.h
///
/// Definition of the PedBase class.
///
//=============================================================================

#ifndef __base_PedBase_h
#define __base_PedBase_h



#include <vector>

#include <bclib/Regression.h>

#include "HapPair.h"
#include "Options.h"
#include "CopyNumberAssocTest.h"
#include <bclib/cvector.h>
#include "RhoType.h"
#include <bclib/pvector.h>

class AdmixOptions; // Too bad, that's what happens when you use horrible design
class AlleleFreqs; // #include "AlleleFreqs.h" // Forward declaration due to circular-dependency
class AffectedsOnlyTest;
class chib;



#define PEDBASE_DEBUG_METHODS	1



namespace genepi { // ----



/** \addtogroup base
 * @{ */



//-----------------------------------------------------------------------------
///
/// Base class for both Individual and AdmixPedigree.  For the moment, all
/// methods throw an exception and must be overridden so that existing code can
/// be simply re-parented here without necessarily implementing every method;
/// later these could be replaced with either pure methods or ones that
/// implement common functionality.
///
/// This allows arrays of objects to be either Pedigree or Individual objects,
/// by containing pointers to PedBase objects.
///
//-----------------------------------------------------------------------------

class PedBase
    {

    // Types:
    public:
	// Was: [NumGametes][NumPops]
	// Now: [NumThetas][NumPops] (???)
	typedef cvector< cvector<double> >	  AlphaType	;
	typedef bclib::pvector<double>		  ThetaElType	;
	typedef cvector< bclib::pvector<double> > ThetaType	; ///< Alternatively, TwoDimArray<FounderIdx,PopIdx,double>
								  ///< Size getNTheta() [effectively getNFounders()] by K
	typedef genepi::cvector< double >	  PopThetaType	;


    public:
	virtual ~PedBase();

	virtual void DeleteGenotypes();
	virtual void HMMIsBad(bool loglikisbad);

	virtual void setOutcome(double*);
	virtual void setCovariates(double*);
	virtual void setGenotypesToMissing();
	virtual void SetMissingGenotypes();

	virtual const double* getAdmixtureProps(bool isXChrom = false) const;
	virtual const std::vector<hapPair> & getPossibleHapPairs( unsigned int locus ) const;
	virtual const int * getSampledHapPair(int locus)const;
	virtual bool GenotypeIsMissing(unsigned int locus)const;///< locus is a comp locus
	virtual bool simpleGenotypeIsMissing(unsigned locus)const;///< locus is a simple locus
	virtual bool isHaploidatLocus(unsigned j)const;
	virtual bool isHaploidIndividual()const;
	virtual bool isPedigree() const = 0;

	virtual double getLogLikelihood( const Options & , bool forceUpdate, bool store );
	virtual double getLogLikelihoodXChr(const Options&, bool forceUpdate, bool store);
	virtual double getLogLikelihoodAtPosteriorMeans(const Options& options);

	/// Called if a Metropolis proposal is accepted.  Public so that can be
	/// called from PopAdmix.  The @A setHMMAsOK parameter is ignored in the
	/// case of a Pedigree, which tries to map this call into calls to
	/// acceptRhoProposal(), and other non-public methods.
	virtual void storeLogLikelihood( const bool setHMMAsOK );

	//--------------------------------------------------------------------------
	// Public rho-proposal methods.  Called from PopAdmix, ignored for individuals.
	virtual void setRho( double nv ) = 0;
	virtual void startRhoProposal () = 0;
	virtual void acceptRhoProposal() = 0;
	virtual void rejectRhoProposal() = 0;
	//--------------------------------------------------------------------------

	virtual void GetLocusAncestry(int locus, int Ancestry[2])const;
	virtual void GetLocusAncestry(int chrm, int locus, int Ancestry[2])const;
	virtual int GetLocusAncestry(int, int, int)const;

	virtual void SampleHiddenStates( const Options & options );

	virtual void SampleHapPair(unsigned chr, unsigned jj, unsigned locus, AlleleFreqs *A,
			  bool skipMissingGenotypes, bool annealthermo, bool UpdateCounts);

	virtual void UpdateAlleleCounts(unsigned j, unsigned jj, unsigned locus, AlleleFreqs *A, bool annealthermo)const;

	virtual void SampleMissingOutcomes(bclib::DataMatrix *Outcome, const std::vector<bclib::Regression*>& R);

	virtual unsigned int getMyNumber() const;
	virtual unsigned int getIndex	() const;
	virtual unsigned int getNumObs	() const = 0;



	// Methods from AdmixedIndividual:
	virtual void drawInitialAdmixtureProps(const AlphaType &alpha);
	virtual void SetGenotypeProbs(int j, int jj, unsigned locus, bool chibindicator);
	virtual void AnnealGenotypeProbs(int j, const double coolness);
	virtual void ResetSufficientStats();
	virtual void SampleJumpIndicators(bool sampleArrivals);
	virtual void SampleTheta( int iteration, double * SumLogTheta, const bclib::DataMatrix * Outcome,
		const DataType * OutcomeType, const std::vector<double> & lambda, int NumCovariates,
		bclib::DataMatrix * Covariates, const std::vector<const double*> & beta, const PopThetaType & poptheta,
		const AdmixOptions & options, const AlphaType & alpha,
		double DInvLink, double dispersion, CopyNumberAssocTest & ancestryAssocTest, bool RW, bool anneal);
	virtual void SampleRho(const AdmixOptions& options, double rhoalpha,
		    double rhobeta, bool updateSumLogRho);
	virtual void SamplePsi(const AdmixOptions& options,
                               const cvector<double>& priormean,
                               const cvector<double>& priorprec,
                               bool updateSumLogPsi);
	virtual void FindPosteriorModes(const AdmixOptions& options, const AlphaType &alpha,
			  double rhoalpha, double rhobeta, AlleleFreqs* A, ofstream &modefile);
	virtual const RhoType & getRho() const;
	virtual double getPsi(int pop) const;
	virtual void resetStepSizeApproximator(int k);
	virtual void UpdateScores(const AdmixOptions& options, bclib::DataMatrix *Outcome, bclib::DataMatrix *Covariates,
		    const vector<bclib::Regression*> & R, AffectedsOnlyTest & affectedsOnlyTest,
		    CopyNumberAssocTest & ancestryAssocTest );
	virtual void setChibNumerator(const AdmixOptions& options, const AlphaType &alpha, double rhoalpha,
		double rhobeta, chib *MargLikelihood, AlleleFreqs *A);
	virtual void updateChib(const AdmixOptions& options, const AlphaType &alpha, double rhoalpha,
		double rhobeta, chib *MargLikelihood, AlleleFreqs *A);
	virtual double getSumrho()const;
	virtual double getLogLikelihoodOnePop();
	virtual double getLogPosteriorTheta()const;
	virtual double getLogPosteriorRho()const;
	virtual double getLogPosteriorAlleleFreqs()const;
	virtual void WritePosteriorMeans(ostream& os, unsigned samples, bool globalrho)const;
	virtual void WritePosteriorMeansXChr(ostream& os, unsigned samples) const;
	virtual void WritePosteriorMeansLoci(ostream& os)const;
	virtual void setOddsRatios(const genepi::cvector<double>& psi);



	// ====== DEBUGGING METHODS: ======
	#if PEDBASE_DEBUG_METHODS
	    virtual void dumpTheta( const char * prefix ) const = 0;
	#endif
    };



/** @} */



} // ---- end namespace genepi



#endif // ! __base_PedBase_h
