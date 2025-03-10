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
/// \file HiddenMarkovModel.new.h
/// Definition of the genepi::HiddenMarkovModel class.
//=============================================================================

#ifndef __base_HiddenMarkovModel_new_h
#define __base_HiddenMarkovModel_new_h



#include "HiddenStateSpace.h"
#include "Pedigree.h"
#include "TransProbCache.h"
#include "HVIterator.h"
#include <bclib/cvector.h>



/// Should the forwards and backwards recursions be done in parallel?  When
/// using the "big cache" model for TransProbCache, this has almost no
/// advantage, since most of the computation is in computing the transaction
/// probabilities, which must be done somewhat atomically.
#define HMM_PARALLELIZE_FWD_BKWD	0



namespace genepi { // ----


/** \addtogroup base
 * @{ */



/// Hidden Markov Model to be used with pedgrees.  This differs from the class
/// used with individuals only in the namespace (genepi:: vs. ::).  This is
/// confusing and should be changed.

class HiddenMarkovModel
    {

    protected:
	typedef double ProbType;
	typedef size_t IntSLocIdx;

	/// Probability array type, for the forward (alpha) and backwards(beta)
	/// probabilities: indexed on locus-index (t); within each locus,
	/// indexed on hidden-state-index.
	typedef cvector< ProbType	  > ProbsAtLocusType ;
	typedef cvector< ProbsAtLocusType > ProbArrType	     ;
	typedef cvector< ProbType	  > PiType	     ; ///< Indexed on "overall" state-index (including 0-prob states)
	typedef Pedigree::AdmixType	    AdmixType	     ;

    private:
	const Pedigree * ped	 ; ///< Ref to ped for hidden state space, etc.
	TransProbCache * tpCache ; ///< Transition probabilities.  These should perhaps
				   ///< be retrieved from the Pedigree object itself

	/// Set to CHR_IS_X if the HMM refers to the X chromosome loci only,
	/// set to CHR_IS_NOT_X otherwise.
	IsXChromType	 isX	 ;

	/// Admixture probabilities for the founder-gametes
	const AdmixType * theta;

	// Range of loci to use for this model: typically, either the X
	// chromosome, or all of the non-X chromosomes combined.
	SLocIdxType firstLocus ;
	SLocIdxType nLoci      ;

	mutable bool dirtyForwards ; ///< Input parameters have changed since last compute; stored alpha is invalid
	mutable bool dirtyBackwards; ///< Input parameters have changed since last compute; stored beta is invalid
	mutable bool dirtyPi	   ; ///< Input parameters have changed since last compute; stored Pi is invalid

	mutable ProbArrType alpha ; ///< Forward probability array: indexed on locus, then hidden-state-non0-index
	mutable ProbArrType beta  ; ///< Forward probability array: indexed on locus, then hidden-state-non0-index
	mutable double norm_log_sum_alpha;
	mutable double norm_log_sum_beta;

	/// Cached conditional state probabilities (these are just the
	/// normalized element-wise product of alpha and beta, could be computed
	/// on-the-fly).  Indexed on locus, then hidden-state-non0-index.
	mutable ProbArrType condStateProbs;
	mutable bool	    dirtyCondStateProbs;

	/// Stationary distribution, calculated in computeForwards() and used
	/// in computeBackwards(), indexed on all states
	mutable PiType Pi;

	/// Recursive function to multiply the transition-probability matrix
	/// by the probability vector alpha
	void transRecursion(const double *alpha, double *res,
			    double f, double g, const ThetaType& h,
			    HVIterator z, size_t sz ) const;

	void computeStationaryDistr( const HiddenStateSpace & hss_0,
				     const ThetaType& th ) const;

	void initStaticAlpha( ProbsAtLocusType &       alpha0	,
			      const HiddenStateSpace & hss	,
			      double &		       norm_sum ) const;

	/// Do the main forwards recursion
	void computeForwards() const;

	/// Do the main backwards recursion
	void computeBackwards() const;

	/// Do both the forwards and backwards computations, in parallel if possible.
	void computeForwardsBackwards() const;

    protected:

	/// This is the main recursion method, which effectively multiplies the
	/// probabilties at one locus by the transition probability matrix
	/// (although it may be implemented by an optimized algorithm which
	/// neither calculates the matrix nor actually multiplies it).  It is
	/// used both for the forwards and backwards recursions.
	///	@param f	The probability of 0 arrivals at the
	///			    destination locus (toProbs)
	///	@param g	"g" factor at the destination locus (toProbs)
	///	@param th	The admixture proportions, adjusted to account
	///			    for the female/male odds ratios if on the
	///			    X chromosome.
	///	@param fr_hss	The hidden-state-space at the "from" locus.
	///			    This is necessary to find the values of the
	///			    ancestry-vector and inheritance-vector for
	///			    each hidden-state.
	///	@param to_hss	The hidden-state-space at the "to" locus.
	///	@param frProbs	The "from" probabilities, indexed on the
	///			    HiddenStateSpace's "overallIndex".
	///	@param toProbs	The "to" probabilities, indexed on the
	///			    HiddenStateSpace's "overallIndex".
	void recursionProbs( double		      f		,
			     double		      g		,
			     const ThetaType &	      th	,
			     const HiddenStateSpace & fr_hss	,
			     const HiddenStateSpace & to_hss	,
			     const ProbsAtLocusType & frProbs	,
			     ProbsAtLocusType &	      toProbs	) const;


	void assureNotDirty() const;

	/// The forwards probability array: indexed on locus, then hidden-state-index.
	const ProbArrType & getAlpha() const
	    {
	    if ( dirtyForwards )
		computeForwardsBackwards();
	    return alpha;
	    }

	/// The backwards probability array: indexed on locus, then hidden-state-index
	const ProbArrType & getBeta () const
	    {
	    if ( dirtyBackwards )
		computeForwardsBackwards();
	    return beta;
	    }


    public:

	/// Construct the model.
	/// @parm theta - see setTheta()
	HiddenMarkovModel( const Pedigree & ped , TransProbCache & tpCache ,
			   IsXChromType isX,
			   const AdmixType * theta = 0 );

	const Pedigree &       getPed() const { return *ped    ; } ///< Get reference to Pedigree
	const TransProbCache & getTPC() const { return *tpCache; } ///< Get reference to TransProbCache


	/// Get the number of loci to model.
	SLocIdxType getNLoci() const { return nLoci; }

	/// Translate our internal (relative) locus-index to "external"
	/// (absolute) locus-index.
	SLocIdxType extLoc( IntSLocIdx idx ) const
	    {
	    gp_assert_lt( idx, nLoci );
	    return firstLocus + idx;
	    }

	/// Translate an "external" (absolute) locus-index our internal
	/// (relative) locus-index.
	IntSLocIdx intLoc( SLocIdxType idx ) const
	    {
	    IntSLocIdx rv = idx - firstLocus;
	    gp_assert_lt( rv, nLoci );
	    return rv;
	    }

	const SimpleLocus & locus( IntSLocIdx t ) const
	    {
	    return getPed().getSLoci()[ extLoc(t) ];
	    }

	const HiddenStateSpace & getHSS( IntSLocIdx t ) const
	    {
	    return getPed().getStateProbs( extLoc(t) );
	    }



	void transProbsChanged() const;

	/// Set the founder-gamete theta probabilities, indexed on
	/// PopIdxType (not Pedigree::FounderIdxType).  The model does @b not copy nor take
	/// ownership; the vector must not be destroyed for as long as the model
	/// continues to exist, or until a new theta vector is set.  We
	/// @i really need a copy-on-write vector template.
	void setTheta( const AdmixType * nv );
	void thetaChanged() const;


	/// Get the natural logarithm of the likelihood (probability of the data
	/// given the model).  We cache the forward/backwards probabilities and
	/// only recompute them when the data has changed; but the dot-product
	/// is recomputed every time that this method is called, so avoid
	/// re-calling with abandon for unchanged input data.
	double getLogLikelihood() const;

	/// Get the conditional probability distribution of the hidden states at locus @a t.
	const cvector<double> & getCondStateProbsAtLocus( SLocIdxType t ) const;

    };



/** @} */

} // ---- end namespace genepi



#endif // ! __base_HiddenMarkovModel_new_h
