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
#include <bclib/cvector.h>


/// Should the state probabilities be normalized "on the fly" (i.e. after each
/// iteration of the forwards-backwards recursion algorithm?  In practice, they
/// underflow almost immediately if this is not turned on.
#define HMM_OTF_RENORM			1

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

	/// Probability array type, for the forward (alpha) and backwards(beta)
	/// probabilities: indexed on locus-index (t); within each locus,
	/// indexed on hidden-state-index.
	typedef cvector< ProbType	  > ProbsAtLocusType ;
	typedef cvector< ProbsAtLocusType > ProbArrType	     ;
	typedef Pedigree::ThetaType	    ThetaType	     ;

    private:
	const Pedigree * ped	 ; ///< Ref to ped for hidden state space, etc.
	TransProbCache * tpCache ; ///< Transition probabilities.  These should perhaps
				   ///< be retrieved from the Pedigree object itself

	/// Theta probabilities for the founder-gametes, indexed on PopIdxType
	/// AKA Theta
	const ThetaType * theta;

	mutable bool dirtyForwards ; ///< Input parameters have changed since last compute; stored alpha is invalid
	mutable bool dirtyBackwards; ///< Input parameters have changed since last compute; stored beta is invalid

	mutable ProbArrType alpha ; ///< Forward probability array: indexed on locus, then hidden-state-non0-index
	mutable ProbArrType beta  ; ///< Forward probability array: indexed on locus, then hidden-state-non0-index
	#if HMM_OTF_RENORM
	    mutable double norm_log_sum_alpha;
	    mutable double norm_log_sum_beta;
	#endif

	/// Cached conditional state probabilities (these are just the
	/// normalized element-wise product of alpha and beta, could be computed
	/// on-the-fly).  Indexed on locus, then hidden-state-non0-index.
	mutable ProbArrType condStateProbs;
	mutable bool	    dirtyCondStateProbs;

	/// Do the main forwards recursion
	void computeForwards() const;

	/// Do the main backwards recursion
	void computeBackwards() const;

	/// Do both the forwards and backwards computations, in parallel if possible.
	void computeForwardsBackwards() const;

    protected:

	void recursionProbs(  );

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
			    const ThetaType * theta = 0 );

	const Pedigree &       getPed() const { return *ped    ; } ///< Get reference to Pedigree
	const TransProbCache & getTPC() const { return *tpCache; } ///< Get reference to TransProbCache


	/// Get the number of loci to model.
	SLocIdxType getNLoci() const { return getPed().getSLoci().size(); }


	void transProbsChanged() const;

	/// Set the founder-gamete theta probabilities, indexed on
	/// PopIdxType (not Pedigree::FounderIdxType).  The model does @b not copy nor take
	/// ownership; the vector must not be destroyed for as long as the model
	/// continues to exist, or until a new theta vector is set.  We
	/// @i really need a copy-on-write vector template.
	void setTheta( const ThetaType * nv );
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
