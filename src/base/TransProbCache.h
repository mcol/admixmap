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
/// \file TransProbCache.h
/// Definition of the TransProbCache class.
//=============================================================================

#ifndef __base_TransProbCache_h
#define __base_TransProbCache_h


#include "config.h"	// AGGRESSIVE_RANGE_CHECK
#include "bclib/exceptions.h"
#include "Pedigree.h"
#include "HiddenStateSpace.h"
#include <bclib/cvector.h>


#define TPC_NO_CACHE    0 ///< Naive implementation: don't cache anything, compute everything every time
#define TPC_BIG_CACHE   1 ///< Slightly less naive implementation: cache absolutely everything as separate matrices

#define TPC_CACHE_MODEL TPC_NO_CACHE



namespace genepi { // ----

/** \addtogroup base
 * @{ */



//-----------------------------------------------------------------------------
//
/// Class to calculate and store transition probabilities for the Hidden Markov
/// Model.
///
/// As currently implemented, we do not cache these probabilities but rather just
/// calculate them on demand, so this class is perhaps more properly called
/// "TransProbCalculator".
///
/// The transition probabilities are retrieved via a <I>from</I> and <I>to</I>
/// state, which are on locus <I>t</I> and <I>t+1</I>, respectively.  The
/// <I>from</I> and <I>to</I> states are referenced as indexes into a particular
/// HiddenStateSpace, i.e. the hidden-state-spaces for a particular Pedigree at
/// locus <I>t</I> and <I>t+1</I>, which are themselves indexes into the
/// SimpleLocusArray.
///
/// The same class both calculates and stores the transition probabilities, and
/// does so in one single object for all pedigrees and loci.  This gives us
/// maximum flexibility in implementing an encapsulated cache model.  At
/// current, two cache models are implemented via conditional compilation:
/// #TPC_BIG_CACHE, which calculates and stores all transition probabilities
/// when the object is constructed, and #TPC_NO_CACHE, which stores no
/// probabilities but (re)calculates them on demand with every request.
///
/// Since we are glomming together the probabilities for multiple pedigrees into
/// a single object while they are likely to be accessed in rapid succession (by
/// a single thread) for a single pedigree, we provide the PedCache nested
/// class, which is an abstraction for (and typically an implementation of) a
/// pedigree-specific cache or table.  This is both for convenience to the
/// calling code and for potential efficiencies in index dereferencing.  Call
/// getPedCache() to get a reference to the PedCache object for a particular
/// pedigree, then make repeated calls to PedCache::getProb() to get transition
/// probabilities.
///
/// If not accessing many probabilities in succession for the same pedigree,
/// TransProbCache::getProb() can be used.
//
//-----------------------------------------------------------------------------

class TransProbCache
    {
    private:
	typedef float CacheProbType;

    public:

	#if 0
	    typedef TwoDimArray<Pedigree::FounderIdxType,PopIdx,double> MuType;
	#else
	    typedef Pedigree::ThetaType MuType;
	#endif

    private:


	#if TPC_CACHE_MODEL == TPC_BIG_CACHE

	    CacheProbType * probs	 ;
	    size_t	    locMultiplier;
	    size_t	    stMultiplier ;
	    #if AGGRESSIVE_RANGE_CHECK
		size_t	    n_probs;
	    #endif

	#endif


	const Pedigree &	 pedigree ;
	const SimpleLocusArray & loci     ;

	/// Ancestry probabilities for the founder-gametes, indexed on PopIdxType
	const cvector<double> * ancestry;


	// For the moment, we are caching transition-probabilities for all pedigrees
	// in the same cache; we therefore get the HSS from each Pedigree object.
	// Alternatively, we could use once TransProbCache for each Pedigree, in
	// which case it might be advantageous to decouple the HiddenStateSpace from
	// the Pedigree object as follows, except that seems like total rubbish
	// since the HSS is dependent on the pedigree's structure anyhow.
	#if 0
	    const HiddenStateSpace * const hss; ///< indexed on locus-idx (t)
	#endif


	#if TPC_CACHE_MODEL == TPC_NO_CACHE

	    struct SLocFacts
		{
		double f; ///< Locus correlation factors
		double g;
		};
	    cvector<SLocFacts>	    factors ; ///< size Nloci-1
	    double		    rho	    ; ///< passed into constructor or re-set later via setRho()

	    const MuType *	    mu	    ; ///< proportion of admixture from population x or gamete y.
					      ///< Passed into constructor or set via setMu().  We neither
					      ///< take ownership nor copy, but just keep a referece,
					      ///< so the calling code must keep the object in
					      ///< existence until no longer needed.
	    const MuType & getMu() const { return *mu; }


	#elif TPC_CACHE_MODEL == TPC_BIG_CACHE

	    //-------------------------------------------------------------------------
	    //
	    /// Big cache model: indexed on pedigree-idx, then
	    /// (locus-idx,from-state-idx,to-state-idx) So overall number of
	    /// floating-point probabilities is:
	    ///
	    /// * If all pedigrees have the same graph topology, then:
	    ///	@code N_pedigrees*(N_loci-1)*(2^M*K^F)^2 @endcode
	    /// But since M and F can vary from one pedigree to another, in fact it is:
	    /// \f$(N_loci-1) \cdot \sum_{each of N pedigrees} 2^M \cdot K^F\f$
	    ///
	    /// Within each pedigree (element of the cvector), one big array
	    /// of FP probs is new'd, and addressing is handled arithmetically
	    /// via PedCache::getProb().  For performance, we cache certain data
	    /// about the pedigree in the Block object, although it could just
	    /// as easily be looked up from the Pedigree class.  We simplify by
	    /// assuming that all HiddenStateSpace's of the same Pedigree
	    /// (i.e. at each locus) have the same size.  This is not strictly
	    /// true since some inheritance vectors may not match the genotype
	    /// data regardless of the founder-haplotype state, but current
	    /// indexing scheme allows for this as "empty" (0
	    /// emission-probability) states.  (In future, must generalize
	    /// this?)  This allows, e.g., one fixed "locus multiplier" for
	    /// addressing, without the need to keep separate pointers within
	    /// the block for each locus.
	    //
	    //-------------------------------------------------------------------------

	    PedCache cache;
	#endif



    protected:
	void computeFactors();



    public:

	TransProbCache( const Pedigree & _pedigree, double _rho, const MuType & _mu );
	~TransProbCache();


	/// Set \f$\rho\f$, the arrival rate parameter.
	void setRho( double nv );

	/// Set \f$\mu\f$, the proportion of admixture from each population on
	/// each gamete.  See also muChanged().
	void setMu( const MuType & nv );

	/// Since we keep a reference to the calling code's copy of mu, it is likely
	/// that the higher-level model will change those values in-place, so rather
	/// than hand us a new array (setMu()), it can just call muChanged() to
	/// signal us that the values changed so that any cached transition
	/// probabilities must be invalidated.
	void muChanged();


	//-----------------------------------------------------------------------------
	// Transition probability computuation:
	//-----------------------------------------------------------------------------

	/// Compute the "f" factor from locus t to t+1
	static double computeF( const SimpleLocusArray & loci, SLocIdxType t, double rho );

	/// Compute the "g" factor from locus t to t+1
	static double computeG( const SimpleLocusArray & loci, SLocIdxType t );

	/// Compute the transition probability between two hidden-states on two "adjacent" loci.
	///
	/// Computes the transition probability from hidden-state @a A on locus
	/// @a t to hidden-state @a B on locus @a t+1.  @a A and @a B are
	/// specified as indexes into the space (@a fromIdx and @a toIdx,
	/// respectively: see HiddenStateSpace::stateAtIdx()).
	static double computeProb( const HiddenStateSpace & frSpace, HiddenStateSpace::StateIdxType frIdx,
				   const HiddenStateSpace & toSpace, HiddenStateSpace::StateIdxType toIdx,
				   double f, double g, const MuType & mu );



	//-----------------------------------------------------------------------------
	/// Retrieve the transition probability from state @a i (as an index
	/// into the HiddenStateSpace at @a t) to state @a j (as an index into
	/// the HiddenStateSpace at @a t + 1).
	//-----------------------------------------------------------------------------

	#if TPC_CACHE_MODEL == TPC_NO_CACHE

	    double getProb( SLocIdxType t ,
			HiddenStateSpace::StateIdxType frStIdx ,
			HiddenStateSpace::StateIdxType toStIdx ) const;

	#elif TPC_CACHE_MODEL == TPC_BIG_CACHE

	    double getProb( SLocIdxType t ,
			HiddenStateSpace::StateIdxType frStIdx ,
			HiddenStateSpace::StateIdxType toStIdx ) const
		{
		const size_t offset = (t * locMultiplier) + (frStIdx * stMultiplier) + toStIdx;
		#if AGGRESSIVE_RANGE_CHECK
		    gp_assert_lt( offset, n_probs );
		#endif
		return probs[ offset ];
		}

	    size_t getNProbs() const; ///< Only meaningful for debugging

	#endif

    };



/** @} */

} // ---- end namespace genepi



#endif // ! __base_TransProbCache_h
