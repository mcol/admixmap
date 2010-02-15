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
/// \file HiddenStateSpace.h
/// Definition of the HiddenStateSpace class.
//=============================================================================

#ifndef __base_HiddenStateSpace_h
#define __base_HiddenStateSpace_h



#include "Pedigree.h"
#include "AncestryVector.h"
#include "InheritanceVector.h"
#include "bclib/exceptions.h"



#define HSS_AV_MOST_SIG		1

#define TRACK_UNVISITED_STATES	0



namespace genepi { // ----

/** \addtogroup base
 * @{ */



//-----------------------------------------------------------------------------
//
/// A representation of the hidden-state-space for one pedigree (at one locus).
/// It is effectively a map or associative-array from the
/// (ancestry-vector,inheritance-vector) tuple to the emission-probability
/// "slot" for that state.  It is implemented here as a simple array containing
/// one element for every possible combination of inheritance-vector and
/// ancestry-vector, the size of which is \f$2^M \cdot K^F\f$, i.e. will grow
/// exponentially with respect to the size of the pedigree but with very fast
/// access times.  This could be replaced by a "sparse" implementation more
/// appropriate for larger pedigrees for which many IVs do not qualify.
///
/// Typically one of these objects exists for each locus.
///
/// Use Iterator to iterate over the states in the space, or indexed-based
/// access via stateAtIdx()
//
//-----------------------------------------------------------------------------

class HiddenStateSpace
    {
    public:
	typedef float ProbType;
	typedef size_t AncestryIdxType	  ;
	typedef size_t InheritanceIdxType ;
	typedef size_t StateIdxType	  ; ///< Indexes all (including 0-EP states)
	typedef size_t Non0IdxType	  ; ///< Indexes only existent (non-0-EP) states

    private:

	/// Had to use a pointer rather than a reference so that could have
	/// default constructor, because ISO C++ does not allow initializers to
	/// array new operator, and we'd like to store these in an array. (see
	/// Pedigree.h).  Likewise, K, N_IVs and N_AVs are not const.
	const Pedigree * ped;
	PopIdx		 K;

	/// These are attributes of the pedigree, but we keep them cached here:
	/// Number of IVs = 2^M where M is the number of meioses in the pedigree
	/// Number of AVs = K^F where K is the number of populations and F is the
	///	number of founder gametes in the pedigree.
	size_t	    N_IVs;
	size_t	    N_AVs; ///< @see { N_IVs }
	Non0IdxType nNon0;

	ProbType *  probs;

	StateIdxType aSize() const { return (N_IVs * N_AVs); }


    protected:

	/// Translate an ancestry-vector and an inheritance-vector into a StateIdxType.
	StateIdxType idxOf( const AncestryVector & av, const InheritanceVector & iv ) const
	    {
	    const unsigned long iv_idx = iv.to_ulong();
	    const unsigned long av_idx = av.to_ulong();

	    gp_assert_lt( iv_idx, N_IVs );
	    gp_assert_lt( av_idx, N_AVs );

	    // We can use either arrangement scheme, but if this is changed,
	    // Iterator::advance() must be reimplemented:
	    #if HSS_AV_MOST_SIG
		return (av_idx * N_IVs) + iv_idx;
	    #else
		return (iv_idx * N_AVs) + av_idx;
	    #endif
	    }


    public:

	HiddenStateSpace( const Pedigree & _ped, PopIdx _K );
	HiddenStateSpace();
	~HiddenStateSpace();


	/// This is distasteful but necessary since we must have a default
	/// constructor (see comment at @a ped).
	void init( const Pedigree & _ped, PopIdx _K );


	const Pedigree & getPed() const { gp_assert(ped!=0); return *ped; }
	PopIdx		 getK  () const { return K; }

	const SimpleLocusArray & getSLoci() const { return getPed().getSLoci(); } ///< Convenience


	/// Get the number of states, including those with 0 probability.  See
	/// also stateAtIdx().
	StateIdxType getNStates () const { return aSize() ; }
	Non0IdxType  getNNon0	() const { return nNon0   ; }


	/// Set all emission-probabilities to zero.
	void resetEmProbsToZero();


	/// Get the emission probabiliy for the state at index @a idx.  It's too
	/// bad that we need to expose the internals in the form of an indexing
	/// scheme, but it's the only easy way to iterate/access for
	/// computations, especially given that we store "parallel arrays".  See
	/// also @link getEProb(const AncestryVector&,const InheritanceVector&)
	/// getEProb(AncestryVector,InheritanceVector) @endlink, which in the
	/// current implementation is less efficient.
	ProbType & getEProb( StateIdxType idx )
	    {
	    ProbType & rv = probs[ idx ];
	    #if TRACK_UNVISITED_STATES
		if ( rv == State::NOT_VISITED )
		    rv = 0.0;
	    #endif
	    return rv;
	    }

	/// const version of @link getEProb(StateIdxType) getEProb() @endlink
	ProbType getEProb( StateIdxType idx ) const
	    {
	    #if TRACK_UNVISITED_STATES
		const ProbType rv = probs[ idx ];
		return (rv == State::NOT_VISITED) ? 0.0 : rv;
	    #else
		return probs[ idx ];
	    #endif
	    }

	/// Add @a amt to the e-prob at @a idx
	void addToEProbAt( StateIdxType idx, double amt )
	    {
	    ProbType & rv = probs[ idx ];
	    #if TRACK_UNVISITED_STATES
		if ( rv == State::NOT_VISITED )
		    rv = 0.0;
	    #endif

	    if ( rv == 0.0 )
		++nNon0;

	    rv += amt;
	    }

	#if TRACK_UNVISITED_STATES
	    bool nodeIsVisited( StateIdxType idx ) { return (probs[idx] != State::NOT_VISITED); }
	#endif


	/// Get the emission probability for the state defined by @a av and @a iv.
	/// See also getEProb(StateIdxType).
	ProbType & getEProb( const AncestryVector & av, const InheritanceVector & iv )
	    {
	    return getEProb( idxOf( av, iv ) );
	    }

	/// Const version; result is not a modifiable reference
	ProbType getEProb( const AncestryVector & av, const InheritanceVector & iv ) const
	    {
	    return getEProb( idxOf( av, iv ) );
	    }

	/// Add @a amt to the e-prob at the state defined by @a av and @a iv.
	/// See also addToEProbAt(StateIdxTyoe).
	void addToEProbAt( const AncestryVector & av, const InheritanceVector & iv, double amt )
	    {
	    addToEProbAt( idxOf(av,iv), amt );
	    }


	//---------------------------------------------------------------------
	/// Representation of a hidden state, used for iteration and index-based access.
	//---------------------------------------------------------------------
	struct State
	    {
	    private:
		#if TRACK_UNVISITED_STATES
		    static const ProbType NOT_VISITED;
		#endif
		friend class HiddenStateSpace;

	    public:
		AncestryVector    av     ;
		InheritanceVector iv     ;
		ProbType	  emProb ; /// Emission probabilitiy
	    };


	//---------------------------------------------------------------------
	//
	/// Iterator class: iterates over the whole state-space, returning only
	/// states with non-zero emission probability.  To use, instantiate for
	/// a space, at which point the iterator will "point to" the first state
	/// in the space with non-zero emission probability.  advance() or
	/// operator++() will move to the next (non-zero) state.  The current
	/// state is returned by getState() or operator*().  Test isFinished()
	/// or operator bool() to check if the iterator has moved past the last
	/// element.
	///
	/// This implementation iterates over AncestryVector first, then
	/// InheritanceVector within that (i.e. InheritanceVector is the inner
	/// loop) because, all other things being equal, current implementation
	/// of AVs are more difficult to "increment".  Perhaps it should be
	/// possible to specify the order of iteration, or implement as two
	/// interators, one "nested" inside the other?
	///
	/// @warning
	/// isFinished() must be checked after instantiation, before the first
	/// call to getState() or operator*(), in case there are no non-zero
	/// states in the space!  Also checked after each call to advance() to
	/// assure that the iterator has not moved "past the end" of the space,
	/// although the return-value from advance() also serves.
	///
	/// Example:
	/// @code
	///	void func( const HiddenStateSpace & space )
	///	    {
	///	    for ( HiddenStateSpace::Iterator it( space ); it; ++it )
	///		do_something( *it );
	///	    }
	/// @endcode
	/// Equivalent but more verbose:
	/// @code
	///	for ( HiddenStateSpace::Iterator it( space ); ! it.isFinished(); it.advance() )
	///	    do_something( *it );
	/// @endcode
	///
	/// Check if the space has any states with non-zero emission probabilities:
	/// @code
	///	if ( HiddenStateSpace::Iterator(space).isFinished() )
	///	    cout << "Mendelian inconsistency.\n";
	/// @endcode
	//
	//---------------------------------------------------------------------

	class Iterator
	    {
	    private:
		const HiddenStateSpace &    space;
		AncestryVector::Iterator    av_it;
		InheritanceVector::Iterator iv_it;
		StateIdxType		    sIdx;
		Non0IdxType		    non0Idx;

		/// Whether the iterator should skip the elements for which
		/// the emission probability is zero.
		bool skipNon0;

		/// Can avoid keeping this finished flag by implementing
		/// isFinished() as "return (idx < space.aSize())";
		bool finished;

	    public:

		/// By default, the iterator skips the elements for which
		/// the emission probability is zero.  If @a sparse is set
		/// to false, then those elements will not be skipped.
		Iterator( const HiddenStateSpace & sp, bool sparse = true );
		Iterator( const Iterator & rhs );

		const HiddenStateSpace & getSpace() const { return space; }


		/// Returns true if there are no more states, i.e. if we've
		/// advanced past the last one.
		bool isFinished() const { return finished; }
		operator bool() const { return (! isFinished()); } ///< Synonym for (!isFinished())


		/// Move to the next state depending on the value of skipNon0.
		/// Returns the new value for isFinished().  Throws an exception
		/// if already past the last state, i.e. if already
		/// isFinished().
		bool advance();
		Iterator & operator++() { advance(); return *this; } ///< Synonym for advance()

		/// Make a copy, advance() the copy, and return it, leaving this iterator unchanged.
		Iterator plusOne() const
		    {
		    Iterator rv( *this );
		    return ++rv;
		    }


		/// Retrieve the value of the state currently "pointed to".
		State getState() const;
		State operator* () const { return getState(); } ///< Synonym for getState()
		Iterator *	 operator->()	    { return this; }
		const Iterator * operator->() const { return this; }

		/// Get the AncestryVector of the state currently "pointed to";
		/// this is typically more efficient than extracting the
		/// components from the State returned by
		/// getState()/operator*().
		const AncestryVector &	  getAV	  () const { return av_it.getAV(); }
		const InheritanceVector & getIV   () const { return iv_it.getIV(); }
		double			  getEProb() const { return space.getEProb(sIdx); }


		/// Provide access to the "index" of the state currently
		/// "pointed to".  This is the index based on the overall
		/// potential number of states at this locus for this pedigree,
		/// *not* index on the states the IV of which is consistent with
		/// the genotyped data (for that, see getNon0Index()).  This is
		/// an inelegant exposure of the internal implementation of the
		/// iterator, but because we @b do use the index-style access
		/// (in part because it enables us to keep "parallel arrays"),
		/// providing this bridge between the two methods of access
		/// allows us to iterate with Iterator at times when we need the
		/// index for access to external parallel arrays; otherwise
		/// would have to do all iteration via integral indices.  The
		/// value returned here can be passed to
		/// HiddenStateSpace::stateAtIdx(), although
		/// getState()/operator*() does the same, but more efficiently.
		StateIdxType getOverallIndex() const { return sIdx; }

		/// Provide access to the "index" of the state currently
		/// "pointed to".  This is the index based on only those states
		/// at this locus for this pedigree the IV of which is
		/// consistent with the genotyped data, *not* indexed on the
		/// potential number of states (for that, see
		/// getOverallIndex()).  This index is really only useful
		/// outside of this class, (e.g. for indexing into the
		/// transition-probability matrix); all of the indexes that are
		/// passed back into this class are StateIdxType.
		Non0IdxType getNon0Index() const { gp_assert( skipNon0 ); return non0Idx; }
	    };



	//-----------------------------------------------------------------------------
	// Index-based access:
	//-----------------------------------------------------------------------------

	AncestryVector	  avFromIdx( AncestryIdxType	aIdx ) const;
	InheritanceVector ivFromIdx( InheritanceIdxType iIdx ) const;

	/// Indexed-based access to states.  See also getNStates().  Consider
	/// using Iterator instead.
	State stateAtIdx( StateIdxType ) const;

    };



/** @} */

} // ---- end namespace genepi



#endif // ! __base_HiddenStateSpace_h
