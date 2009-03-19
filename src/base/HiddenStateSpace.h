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
#include "exceptions.h"



#define HSS_AV_MOST_SIG		1



namespace genepi { // ----



/** \addtogroup base
 * @{ */



//-----------------------------------------------------------------------------
//
/// A representation of the hidden-state-space for one pedigree (at one locus).
/// It is effectively an map or associative-array from the
/// (ancestry-vector,inheritance-vector) tuple to the emission-probability
/// "slot" for that state.  It is implemented here as a simple array containing
/// one element for every possible combination of inheritance-vector and
/// ancestry-vector, the size of which is (2^M)*(K^F), i.e. will grow
/// exponentially with respect to the size of the pedigree but with very fast
/// access times.  This could be replaced by a "sparse" implementation more
/// appropriate for larger pedigrees for which many IVs do not qualify.
///
/// Typically one of these objects exists for each locus.
///
/// Use Iterator to iterate over the states in the space.
//
//-----------------------------------------------------------------------------

class HiddenStateSpace
    {
    public:
	typedef float ProbType;

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
	size_t N_IVs;
	size_t N_AVs; ///< @see { N_IVs }

	ProbType * probs;

	size_t aSize() const { return (N_IVs * N_AVs); }


    public:

	HiddenStateSpace( const Pedigree & _ped, PopIdx _K );
	HiddenStateSpace();
	~HiddenStateSpace();

	/// This is distasteful but necessary since we must have a default
	/// constructor (see comment at @a ped).
	void init( const Pedigree & _ped, PopIdx _K );


	const Pedigree & getPed() const { gp_assert(ped!=0); return *ped; }
	PopIdx		 getK  () const { return K; }


	/// Get the number of states, including those with 0 probability:
	size_t getNStates() const { return aSize(); }

	void resetEmProbsToZero();


	ProbType & getProb( const AncestryVector & av, const InheritanceVector & iv )
	    {
	    const unsigned long iv_idx = iv.to_ulong();
	    const unsigned long av_idx = av.to_ulong();

	    gp_assert_lt( iv_idx, N_IVs );
	    gp_assert_lt( av_idx, N_AVs );

	    // We can use either arrangement scheme; if this is changed,
	    // Iterator::advance() must be reimplemented:
	    #if HSS_AV_MOST_SIG
		return probs[ (av_idx * N_IVs) + iv_idx ];
	    #else
		return probs[ (iv_idx * N_AVs) + av_idx ];
	    #endif
	    }


	/// Const version; result is not a modifiable reference
	ProbType getProb( const AncestryVector & av, const InheritanceVector & iv ) const
	    {
	    return const_cast<HiddenStateSpace*>(this)->getProb( av, iv );
	    }


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
	//---------------------------------------------------------------------

	class Iterator
	    {
	    private:
		const HiddenStateSpace & space;
		AncestryVector		 av;
		InheritanceVector	 iv;
		size_t			 idx;

		/// Can avoid keeping this finished flag by implementing
		/// isFinished() as "return (idx < space.aSize())";
		bool finished;

	    public:

		Iterator( const HiddenStateSpace & sp );


		struct State
		    {
		    AncestryVector    av     ;
		    InheritanceVector iv     ;
		    ProbType	      emProb ; /// Emission probabilitiy
		    };


		/// Returns true if there are no more states, i.e. if we've
		/// advanced past the last one.
		bool isFinished() const { return finished; }
		operator bool() const { return (! isFinished()); } ///< Synonym for (!isFinished())


		/// Move to the next state with non-zero emission probability.
		/// Returns the new value for isFinished().
		bool advance();
		Iterator & operator++() { advance(); return *this; } ///< Synonym for advance()


		/// Retrieve the value of the state currently "pointed to".
		State getState() const;
		State operator*() const { return getState(); } ///< Synonym for getState()
	    };

    };



/** @} */



} // ---- end namespace genepi



#endif // ! __base_HiddenStateSpace_h
