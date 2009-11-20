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
/// \file ProposalRingBuffer.h
/// Definition of the ProposalRingBuffer class.
//=============================================================================


#ifndef __base_ProposalRingBuffer_h
#define __base_ProposalRingBuffer_h



namespace genepi { // ----

/** \addtogroup base
 * @{ */



//-----------------------------------------------------------------------------
/// A simple ring buffer to help manage proposing, accepting, and rejecting new
/// parameters and the resulting computions.  Used internally by Pedigree.  As
/// currently implemented, the ring-buffer has only 2 elements, so there can be
/// at most 1 active proposal.
//-----------------------------------------------------------------------------

template < typename T > class ProposalRingBuffer
    {
    private:
	static const size_t MAX_VALS = 2;

	T values[ MAX_VALS ];

	T * curVal  ; ///< proposed if a proposal is active, last accepted otherwise
	T * accepted; ///< last accepted

    public:

	ProposalRingBuffer() :
		curVal	( values ) ,
		accepted( values ) {}

	/// Copy constructor; copies the contents (both accepted and proposed)
	/// from the rhs.  In practice, these only get copied when they are not
	/// yet used, so this could be blank-initialized rather than copied.
	ProposalRingBuffer( const ProposalRingBuffer & rhs ) :
		curVal	( values + (rhs.curVal	 - rhs.values) ) ,
		accepted( values + (rhs.accepted - rhs.values) )
	    {
	    for ( size_t idx = MAX_VALS ; idx-- != 0 ; )
		values[idx] = rhs.values[idx];
	    }


	/// Assignment operator; copies the contents (both accepted and
	/// proposed) from the rhs.  In practice, these only get copied when
	/// they are not yet used, so this could be blank-initialized rather
	/// than copied.
	ProposalRingBuffer & operator=( const ProposalRingBuffer & rhs )
	    {
	    for ( size_t idx = MAX_VALS ; idx-- != 0 ; )
		values[idx] = rhs.values[idx];

	    curVal   = values + (rhs.curVal   - rhs.values);
	    accepted = values + (rhs.accepted - rhs.values);

	    return *this;
	    }


	/// "Allocates" a new slot on the ring-buffer for proposal and returns
	/// it.  The contents of the new slot are left indeterminate, and in
	/// reality will often be a previous proposal.
	T & startNewProposal()
	    {
	    gp_assert( ! proposalActive() );
	    curVal = (accepted == values) ? (values + 1) : values;
	    return *curVal;
	    }

	/// Is a proposal currently underway?
	bool proposalActive() const { return (curVal != accepted); }

	/// Return the proposed slot if a proposal is currently being
	/// considered; otherwise the last accepted slot (non-const version).
	T &	  getCurrent()	     { return *curVal; }

	/// Return the proposed slot if a proposal is currently being
	/// considered; otherwise the last accepted slot (const version).
	const T & getCurrent() const { return *curVal; }

	/// Assert that a proposal is currently being considered and
	/// return the proposed value.
	T &	  getProposed()	      { gp_assert(proposalActive()); return *curVal; }
	const T & getProposed() const { gp_assert(proposalActive()); return *curVal; }

	/// Return the last accepted value (non-const version).
	T &	  getAccepted()	      { return *accepted; }
	/// Return the last accepted value (const version).
	const T & getAccepted() const { return *accepted; }

	/// Accept a proposal (one must be active): makes previously-alloc'd
	/// proposal the new accepted-slot; "frees" the old accepted-slot (in
	/// reality, just keeps it for the next proposal).
	void acceptProposal()
	    {
	    gp_assert( proposalActive() );
	    accepted = curVal;
	    }

	/// Reject a proposal (one must be active): "frees" the
	/// previously-allocated proposal (in reality, just keeps it for the
	/// next proposal); the current accepted slot remains accepted.
	void rejectProposal()
	    {
	    gp_assert( proposalActive() );
	    curVal = accepted;
	    }


	/// "Iterator" style access to the raw array, might be used for
	/// initialisation, otherwise should not ordinarily be used.
	typedef T * Iterator;
	Iterator begin() { return values	     ; }
	Iterator end  () { return (values + MAX_VALS); }
    };



/** @} */

} // ---- end namespace genepi



#endif // ! __base_ProposalRingBuffer_h
