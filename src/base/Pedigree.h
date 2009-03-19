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
/// \file Pedigree.h
/// Definition of the Pedigree class
//=============================================================================

#ifndef __base_Pedigree_h
#define __base_Pedigree_h


#include <cstddef>  // size_t

#include "config.h" // AGGRESSIVE_RANGE_CHECK

#include "Genotype.h"
#include "OrganismArray.h"
#include "SimpleLocusArray.h"

#include "AlleleArray.h" // AlleleProbTable, AlleleProbVect (for emission probabilities)

#include <vector>



namespace genepi { // ----



/** \addtogroup base
 * @{ */



// Needed for emission probability computation:
class AncestryVector;
class InheritanceVector;
class HiddenStateSpace;




//-----------------------------------------------------------------------------
//
/// Pedigree class.
///
/// Need more documentation here.
///
/// @warning
/// <SPAN STYLE="font-weight: bold; color: red;">IMPORTANT!</SPAN>: see
/// <A HREF="Pedigree_8cc.html#note-1"><B>NOTE *1*</B> in Pedigree.cc</A> regarding
/// copy constructors and assignment operators.
///
/// <A name="note-2"></A>
/// <TABLE STYLE="border: groove 3pt aqua;">
///
///  <TR>
///	<TD><B>NOTE *2*</B></TD>
///	<TD>
///	We implement an ordering of the Organism within the Pedigree for
///	iteration such that no member of a family will be visited prior to both
///	of its parents being visited.  This ordering could be achieved by
///	assigning each non-founder (or even founders too) a "depth" which will
///	be the <I>largest</I> number of steps required to traverse the
///	parent-tree starting at that node.  It is similar to "generation",
///	except that all founders of a given pedigree have the same depth of 0
///	(despite the fact that some may be of a different generation than
///	others), and a given individual's depth is then the <I>maximum</I>
///	path-length to a founder.  Another way to say this is that every sib's
///	depth is the maximum of its two parents' depths, plus 1.  Thus, all
///	founders have depth 0; a child both of whose parents are founders has
///	depth 1; and a child one of whose parents is a founder and the other of
///	whom is the child of two founders has depth 2.
///	<P>All of the iterative and array-style access to the container is via
///	this ordering.
///	<P>This ordering by "generational depth" has several desirable
///	properties: all founders are at the beginning of the list, so we can
///	easily define a "member index" within the pedigree and a "founder
///	index", both of which are the same number; and a "sib index" which is
///	the member index minus the number of founders; all of which are
///	contiguous 0-based ranges, useful for "parallel arrays," e.g. indexing
///	of the segregation indicators within the inheritance vector.
///	</TD>
///  </TR>
///
/// </TABLE>
///
/// <A name="note-3"></A>
/// <TABLE STYLE="border: groove 3pt aqua;">
///
///  <TR>
///	<TD><B>NOTE *3*</B></TD>
///	<TD>
///	The Pedigree does not "own" the Organism s within it (i.e. will not
///	delete them when it is deleted).  They remain controlled by the
///	OrganismArray, a reference to which the Pedigree keeps.
///	</TD>
///  </TR>
///
/// </TABLE>
//
//-----------------------------------------------------------------------------

class Pedigree
    {
    public:
	typedef Organism	 Member	    ;
	typedef Member * const * Iterator   ;
	typedef FamIdType	 IdType	    ; ///< Shorthand for FamIdType
	typedef size_t		 MemberIdx  ; ///< Index into sorted-organism-array for any member
	typedef size_t		 FounderIdx ; ///< Index into sorted-organism-array for founders
	typedef size_t		 SibIdx	    ; ///< Index into sorted-organism-array for non-founders

	// This is one way to receive the generated "states":
	typedef void (*StateReceiver)(	const Pedigree &	  ped		  ,
					size_t			  sLocIdx	  ,
					const AncestryVector &	  av		  ,
					const InheritanceVector & iv		  ,
					const Haplotype *	  founderHapState ,
					double			  emProb	  );


    private:
	const OrganismArray & memberPool;

	IdType	   id		 ;
	size_t	   nMembers	 ;
	size_t	   nFounders	 ;
	Member * * sortedMembers ;

	/// Array of hidden-state-spaces, one for each locus.  This belongs in a
	/// subclass.  We use a pointer-to-array rather than std::vector to
	/// avoid including the full class definition here.
	mutable HiddenStateSpace * stateProbs;


	void recurseSib(   SLocIdxType		sLocIdx		,
			   Haplotype *		memberHapState	,
			   const AncestryVector&ancestry	,
			   InheritanceVector &	iv		,
			   MemberIdx		memDepth	,
			   double		emProbTerm	,
			   StateReceiver	receiver	) const;

	void recurseFounder(	SLocIdxType		sLocIdx		,
				PopIdx			K		,
				const AlleleProbTable & alProbTab	,
				Haplotype *		memberHapState	,
				AncestryVector&		ancestry	,
				MemberIdx		memDepth	,
				double			probProdSoFar	,
				StateReceiver		receiver	) const;


    protected:
	Pedigree( const OrganismArray &	      pool   ,
		  OrganismArray::ConstPedIter firstM ,
		  OrganismArray::ConstPedIter endM   );


	// Range-checking:
	void throwFRange( size_t fIdx ) const;
	void throwMRange( size_t mIdx ) const;


	//-------------------------------------------------------
	// Generation of emission probabilities:
	//-------------------------------------------------------

	static void accumStateInArray(
			const Pedigree &	  ped		  ,
			size_t			  sLocIdx	  ,
			const AncestryVector &	  av		  ,
			const InheritanceVector & iv		  ,
			const Haplotype *	  founderHapState ,
			double			  emProb	  );


    public:
	/// Do not use: <A HREF="Pedigree_8cc.html#note-1"><B>NOTE *1*</B> in Pedigree.cc</A>
	Pedigree( const Pedigree & rhs );		// vector apparently requires this for push_back()
	Pedigree & operator=( const Pedigree & rhs );	// vector apparently requires this for push_back()
	~Pedigree();

	//---------------------------------------------------------------
	/// Create Pedigree's from raw genotype data (pedfile)
	/// The Pedigree objects are created in @a rv (output parameter)
	//---------------------------------------------------------------
	static void generatePedigrees( const OrganismArray &   organisms ,
				       std::vector<Pedigree> & rv	 );


	//---------------------------------------------------------------
	// Access to data members:
	//---------------------------------------------------------------

	/// Return the pool (pedfile) from which members are drawn
	const OrganismArray &	 getMemberPool() const { return memberPool; }
	const SimpleLocusArray & getSLoci     () const { return memberPool.getSLoci(); } ///< Convenience
	const IdType &		 getId	      () const { return id; }
	///< The "family ID" from the pedfile.  Every Organism in this Pedigree
	///< returns the same value of Organism::getFamId() (which is also returned
	///< here).


	// Iterative access:
	Iterator getFirstMember	 () const { return sortedMembers; }		// Docs below
	Iterator getEndMember	 () const { return (sortedMembers+nMembers); }	// Docs below; cache ptr?
	Iterator getFirstFounder () const { return sortedMembers; }		// Docs below
	Iterator getEndFounder	 () const { return (sortedMembers+nFounders); } // Docs below; cache ptr?
	Iterator getFirstNonFndr () const { return (sortedMembers+nFounders); }	// Docs below
	Iterator getEndNonFndr	 () const { return (sortedMembers+nMembers ); } // Docs below; cache ptr?

	size_t getNMembers () const { return nMembers ; }		///< Number of members
	size_t getNFounders() const { return nFounders; }		///< Number of founders
	size_t getNSibs	   () const { return (nMembers - nFounders); }	///< Number of non-founders


	/// Array-style access: @a mIdx must be between 0 and getNMembers().
	const Member & memberAt( MemberIdx mIdx ) const
	    {
	    #if AGGRESSIVE_RANGE_CHECK
		if ( mIdx >= getNMembers() )
		    throwMRange( mIdx );
	    #endif

	    return *(sortedMembers[mIdx]);
	    }


	/// Array-style access: @a fIdx must be between 0 and getNFounders().
	const Member & founderAt( FounderIdx fIdx ) const
	    {
	    #if AGGRESSIVE_RANGE_CHECK
		if ( fIdx >= getNFounders() )
		    throwFRange( fIdx );
	    #endif

	    return *(sortedMembers[fIdx]);
	    }



	//---------------------------------------------------------------
	// Generation of emission probabilities: these methods are implemented
	// in the PedigreeGenStates.cc
	//---------------------------------------------------------------

	// Documentation in source file:
	void genPossibleStates( StateReceiver receiver, PopIdx K, const AlleleProbTable & alProbTab, SLocIdxType sLocIdx ) const;
	void genPossibleStates( StateReceiver receiver, PopIdx K, const AlleleProbVect & alProbVect ) const;


	/// Same as genPossibleStates(), but accumulates the hidden states'
	/// emission probabilities in an internal structure that can be
	/// retrieved by getStateProbs().
	void genPossibleStatesInternal( PopIdx K, const AlleleProbVect & alProbVect ) const;


	/// Self-contained structure of states' probabilities, must be created
	/// using genPossibleStatesInternal() prior to calling.
	const HiddenStateSpace & getStateProbs( SLocIdxType sLocIdx ) const;


	// These are not guaranteed to produce any output if ostream support is
	// not compiled into Pedigree, InheritanceVector, AncestryVector, etc.:
	static void dbgRecursion( bool nv ); ///< Request that recursion debugging information output to cout
	static void dbgEmission ( bool nv ); ///< Request that emission-probability debugging output to cout
    };



/** @} */



// ============== ADDITIONAL DOCUMENTATION: ================
/**
 * \fn Pedigree::Iterator Pedigree::getFirstMember() const
 *	Returns an Iterator pointing to first member in
 *	<A HREF="#note-2"><I>generational-depth ordering</I></A>;
 *	use getEndMember() to test for end-of-range.
 *
 *\fn Pedigree::Iterator Pedigree::getEndMember() const
 *	Returns an Iterator pointing to the "last member plus one" (see getFirstMember()).
 *
 *\fn Pedigree::Iterator Pedigree::getFirstFounder() const
 *	Returns an Iterator pointing to first founder;
 *	use getEndFounder() to test for end-of-range.
 *
 *\fn Pedigree::Iterator Pedigree::getEndFounder() const
 *	Returns an Iterator pointing to the "last founder plus one" (see getFirstFounder()).
 *	This is equivalent to getFirstNonFndr()
 *
 *\fn Pedigree::Iterator Pedigree::getFirstNonFndr() const
 *	Returns an Iterator pointing to first non-founder;
 *	use getEndNonFndr() to test for end-of-range.
 *	This is equivalent to getEndFounder()
 *
 *\fn Pedigree::Iterator Pedigree::getEndNonFndr() const
 *	Returns an Iterator pointing to the "last non-founder plus one" (see getFirstNonFndr()).
 *	This is equivalent to getEndMember()
 *
 */



} // ---- end namespace genepi



#endif // ! __base_Pedigree_h
