//=============================================================================
//
// Copyright (C) 2009  David D. Favro  gpl@meta-dynamic.com
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
//
/// \file Pedigree.cc
/// Implementation of the Pedigree class (except state-generation is in PedigreeGenStates.cc).
///
/// <A name="note-1"></A>
/// <TABLE STYLE="border: groove 3pt aqua;">
///
///  <TR>
///	<TD><B>NOTE *1*</B></TD>
///	<TD>
///	The copy-constructors and assignment-operators of State and Pedigree
///	only exist for the purposes of inserting in and/or moving around in a
///	container; they are public but should not be used to make multiple
///	copies of an object because they are <B>extremely</B> dangerous: we
///	presume that the right-hand-side is about to be destroyed; rather than
///	reference-count or allocate-and-copy, the left-hand-side "takes
///	ownership" of the dynamically-allocated objects, and we therefore remove
///	the RHS's references, so that it does not destroy (i.e.
///	<CODE>delete[]</CODE>) them when it is destroyed.
///	</TD>
///  </TR>
///
/// </TABLE>
//
//=============================================================================

#include "Pedigree.h"

#include <cstring>	// memcpy()

#include "InheritanceVector.h"
#include "SimpleLocus.h"
#include "SimpleLocusParser.h"


#define USE_QSORT	1
#define DEBUG_GH_COMP	0


#if DEBUG_GH_COMP // ****** DEBUG: ******
    #include <iostream>
#endif


// We really must find a way to auto-box with STL:
#if 0
    #define PED_PB(rv, organisms, firstInd, lastInd ) \
	rv.push_back( organisms, firstInd, lastInd )
#else
    #define PED_PB(rv, organisms, firstInd, lastInd ) \
	rv.push_back( Pedigree( organisms, firstInd, lastInd ) )
#endif



// Using qsort() rather than std::sort() because it wasn't working for me and I
// can't be bothered to figure out why.
#if USE_QSORT
    #include <cstdlib> // qsort()
    template < typename T >
	void gqsort( T * array, size_t n_el, int(*compFunc)(const T *, const T *) )
	    {
	    ::qsort( array, n_el, sizeof(T), (int(*)(const void *, const void *))compFunc );
	    }
    #define SORT(array,n_el,compFunc) gqsort(array,n_el,compFunc)
#else
    #include <algorithm> // std::sort()
    #define SORT(array,n_el,compFunc) std::sort(array,array+n_el,compFunc)
#endif



using namespace std;



namespace genepi { // ----



State::State( const State & rhs ) :
	iv( rhs.iv )
    {
    #if 0 // See NOTE *1*
	founderHaps = new Haplotype[ iv.getPedigree().getNFounders() ];

	// See NOTE *3* in Genotype.h regarding safety of memcpy():
	memcpy( founderHaps, rhs.founderHaps, iv.getPedigree().getNFounders() * sizeof(*founderHaps) );
    #else
	// !!!WARNING!!! -- see NOTE *1*
	founderHaps = rhs.founderHaps;
	const_cast<State&>(rhs).founderHaps = 0;
    #endif
    }


State::State( const InheritanceVector & _iv, const Haplotype * _founderHaps ) :
	iv( _iv )
    {
    founderHaps = new Haplotype[ iv.getNFounders() ];

    // See NOTE *3* in Genotype.h regarding safety of memcpy():
    memcpy( founderHaps, _founderHaps, getIV().getNFounders() * sizeof(*founderHaps) );
    }


State::~State()
    {
    delete[] founderHaps;
    }


State & State::operator=( const State & rhs )
    {
    iv = rhs.iv;

    #if 0 // See NOTE *1*
	delete[] founderHaps;
	founderHaps = new Haplotype[ iv.getPedigree().getNFounders() ];

	// See NOTE *3* in Genotype.h regarding safety of memcpy():
	memcpy( founderHaps, rhs.founderHaps, iv.getPedigree().getNFounders() * sizeof(*founderHaps) );
    #else
	// !!!WARNING!!! -- see NOTE *1*
	founderHaps = rhs.founderHaps;
	const_cast<State&>(rhs).founderHaps = 0;
    #endif

    return *this;
    }



void State::throwRange( size_t fIdx ) const
    {
    throw std::runtime_error( estr("Founder-index ") + fIdx +
		" out of range (" + getIV().getNFounders() + ')' );
    }



//-----------------------------------------------------------------------------
// traverseChildTree()
//
/// Static helper function (recursive).  Used to compute the "depth" for each node.
//-----------------------------------------------------------------------------

static void traverseChildTree( Organism & ind )
    {
    typedef std::list<Organism*> ChList;


    const int depth = ind.getDepth() + 1;

    // Traverse the children (depth-first):
    const ChList & children = ind.getChildren();
    const ChList::const_iterator limit = children.end();
    for ( ChList::const_iterator iter = children.begin() ; iter != limit ; ++iter )
	{
	Organism & child = **iter;

	// Has been visited already:
	const bool visited = (child.getDepth() != Organism::UNKNOWN_DEPTH);

	// We avoid re-traversing sub-trees that we've already seen (e.g. via
	// the other parent), provided that the depth is the same as or less
	// than the last time; if the depth is now greater, we must fix this
	// node and the sub-tree below it:
	if ( (! visited) || (child.getDepth() < depth) )
	    {
	    child.setDepth( depth );
	    traverseChildTree( child );
	    }
	}
    }



//-----------------------------------------------------------------------------
// depthCompare()
//
/// Comparison function for sort(), sorts the members of a pedigree by depth
/// (see <B>NOTE *1*</B> in OrganismArray.h).
//-----------------------------------------------------------------------------

static int depthCompare( Pedigree::Member * const * lhs_ptr, Pedigree::Member * const * rhs_ptr )
    {
    Pedigree::Member * const lhs = *lhs_ptr;
    Pedigree::Member * const rhs = *rhs_ptr;

    #if 0 // ****** DEBUG: ******
	if ( lhs->getFamId() == "10" )
	    fprintf( stderr, " compare: %p(%d) vs. %p(%d): %d\n", lhs, lhs->getDepth(),
				rhs, rhs->getDepth(), (lhs->getDepth() - rhs->getDepth()) );
    #endif
    return (lhs->getDepth() - rhs->getDepth());
    }



//-----------------------------------------------------------------------------
// Constructor [protected; use generatePedigrees()]
//-----------------------------------------------------------------------------

Pedigree::Pedigree( const OrganismArray & pool,
		    OrganismArray::ConstPedIter firstM, OrganismArray::ConstPedIter endM ) :
	memberPool	( pool			   ) ,
	id		( (*firstM)->getFamId()	   ) ,
	nMembers	( endM - firstM		   ) ,
	sortedMembers	( new Member* [ nMembers ] ) ,
	consistentStates( new vector<State>[ pool.getNSimpleLoci() ] )
    {
    // Initialize the array of pointers-to-members, while simultaneously
    // counting the number of founders and traversing the parent-tree to compute
    // the depth value (maximum depth-to-founder) for each member:
    nFounders = 0;
    Member * * mPtr = sortedMembers;
    while ( firstM != endM )
	{
	*mPtr = *firstM;
	Member & m = **firstM;
	if ( m.isFounder() )
	    {
	    gp_assert( m.getDepth() == 0 );
	    ++nFounders;
	    traverseChildTree( m );
	    }
	++firstM;
	++mPtr;
	}

    // Once we've traversed every founder's child-tree, we should have the
    // correct depth for every member of the pedigree since there should be no
    // disjoint connected sub-graphs per the check in OrganismArray.

    #if 0 // ****** DEBUG: ******
	if ( id == "10" )
	    {
	    fprintf( stderr, "ID: %s  N-members: %lu  N-founders; %lu\n", id.c_str(), getNFounders(), getNMembers() );
	    for ( size_t idx = 0 ; idx < nMembers ; ++idx )
		fprintf( stderr, "  %lu %p %hd %s\n", idx, sortedMembers[idx], sortedMembers[idx]->getDepth(),
		    sortedMembers[idx]->isFounder() ? "founder" : "non-founder" );
	    }
    #endif

    SORT( sortedMembers, getNMembers(), depthCompare );

    #if 0 // ****** DEBUG: ******
	if ( id == "10" )
	    {
	    fprintf( stderr, "ID: %s  N-members: %lu  N-founders; %lu\n", id.c_str(), getNFounders(), getNMembers() );
	    for ( size_t idx = 0 ; idx < nMembers ; ++idx )
		fprintf( stderr, "  %lu %p %hd %s\n", idx, sortedMembers[idx], sortedMembers[idx]->getDepth(),
		    sortedMembers[idx]->isFounder() ? "founder" : "non-founder" );
	    }
    #endif


    // At this point, the depth field is effectively no longer needed; we now
    // re-use it as the index within the sorted-list (see NOTE *2* in
    // Organism.h):
    for ( int idx = getNMembers() ; idx-- != 0 ; )
	{
	// Assure that the child-tree-traversals of the founders above did
	// indeed reach every member of the pedigree (this should have been
	// guaranteed by the connected-test in GenotypeParser):
	gp_assert( sortedMembers[idx]->getDepth() != Member::UNKNOWN_DEPTH );
	sortedMembers[idx]->setDepth( idx );
	}


    // Cache pointers to the first and ending founders, which should all be at
    // the start of the sortedMembers array:
    #if 0
	firstFounder = sortedMembers;
	endFounder = sortedMembers + nFounders;
    #endif

    // Sanity checking:
    gp_assert( nFounders != 0 );
    gp_assert_le( nFounders, nMembers );
    gp_assert_le( getNSibs(), InheritanceVector::MAX_ORGANISMS );
    gp_assert( sortedMembers[0]->isFounder() );
    gp_assert( getEndFounder()[-1]->isFounder() );

    if ( getEndFounder() != getEndMember() ) // If we have any non-founders at all
	gp_assert( ! (*getEndFounder())->isFounder() );
    }



//-----------------------------------------------------------------------------
// Copy constructor and assignment operator: these are required for
// std::vector<Pedigree> and I believe are safe if never used from anywhere but
// there; they are public, but should never be used except for immediate
// destruction of the copied object.
//-----------------------------------------------------------------------------

Pedigree::Pedigree( const Pedigree & rhs ) :
	memberPool	( rhs.memberPool	) ,
	id		( rhs.id		) ,
	nMembers	( rhs.nMembers		) ,
	nFounders	( rhs.nFounders		) ,
	sortedMembers	( rhs.sortedMembers	) ,
	consistentStates( rhs.consistentStates	)
    {
    // !!!WARNING!!! -- see NOTE *1*
    const_cast<Pedigree&>(rhs).sortedMembers	= 0;
    const_cast<Pedigree&>(rhs).consistentStates = 0;
    }

Pedigree & Pedigree::operator=( const Pedigree & rhs )
    {
    gp_assert( &memberPool == &rhs.memberPool );

    id		     = rhs.id		    ;
    nMembers	     = rhs.nMembers	    ;
    nFounders	     = rhs.nFounders	    ;
    sortedMembers    = rhs.sortedMembers    ;
    consistentStates = rhs.consistentStates ;

    // !!!WARNING!!! -- see NOTE *1*
    const_cast<Pedigree&>(rhs).sortedMembers	= 0;
    const_cast<Pedigree&>(rhs).consistentStates = 0;

    return *this;
    }



//-----------------------------------------------------------------------------
// Destructor
//-----------------------------------------------------------------------------

Pedigree::~Pedigree()
    {
    delete[] sortedMembers;
    delete[] consistentStates;
    }



//---------------------------------------------------------------
// generatePedigrees() [static]
//
// Create Pedigree's from raw genotype data (pedfile)
//---------------------------------------------------------------

void Pedigree::generatePedigrees( const OrganismArray & organisms, vector<Pedigree> & rv )
    {
    FamIdType curFamId;

    const OrganismArray::ConstPedIter limit = organisms.endByPed();
    for ( OrganismArray::ConstPedIter iter = organisms.beginByPed(); iter != limit; ++iter )
	{
	const Organism & ind = **iter;
	if ( ind.getFamId() != curFamId )
	    {
	    if ( ! curFamId.empty() )
		{
		OrganismArray::ConstPedIter firstInd;
		OrganismArray::ConstPedIter lastInd;
		organisms.findOrgsInPed( curFamId, firstInd, lastInd );

		PED_PB( rv, organisms, firstInd, lastInd );
		}

	    curFamId = ind.getFamId();
	    }
	}

    gp_assert( ! curFamId.empty() );

    // Ugly data-type but pretty method-call:
    const std::pair<OrganismArray::ConstPedIter,OrganismArray::ConstPedIter> &
	interval = organisms.findOrgsInPed( curFamId );

    PED_PB( rv, organisms, interval.first, interval.second );
    }



//-----------------------------------------------------------------------------
// throwFRange(), throwMRange() [protected]
//-----------------------------------------------------------------------------

void Pedigree::throwFRange( size_t fIdx ) const
    {
    throw std::runtime_error( estr("Founder-index ") + fIdx +
		" out of range (" + getNFounders() + ')' );
    }

void Pedigree::throwMRange( size_t mIdx ) const
    {
    throw std::runtime_error( estr("Member-index ") + mIdx +
		" out of range (" + getNMembers() + ')' );
    }



} // ---- end namespace genepi
