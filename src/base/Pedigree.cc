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
///	The copy-constructors and assignment-operators of AlleleArray and
///	Pedigree only exist for the purposes of inserting in and/or moving
///	around in a container; they are public but should not be used to make
///	multiple copies of an object because they are <B>extremely</B>
///	dangerous: we presume that the right-hand-side is about to be destroyed;
///	rather than reference-count or allocate-and-copy, the left-hand-side
///	"takes ownership" of the dynamically-allocated objects, and we therefore
///	remove the RHS's references, so that it does not destroy (i.e.
///	<CODE>delete[]</CODE>) them when it is destroyed.
///	</TD>
///  </TR>
///
///  <TR>
///	<TD><B>NOTE *2*</B></TD>
///	<TD>
///	This is really two classes, the general-purpose Pedigree class and the
///	specialization for Admixmap, but combined into one.  These can and
///	should be separated.  Most of the code which belongs in the derived
///	class is in the separate source file AdmixPedigree.cc, but a few things
///	(e.g. constructor initialisation of data members) must be mixed with the
///	base-class code.  It is identified with the comment "NOTE *2*" to
///	facilitate later separation.
///	</TD>
///  </TR>
///
/// </TABLE>
//
//=============================================================================

#include "Pedigree.h"

#include <cstring>	// memcpy()

#include "SimpleLocus.h"
#include "SimpleLocusParser.h"
#include "InheritanceVector.h"
#include "HiddenStateSpace.h"

#include "HiddenMarkovModel.new.h"  // See NOTE *2*
#include "TransProbCache.h"	    // See NOTE *2*


#define USE_QSORT	1 ///< Whether to use ::qsort() or std::sort() for sorting arrays


// We really must find a way to auto-box with STL:
#if 0
    #define PED_PB(rv, organisms, firstInd, lastInd ) \
	rv.push_back( organisms, firstInd, lastInd ) , \
	rv.back().myNumber = rv.size()
#else
    #define PED_PB(rv, organisms, firstInd, lastInd ) \
	rv.push_back( Pedigree( organisms, firstInd, lastInd ) ) , \
	rv.back().myNumber = rv.size()
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

    int rv = lhs->getDepth() - rhs->getDepth();

    // Tiebreaker: organisms with genotyped data will be higher in the list:
    if ( rv == 0 )
	rv = int(rhs->isGenotyped()) - int(lhs->isGenotyped());

    #if 0 // ****** DEBUG: ******
	if ( lhs->getFamId() == "10" )
	    fprintf( stderr, " compare: %p(%d-%s) vs. %p(%d-%s): %d\n",
		lhs, lhs->getDepth(), lhs->isGenotyped() ? "gtype" : "missing",
		rhs, rhs->getDepth(), rhs->isGenotyped() ? "gtype" : "missing", rv );
    #endif

    return rv;
    }



//-----------------------------------------------------------------------------
// Constructor [protected]
//-----------------------------------------------------------------------------

Pedigree::Pedigree( const OrganismArray & pool,
		    OrganismArray::ConstPedIter firstM, OrganismArray::ConstPedIter endM ) :
	memberPool	( pool			   ) ,
	id		( (*firstM)->getFamId()	   ) ,
	nMembers	( endM - firstM		   ) ,
	sortedMembers	( new Member* [ nMembers ] ) ,
	nMendelErrs	( 0			   ) ,
	mendelErrsByLocus( 0			   ) ,
	stateProbs	( 0			   ) ,
	nAffected	( 0			   ) ,
	step		( 0.3			   ) ,
	NumberOfUpdates ( 0			   ) ,
	w		( 1			   ) ,
	NumGametes	( 2			   ) ,
	tpCache		( 0			   ) ,
	hmm		( 0			   )
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

	if ( m.getOutcome() != 0 )
	    ++nAffected;
#if 0
	else
	    {
	    // If requested, exclude unaffected siblings.  This should probably take
	    // place earlier, in GenotypeParser.
	    if ( getOptions().getExcludeUnaffectedSibs() && (! m.isFounder()) && m.getChildren().empty() )
		{
		++firstM;
		continue;
		}
	    }
#endif

	++firstM;
	++mPtr;
	}


    // See NOTE *2*:
    // We can't allocate these prior to knowing the number of founders.  Of
    // course they should be in a subclass, in which case they wouldn't be
    // allocated until much later anyhow.
    bclib::pvector<double> Kzero;
    Kzero.resize( getK(), 0.0 );
    Theta	    .resize( getNTheta()	);
    ThetaProposal   .resize( getNTheta()	);
    SumSoftmaxTheta .resize( getNTheta(), Kzero );
    thetahat	    .resize( getNTheta(), Kzero );
    dirparams	    .resize( getK()		);


    // Set the initial values for Theta:
    SetUniformAdmixtureProps();


    #if 0 // ****** DEBUG: ******
	if ( id == "10" )
	    {
	    fprintf( stderr, "ID: %s  N-members: %lu  N-founders; %lu\n", id.c_str(), getNFounders(), getNMembers() );
	    for ( size_t idx = 0 ; idx < nMembers ; ++idx )
		fprintf( stderr, "  %lu %p %hd %s\n", idx, sortedMembers[idx], sortedMembers[idx]->getDepth(),
		    sortedMembers[idx]->isFounder() ? "founder" : "non-founder" );
	    }
    #endif

    // Once we've traversed every founder's child-tree, we should have the
    // correct depth for every member of the pedigree since there should be no
    // disjoint connected sub-graphs per the check in OrganismArray.

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
    gp_assert_le( getNNonFndrs(), InheritanceVector::MAX_ORGANISMS );
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
	memberPool	    ( rhs.memberPool	    ) ,
	id		    ( rhs.id		    ) ,
	nMembers	    ( rhs.nMembers	    ) ,
	nFounders	    ( rhs.nFounders	    ) ,
	sortedMembers	    ( rhs.sortedMembers	    ) ,
	nMendelErrs	    ( rhs.nMendelErrs	    ) ,
	mendelErrsByLocus   ( rhs.mendelErrsByLocus ) ,
	stateProbs	    ( rhs.stateProbs	    ) ,
	nAffected	    ( rhs.nAffected	    ) ,
	Theta		    ( rhs.Theta		    ) , // See NOTE *2*
	ThetaProposal	    ( rhs.ThetaProposal	    ) , // See NOTE *2*
	SumSoftmaxTheta	    ( rhs.SumSoftmaxTheta   ) , // See NOTE *2*
	thetahat	    ( rhs.thetahat	    ) , // See NOTE *2*
	dirparams	    ( rhs.dirparams	    ) , // See NOTE *2*
	_rho		    ( rhs._rho		    ) , // See NOTE *2*
	sumlogrho	    ( rhs.sumlogrho	    ) , // See NOTE *2*
	rhohat		    ( rhs.rhohat	    ) , // See NOTE *2*
	logLikelihood	    ( rhs.logLikelihood	    ) , // See NOTE *2*
	step		    ( rhs.step		    ) , // See NOTE *2*
	NumberOfUpdates	    ( rhs.NumberOfUpdates   ) , // See NOTE *2*
	w		    ( rhs.w		    ) , // See NOTE *2*
	ThetaTuner	    ( rhs.ThetaTuner	    ) , // See NOTE *2*
	NumGametes	    ( rhs.NumGametes	    ) , // See NOTE *2*
	myNumber	    ( rhs.myNumber	    ) , // See NOTE *2*
	tpCache		    ( rhs.tpCache	    ) , // See NOTE *2*
	hmm		    ( rhs.hmm		    )	// See NOTE *2*
    {
    // !!!WARNING!!! -- see NOTE *1*
    const_cast<Pedigree&>(rhs).sortedMembers	 = 0;
    const_cast<Pedigree&>(rhs).mendelErrsByLocus = 0;
    const_cast<Pedigree&>(rhs).stateProbs	 = 0;
    const_cast<Pedigree&>(rhs).tpCache		 = 0; // See NOTE *2*
    const_cast<Pedigree&>(rhs).hmm		 = 0; // See NOTE *2*
    }

Pedigree & Pedigree::operator=( const Pedigree & rhs )
    {
    gp_assert( &memberPool == &rhs.memberPool );

    id			= rhs.id		;
    nMembers		= rhs.nMembers		;
    nFounders		= rhs.nFounders		;
    sortedMembers	= rhs.sortedMembers	;
    nMendelErrs		= rhs.nMendelErrs	;
    mendelErrsByLocus	= rhs.mendelErrsByLocus ;
    stateProbs		= rhs.stateProbs	;
    nAffected		= rhs.nAffected		;

    // See NOTE *2*:
    Theta		= rhs.Theta		;
    ThetaProposal	= rhs.ThetaProposal	;
    SumSoftmaxTheta	= rhs.SumSoftmaxTheta	;
    thetahat		= rhs.thetahat		;
    dirparams		= rhs.dirparams		;
    _rho		= rhs._rho		;
    sumlogrho		= rhs.sumlogrho		;
    rhohat		= rhs.rhohat		;
    logLikelihood	= rhs.logLikelihood	;
    step		= rhs.step		;
    NumberOfUpdates	= rhs.NumberOfUpdates	;
    w			= rhs.w			;
    ThetaTuner		= rhs.ThetaTuner	;
    NumGametes		= rhs.NumGametes	;
    myNumber		= rhs.myNumber		;


    tpCache		= rhs.tpCache		;
    hmm			= rhs.hmm		;

    // !!!WARNING!!! -- see NOTE *1*
    const_cast<Pedigree&>(rhs).sortedMembers	 = 0;
    const_cast<Pedigree&>(rhs).mendelErrsByLocus = 0;
    const_cast<Pedigree&>(rhs).stateProbs	 = 0;
    const_cast<Pedigree&>(rhs).tpCache		 = 0; // See NOTE *2*
    const_cast<Pedigree&>(rhs).hmm		 = 0; // See NOTE *2*

    return *this;
    }



//-----------------------------------------------------------------------------
// Destructor
//-----------------------------------------------------------------------------

Pedigree::~Pedigree()
    {
    delete[] sortedMembers;
    delete[] stateProbs;
    delete[] mendelErrsByLocus;
    delete tpCache ; // See NOTE *2*
    delete hmm	   ; // See NOTE *2*
    }



//---------------------------------------------------------------
// generatePedigrees() [static]
//
/// Create Pedigree's from raw genotype data (pedfile)
//---------------------------------------------------------------

void Pedigree::generatePedigrees( const OrganismArray & organisms, cvector<Pedigree> & rv )
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

    // Ugly data-type but pretty method-call (no "output"/"reference" parameters):
    const std::pair<OrganismArray::ConstPedIter,OrganismArray::ConstPedIter> &
	interval = organisms.findOrgsInPed( curFamId );

    PED_PB( rv, organisms, interval.first, interval.second );
    }



//-----------------------------------------------------------------------------
// throwFRange(), throwMRange() [protected]
//-----------------------------------------------------------------------------

void Pedigree::throwFRange( size_t fIdx ) const
    {
    throw std::out_of_range( estr("Founder-index ") + fIdx +
		" out of range (" + getNFounders() + ')' );
    }

void Pedigree::throwMRange( size_t mIdx ) const
    {
    throw std::out_of_range( estr("Member-index ") + mIdx +
		" out of range (" + getNMembers() + ')' );
    }



//-----------------------------------------------------------------------------
// haveMendelErrAt()
//-----------------------------------------------------------------------------

bool Pedigree::haveMendelErrAt( SLocIdxType t ) const
    {
    #if AGGRESSIVE_RANGE_CHECK
	if ( t >= getNSLoci() )
	    throw std::out_of_range( estr("Simple-locus-index ") + t +
		    " out of range (" + getNSLoci() + ')' );
    #endif

    return (nMendelErrs != 0) && mendelErrsByLocus[t];
    }



} // ---- end namespace genepi
