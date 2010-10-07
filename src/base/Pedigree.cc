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

#include "SimpleLocus.h"
#include "SimpleLocusParser.h"
#include "InheritanceVector.h"
#include "HiddenStateSpace.h"

#include "HiddenMarkovModel.new.h"	// See NOTE *2*
#include "TransProbCache.h"		// See NOTE *2*
#include "../admixmap/AdmixOptions.h"	// See NOTE *2* (for --exclude-unaffected-sibs)


#define USE_QSORT	1 ///< Whether to use ::qsort() or std::sort() for sorting arrays



#if PED_USE_EMPLACE
  //#warning Pedigree is using c++0x vector::emplace() (and has a public constructor)
  #define PED_PB(rv, organisms, firstInd, lastInd, max_n ) \
	rv.emplace_back( organisms, firstInd, lastInd, max_n ) , \
	rv.back().setMyNumber( rv.size() );
#else
  #define PED_PB(rv, organisms, firstInd, lastInd, max_n ) \
	rv.push_back( Pedigree( organisms, firstInd, lastInd, max_n ) ) , \
	rv.back().setMyNumber( rv.size() );
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

typedef std::list<genepi::Organism*> ChList;



namespace genepi { // ----



//-----------------------------------------------------------------------------
// traverseChildTree()
//
/// Static helper function (recursive).  Used to compute the "depth" for each node.
//-----------------------------------------------------------------------------

static void traverseChildTree( Organism & ind )
    {

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

    // Tiebreaker: organisms with genotyped data will be higher in the list.
    // Since we value affected individuals more highly than non-affected, we
    // want the --max-pedigree-size option will drop unaffected before affected
    // individuals, so if we still have a tie, we sort affected individuals
    // before unaffected.
    if ( rv == 0 )
	{
	rv = int(rhs->isGenotyped()) - int(lhs->isGenotyped());
	if ( rv == 0 )
	    {
	    const bool lhs_affected = (lhs->getOutcome() == Organism::OUTCOME_AFFECTED);
	    const bool rhs_affected = (rhs->getOutcome() == Organism::OUTCOME_AFFECTED);
	    rv = int(rhs_affected) - int(lhs_affected);
	    }
	}

    return rv;
    }



//-----------------------------------------------------------------------------
// dumpFGMappings()
//-----------------------------------------------------------------------------

#if DEBUG_PRINT_FOUNDER_GAMETE_MAPS

    void Pedigree::dumpFGMappings() const
	{

	fprintf( stderr, "\n==== Founder-to-gamete map for pedigree %s at %p (non-X chromosomes)\n",
			getId().c_str(), this );
	for ( FounderIdx fIdx = 0 ; fIdx < getNFounders() ; ++fIdx )
	    {
	    const Organism & founder = founderAt( fIdx );
	    if ( founder.isHaploid( CHR_IS_NOT_X /*FIXME-PED-XCHR*/ ) )
		fprintf( stderr, "    fndr-idx %02zd (org #%d) ->   single-gamete %zd\n",
		    fIdx, founder.getOrgId(), founderGameteOfFounder( fIdx, GT_SINGLE, onX ) );
	    else
		fprintf( stderr, "    fndr %02zd (org #%d) -> paternal-gamete %zd\n"
				 "            (org #%d) -> maternal-gamete %zd\n",
		    fIdx,
		    founder.getOrgId(), founderGameteOfFounder( fIdx, GT_PATERNAL, onX ),
		    founder.getOrgId(), founderGameteOfFounder( fIdx, GT_MATERNAL, onX ) );
	    }

	fprintf( stderr, "\n  == Gamete-to-founder map:\n" );
	GameteType whichOne;
	for ( FGameteIdx fgIdx = 0 ; fgIdx < getNFounderGametes() ; ++fgIdx )
	    {
	    const FounderIdx fIdx = founderOfGameteIdx( fgIdx, whichOne );
	    fprintf( stderr, "    gamete-idx %zd -> fndr %zd (org #%d) [%s]\n", fgIdx,
		fIdx, founderAt(fIdx).getOrgId(), gameteTypeDesc(whichOne) );
	    }

	putc( '\n', stderr );

	}

#endif



//-----------------------------------------------------------------------------
// Constructor [protected]
//-----------------------------------------------------------------------------

Pedigree::Pedigree( const OrganismArray &	pool	,
		    OrganismArray::ConstPedIter firstM	,
		    OrganismArray::ConstPedIter endM	,
		    MemberIdx			max_n	) :
	memberPool	  ( pool		     ) ,
	id		  ( (*firstM)->getFamId()    ) ,
	nMembers	  ( endM - firstM	     ) ,
	sortedMembers	  ( new Member* [ nMembers ] ) ,
	nMendelErrs	  ( 0			     ) ,
	mendelErrsByLocus ( 0			     ) ,
	stateProbs	  ( 0			     ) ,
	nAffected	  ( 0			     ) ,
	nAffNonFndr	  ( 0			     ) ,
	step		  ( 0.3			     ) ,
	NumberOfUpdates	  ( 0			     ) ,
	w		  ( 1			     ) ,
	NumGametes	  ( 2			     ) ,
	tpCache		  ( 0			     ) ,
	hmm		  ( 0			     ) ,
	llCache		  ( *this		     ) ,
	aoCache		  ( 0			     )
    {

    // Initialize the array of pointers-to-members, while simultaneously
    // counting the number of founders and traversing the parent-tree to compute
    // the depth value (maximum depth-to-founder) for each member:
    nFounders = 0;
    FGameteIdx nFGametes  = 0; // Counter for number of founder-gametes, non-X chroms
    FGameteIdx nFGametesX = 0; // Counter for number of founder-gametes, on the X chrom
    Member * * mPtr = sortedMembers;
    while ( firstM != endM )
	{
	*mPtr = *firstM;
	Member & m = **firstM;
	if ( m.isFounder() )
	    {
	    gp_assert( m.getDepth() == 0 );
	    ++nFounders;
	    nFGametes  += m.isHaploid( CHR_IS_NOT_X ) ? 1 : 2;
	    nFGametesX += m.isHaploid( CHR_IS_X	    ) ? 1 : 2;
	    traverseChildTree( m );
	    }

	if ( m.getOutcome() != Organism::OUTCOME_AFFECTED )
	    {
	    // If requested, exclude unaffected siblings.  This should probably take
	    // place earlier, in GenotypeParser.
	    if ( getOptions().getExcludeUnaffectedSibs() && (! m.isFounder()) && m.getChildren().empty() )
		{
		if ( getOptions().getDisplayLevel() >= 3 )
		    cout << "Pedigree " << getId() <<
			": excluding unaffected sibling " << m.getOrgId() << '\n';
		++firstM;
		--nMembers;
		continue;
		}
	    }

	++firstM;
	++mPtr;
	}


    // Once we've traversed every founder's child-tree, we should have the
    // correct depth for every member of the pedigree since there should be no
    // disjoint connected sub-graphs per the check in OrganismArray.

    SORT( sortedMembers, getNMembers(), depthCompare );



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


    // If max_n is specified and this pedigree's size exceeds it, remove enough
    // non-founders to reduce the size to max_n:
    if ( (max_n != 0) && (getNMembers() > max_n) )
	{

	if ( getNFounders() > max_n )
	    throw runtime_error( estr("Pedigree ") + getId() + " has " + getNFounders() +
			" founders, which exceeds the maximum number of members (" +
			max_n + ')' );

	if ( getOptions().getDisplayLevel() >= 3 )
	    {
	    cout << "Pedigree " << getId() << " size reduced from " << getNMembers() << " to "
		<< max_n << " per parameters (excluding organism IDs:";
	    for ( Member * * mPtr = sortedMembers + getNMembers() ; mPtr-- != (sortedMembers + max_n) ; )
		cout << ' ' << (*mPtr)->getOrgId();
	    cout << ").\n";
	    }

	nMembers = max_n;

	// We don't bother to re-allocate sortedMembers; just leave the
	// abandoned ones hanging on the end of the too-long array.

	// "Traverse" the graph, removing child-pointers to the removed
	// non-founders.  This is super-hack.
	for ( Member * * mPtr = sortedMembers + getNMembers() ; mPtr-- != sortedMembers ; )
	    {
	    ChList & children = (*mPtr)->getChildren();
	    for ( ChList::iterator chIter = children.begin() ; chIter != children.end() ; ++chIter )
		while ( ((*chIter)->getPIdx() >= max_n) &&
			((chIter = children.erase(chIter)) == children.end()) )
		    ; // Do nothing
	    }

	}


    // ----- SINGLE-GAMETE-FOUNDER MODEL -----
    // Since organisms now know their children, we can reduce certain founders
    // to be modeled as a single gamete: those without genotyped data and only a
    // single offspring in the pedigree.
    if ( getOptions().getSglGameteFounder() )
	for ( Iterator it = getFirstFounder() ; it != getEndFounder() ; ++it )
	    {

	    Member & f = **it;

	    if ( (! f.isGenotyped()) && (f.getChildren().size() == 1) &&
		    (! f.isHaploid( CHR_IS_NOT_X )) )
		{
		if ( getOptions().getDisplayLevel() > 3 )
		    cout << "Pedigree " << getId()
			<< ": using single-gamete model for founder "
			<< f.getOrgId();
		--nFGametes;
		if ( ! f.isHaploid( CHR_IS_X ) )
		    {
		    if ( getOptions().getDisplayLevel() > 3 )
			cout << " (both X & non-X chromosomes).\n";
		    --nFGametesX;
		    }
		else
		    {
		    if ( getOptions().getDisplayLevel() > 3 )
			cout << " (X chromosome was already haploid).\n";
		    }
		f.setSingleGameteModel( true );
		}

	    }



    // Build the founder-index to founder-gamete-index and founder-gamete-index
    // to founder-index mappings:
    fgMap.resize( nFounders );
    gfMap.resize( nFGametes );
    gfMapX.resize( nFGametesX );
    FGMapEl fg = { 0, 0 };
    for ( FounderIdx fIdx = 0 ; fIdx < nFounders ; ++fIdx )
	{

	fgMap[ fIdx ] = fg;
	const Organism & founder = founderAt( fIdx );

	gfMap[fg.nonX].fIdx = fIdx;
	if ( founder.isHaploid( CHR_IS_NOT_X ) )
	    gfMap[fg.nonX++].whichOne = Pedigree::GT_SINGLE;
	else
	    {
	    gfMap[fg.nonX++].whichOne = Pedigree::GT_PATERNAL;
	    gfMap[fg.nonX].fIdx = fIdx;
	    gfMap[fg.nonX++].whichOne = Pedigree::GT_MATERNAL;
	    }

	gfMapX[fg.X].fIdx = fIdx;
	if ( founder.isHaploid( CHR_IS_X ) )
	    gfMapX[fg.X++].whichOne = Pedigree::GT_SINGLE;
	else
	    {
	    gfMapX[fg.X++].whichOne = Pedigree::GT_PATERNAL;
	    gfMapX[fg.X].fIdx = fIdx;
	    gfMapX[fg.X++].whichOne = Pedigree::GT_MATERNAL;
	    }

	}

    gp_assert_eq( fg.nonX     , nFGametes  );
    gp_assert_eq( fg.X	      , nFGametesX );
    gp_assert_eq( fgMap.size(), nFounders  );


    // Count the number of affected members (with separate count of
    // non-founders) for the affected-only test (NOTE *2*):
    for ( Pedigree::Iterator it = Pedigree::getFirstMember() ; it != Pedigree::getEndMember() ; ++it )
	if ( (*it)->getOutcome() == Organism::OUTCOME_AFFECTED )
	    {
	    ++nAffected;
	    if ( ! (*it)->isFounder() )
		++nAffNonFndr;
	    }


    #if 0
	// Cache pointers to the first and ending founders, which should all be at
	// the start of the sortedMembers array:
	firstFounder = sortedMembers;
	endFounder = sortedMembers + nFounders;
    #endif


    // Sanity checking:
    gp_assert( nFounders != 0 );
    gp_assert_ge( nFGametes , nFounders );
    gp_assert_ge( nFGametesX, nFounders );
    gp_assert_le( nFounders, nMembers );
    gp_assert_le( getNNonFndrs(), InheritanceVector::MAX_ORGANISMS );
    gp_assert( sortedMembers[0]->isFounder() );
    gp_assert( getEndFounder()[-1]->isFounder() );

    if ( getEndFounder() != getEndMember() ) // If we have any non-founders at all
	gp_assert( ! (*getEndFounder())->isFounder() );

    ivSpace = new InheritanceSpace( *this );

    }



//-----------------------------------------------------------------------------
// Copy constructor and assignment operator: these are required for
// std::vector<Pedigree> and I believe are safe if never used from anywhere but
// there; they are public, but should never be used except for immediate
// destruction of the copied object.
//-----------------------------------------------------------------------------

Pedigree::Pedigree( const Pedigree & rhs ) :
	memberPool	  ( rhs.memberPool	  ) ,
	id		  ( rhs.id		  ) ,
	nMembers	  ( rhs.nMembers	  ) ,
	nFounders	  ( rhs.nFounders	  ) ,
	sortedMembers	  ( rhs.sortedMembers	  ) ,
	fgMap		  ( rhs.fgMap		  ) ,
	gfMap		  ( rhs.gfMap		  ) ,
	gfMapX		  ( rhs.gfMapX		  ) ,
	ivSpace		  ( rhs.ivSpace		  ) ,
	nMendelErrs	  ( rhs.nMendelErrs	  ) ,
	mendelErrsByLocus ( rhs.mendelErrsByLocus ) ,
	stateProbs	  ( rhs.stateProbs	  ) ,
	nAffected	  ( rhs.nAffected	  ) ,
	nAffNonFndr	  ( rhs.nAffNonFndr	  ) ,
	thetas		  ( rhs.thetas		  ) , // See NOTE *2*
	SumSoftmaxTheta	  ( rhs.SumSoftmaxTheta	  ) , // See NOTE *2*
	rhos		  ( rhs.rhos		  ) , // See NOTE *2*
	sumlogrho	  ( rhs.sumlogrho	  ) , // See NOTE *2*
	step		  ( rhs.step		  ) , // See NOTE *2*
	NumberOfUpdates	  ( rhs.NumberOfUpdates	  ) , // See NOTE *2*
	w		  ( rhs.w		  ) , // See NOTE *2*
	ThetaTuner	  ( rhs.ThetaTuner	  ) , // See NOTE *2*
	NumGametes	  ( rhs.NumGametes	  ) , // See NOTE *2*
	tpCache		  ( rhs.tpCache		  ) , // See NOTE *2*
	hmm		  ( rhs.hmm		  ) , // See NOTE *2*
	llCache		  ( *this, rhs.llCache	  ) ,
	aoCache		  ( rhs.aoCache		  )
    {
    // !!!WARNING!!! -- see NOTE *1*
    const_cast<Pedigree&>(rhs).sortedMembers	 = 0;
    const_cast<Pedigree&>(rhs).ivSpace		 = 0;
    const_cast<Pedigree&>(rhs).mendelErrsByLocus = 0;
    const_cast<Pedigree&>(rhs).stateProbs	 = 0;
    const_cast<Pedigree&>(rhs).tpCache		 = 0; // See NOTE *2*
    const_cast<Pedigree&>(rhs).hmm		 = 0; // See NOTE *2*

    setMyNumber( rhs.myNumber ); // See NOTE *2*
    }


Pedigree & Pedigree::operator=( const Pedigree & rhs )
    {
    gp_assert( &memberPool == &rhs.memberPool );

    id			= rhs.id		;
    nMembers		= rhs.nMembers		;
    nFounders		= rhs.nFounders		;
    sortedMembers	= rhs.sortedMembers	;
    fgMap		= rhs.fgMap		;
    gfMap		= rhs.gfMap		;
    gfMapX		= rhs.gfMapX		;
    ivSpace		= rhs.ivSpace		;
    nMendelErrs		= rhs.nMendelErrs	;
    mendelErrsByLocus	= rhs.mendelErrsByLocus ;
    stateProbs		= rhs.stateProbs	;
    nAffected		= rhs.nAffected		;
    nAffNonFndr		= rhs.nAffNonFndr	;

    // See NOTE *2*:
    thetas		= rhs.thetas		;
    SumSoftmaxTheta	= rhs.SumSoftmaxTheta	;
    rhos		= rhs.rhos		;
    sumlogrho		= rhs.sumlogrho		;
    step		= rhs.step		;
    NumberOfUpdates	= rhs.NumberOfUpdates	;
    w			= rhs.w			;
    ThetaTuner		= rhs.ThetaTuner	;
    NumGametes		= rhs.NumGametes	;
    setMyNumber		( rhs.myNumber		);


    tpCache		= rhs.tpCache		;
    hmm			= rhs.hmm		;

    llCache		= rhs.llCache		;

    aoCache		= rhs.aoCache		;


    // !!!WARNING!!! -- see NOTE *1*
    const_cast<Pedigree&>(rhs).sortedMembers	 = 0;
    const_cast<Pedigree&>(rhs).ivSpace		 = 0;
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
    delete ivSpace;
    delete[] stateProbs;
    delete[] mendelErrsByLocus;
    delete tpCache; // See NOTE *2*
    delete hmm;	    // See NOTE *2*
    }



//---------------------------------------------------------------
// generatePedigrees() [static]
//
/// Create Pedigree's from raw genotype data (pedfile)
//---------------------------------------------------------------

void Pedigree::generatePedigrees( const OrganismArray & organisms, cvector<Pedigree> & rv, size_t max_n )
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

		PED_PB( rv, organisms, firstInd, lastInd, max_n );
		}

	    curFamId = ind.getFamId();
	    }
	}

    gp_assert( ! curFamId.empty() );

    // Ugly data-type but pretty method-call (no "output"/"reference" parameters):
    const pair<OrganismArray::ConstPedIter,OrganismArray::ConstPedIter> &
	interval = organisms.findOrgsInPed( curFamId );

    PED_PB( rv, organisms, interval.first, interval.second, max_n );
    }



//-----------------------------------------------------------------------------
// throwFRange(), throwMRange() [protected]
//-----------------------------------------------------------------------------

void Pedigree::throwFRange( size_t fIdx ) const
    {
    throw out_of_range( estr("Founder-index ") + fIdx +
		" out of range (" + getNFounders() + ')' );
    }

void Pedigree::throwMRange( size_t mIdx ) const
    {
    throw out_of_range( estr("Member-index ") + mIdx +
		" out of range (" + getNMembers() + ')' );
    }



//-----------------------------------------------------------------------------
// gameteTypeDesc() [static]
//-----------------------------------------------------------------------------

const char * Pedigree::gameteTypeDesc( GameteType gt )
    {
    const char * rv;

    if ( gt == GT_PATERNAL )
	rv = "paternal";
    else if ( gt == GT_MATERNAL )
	rv = "maternal";
    else if ( gt == GT_SINGLE )
	rv = "single";
    else
	throw std::runtime_error( "invalid gamete-type" );

    return rv;
    }



//-----------------------------------------------------------------------------
// getGTypeOfFndrGamete()
//-----------------------------------------------------------------------------

Genotype::AlleleType Pedigree::getGTypeOfFndrGamete( FGameteIdx fgIdx, SLocIdxType sLocIdx ) const
    {

    Genotype::AlleleType rv;

    const IsXChromType is_xchrom = getSLoci()[sLocIdx].isXChrom();

    GameteType whichOne;
    const FounderIdx fIdx    = founderOfGameteIdx( fgIdx, whichOne, is_xchrom );
    const Member &   founder = founderAt( fIdx );
    const Genotype & gType   = founder.getGType( sLocIdx );

    if ( whichOne == GT_SINGLE )
	{
	gp_assert( founder.isHaploid( is_xchrom ) );
	gp_assert( gType.isHaploid() );
	rv = gType.getVal1();
	}
    else
	{
	gp_assert( ! founder.isHaploid( is_xchrom ) );
	gp_assert( ! gType.isHaploid() );

	// Is this right?  Arbitrarily call the first allele in the unphased
	// genotype the "paternal" and the second the "maternal"?
	rv = (whichOne == GT_PATERNAL) ? gType.getVal1() : gType.getVal2();
	}

    return rv;

    }



//-----------------------------------------------------------------------------
// haveMendelErrAt()
//-----------------------------------------------------------------------------

bool Pedigree::haveMendelErrAt( SLocIdxType t ) const
    {
    #if AGGRESSIVE_RANGE_CHECK
	if ( t >= getNSLoci() )
	    throw out_of_range( estr("Simple-locus-index ") + t +
		    " out of range (" + getNSLoci() + ')' );
    #endif

    return (nMendelErrs != 0) && mendelErrsByLocus[t];
    }



} // ---- end namespace genepi
