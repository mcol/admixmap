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
/// \file PedigreeGenStates.cc
/// Implementation of the generation of founder-haplotype-set +
/// inheritance-vector states portion of the Pedigree class.  See also
/// Pedigree.cc.
//=============================================================================


#include "Pedigree.h"



using namespace std;



namespace genepi { // ----



//==============================================================================
// INTERNAL ARRAY OF STATES
//==============================================================================


//-----------------------------------------------------------------------------
// getConsistentStates()
//-----------------------------------------------------------------------------

const vector<State> & Pedigree::getConsistentStates( size_t sLocIdx ) const
    {
    #if AGGRESSIVE_RANGE_CHECK
	if ( sLocIdx >= getSLoci().size() )
	    throw std::runtime_error( estr("Simple-locus-index ") + sLocIdx +
			" out of range (" + getSLoci().size() );
    #endif

    return consistentStates[ sLocIdx ];
    }



//-----------------------------------------------------------------------------
//
// accumStateInArray() [static]
//
/// This is just here to receive qualifying states and store them in an array
/// (we pass a pointer-to-method into the recursive iterator/checker).  If we
/// want to hard-code the store-in-array action, this can go away (as well as
/// the receiver method passed into the iterator/checker).
//
//-----------------------------------------------------------------------------

void Pedigree::accumStateInArray( const Pedigree &	    ped		    ,
				  size_t		    sLocIdx	    ,
				  const InheritanceVector & iv		    ,
				  const Haplotype *	    founderHapState )
    {
    #if 1
	const_cast<Pedigree&>(ped).consistentStates[sLocIdx].push_back( State(iv,founderHapState) );
    #else // Why doesn't this autobox work?
	const_cast<Pedigree&>(ped).consistentStates[sLocIdx].push_back( iv, founderHapState );
    #endif
    }



//==============================================================================
// GENERATION OF STATES
//==============================================================================


//-----------------------------------------------------------------------------
/// Compute inherited haplotype from the relevant bits of an inheritance vector
/// and the two parents' phased haplotypes:
//-----------------------------------------------------------------------------

static inline Haplotype inheritHap(
		   const InheritanceVector::Bits & ivEl		,
		   const Haplotype &		   paternalHap	,
		   const Haplotype &		   maternalHap	)
    {
    return Haplotype(
		(ivEl.paternal() == InheritanceVector::SI_PATERNAL) ?
		    paternalHap.getPaternalAllele()		    :
		    paternalHap.getMaternalAllele()			,
		(ivEl.maternal() == InheritanceVector::SI_PATERNAL) ?
		    maternalHap.getPaternalAllele()		    :
		    maternalHap.getMaternalAllele()			);
    }



//-----------------------------------------------------------------------------
//
// recurseSib() [private]
//
/// Need more documentation here.
///
/// Perhaps all these parameters should be wrapped together into one
/// "founder-state-enumeration recursion state" object so we can just
/// pass a pointer to that into here.
//
//-----------------------------------------------------------------------------

void Pedigree::recurseSib( size_t		sLocIdx		,
			   Haplotype *		memberHapState	,
			   InheritanceVector &	iv		,
			   MemberIdx		memDepth	,
			   StateReceiver	receiver	) const
    {

    if ( memDepth == getNMembers() )
	{
	(*receiver)( *this, sLocIdx, iv, memberHapState );
	return; // **** RETURN HERE ****
	}

    Haplotype &	     thisMemberHap   = memberHapState[ memDepth ];
    const Member &   member	     = memberAt( memDepth );
    Haplotype &	     patHapState     = memberHapState[ member.getFather()->getPIdx() ];
    Haplotype &	     matHapState     = memberHapState[ member.getMother()->getPIdx() ];
    const Genotype & myUnphasedGType = member.getGType( sLocIdx );

    // Generate every possible (all 4) segregation-indicator combination for
    // this node, check them against the observed data, and recurse each that is
    // consistent.

    // These two odd-looking for loops each iterate over the two possible
    // segregation-indicators:

    for ( InheritanceVector::SegInd patSI = InheritanceVector::SI_PATERNAL ; ;
	    patSI = InheritanceVector::SI_MATERNAL )
	{
	for ( InheritanceVector::SegInd matSI = InheritanceVector::SI_PATERNAL ; ;
		matSI = InheritanceVector::SI_MATERNAL )
	    {
	    const InheritanceVector::Bits ivSeg( patSI, matSI );
	    thisMemberHap = inheritHap( ivSeg, patHapState, matHapState );

	    // One possible added optimization: if two different IVs produce the
	    // same haplotype here, we don't need to recurse each separately,
	    // but we may need a way to express the set of possible IVs as an
	    // equivalence class, so for the moment we'll use all separately.

	    if ( myUnphasedGType.consistent( thisMemberHap ) )
		{
		// Store the bits in the IV so that if/when we reach the bottom
		// of the recursion, we'll have accumulated the whole qualifying
		// IV for emission:
		iv.setMember( memDepth, ivSeg );
		recurseSib( sLocIdx, memberHapState, iv, memDepth + 1, receiver );
		}

	    if ( matSI == InheritanceVector::SI_MATERNAL )
		break;
	    }
	if ( patSI == InheritanceVector::SI_MATERNAL )
	    break;
	}

    }



//-----------------------------------------------------------------------------
// recurseFounder()
/// Need documentation here.
//-----------------------------------------------------------------------------

void Pedigree::recurseFounder(	size_t		sLocIdx		,
				Haplotype *	memberHapState  ,
				MemberIdx	memDepth	,
				StateReceiver	receiver	) const
    {

    if ( memDepth == getNFounders() )
	{
	InheritanceVector iv( *this );
	recurseSib( sLocIdx, memberHapState, iv, memDepth, receiver );
	return; // **** RETURN HERE ****
	}

    gp_assert_lt( memDepth, getNFounders() );

    Haplotype &		thisMemberHap	= memberHapState[ memDepth ];
    const Member &	member		= memberAt( memDepth++ );
    const Genotype &	thisUnphasedGT	= member.getGType( sLocIdx );
    const SimpleLocus & locus		= getSLoci()[ sLocIdx ];
    const size_t	nAlleles	= locus.getNumAlleles();

    gp_assert( member.isFounder() );


    // Iterate over every possible haplotype for this node, and recurse each
    // one:

    // ---- How to deal with haploid? ----

    if ( thisUnphasedGT.hasOneVal() )
	{
	// Only one observed allele in genotype: haploid or one missing:
	const Genotype::AlleleType obsVal = thisUnphasedGT.getVal1();
	for ( Genotype::AlleleType aIdx = 1; aIdx <= nAlleles; ++aIdx )
	    {
	    thisMemberHap.setVals( aIdx, obsVal );
	    recurseFounder( sLocIdx, memberHapState, memDepth, receiver );
	    if ( aIdx != obsVal )
		{
		thisMemberHap.setVals( obsVal, aIdx );
		recurseFounder( sLocIdx, memberHapState, memDepth, receiver  );
		}
	    }
	}
    else if ( thisUnphasedGT.hasTwoVals() )
	{
	// Two observed alleles in genotype:
	const Genotype::AlleleType a1 = thisUnphasedGT.getVal1();
	const Genotype::AlleleType a2 = thisUnphasedGT.getVal2();
	thisMemberHap.setVals( a1, a2 );
	recurseFounder( sLocIdx, memberHapState, memDepth, receiver  );
	if ( a1 != a2 )
	    {
	    thisMemberHap.setVals( a2, a1 );
	    recurseFounder( sLocIdx, memberHapState, memDepth, receiver	 );
	    }
	}
    else
	{
	// Untyped (both alleles missing):
	for ( Genotype::AlleleType pIdx = 1; pIdx <= nAlleles; ++pIdx )
	    for ( Genotype::AlleleType mIdx = 1; mIdx <= nAlleles; ++mIdx )
		{
		thisMemberHap.setVals( pIdx, mIdx );
		recurseFounder( sLocIdx, memberHapState, memDepth, receiver  );
		}
	}

    }



//-----------------------------------------------------------------------------
// genPossibleStates()
/// This is typically the public entry-point for generation the space of
/// possible (founder-haplotype-set,inheritance-vector) pairs for each locus
/// ("possible" meaning consistent with the observed unphased genotype data).
/// For each consistent pair, it will call the pointer-to-method receiver .  See
/// also genPossibleStatesInt() which does not require this parameter.
//-----------------------------------------------------------------------------

void Pedigree::genPossibleStates( StateReceiver receiver, size_t sLocIdx ) const
    {
    Haplotype memberHapState[ getNMembers() ];
    recurseFounder( sLocIdx, memberHapState, 0, receiver );
    }



//-----------------------------------------------------------------------------
// genPossibleStates()
/// Same as genPossibleStates(StateReceiver,size_t)const but does all simple-loci.
//-----------------------------------------------------------------------------

void Pedigree::genPossibleStates( StateReceiver receiver ) const
    {
    Haplotype memberHapState[ getNMembers() ];

    for ( size_t sLocIdx = 0 ; sLocIdx < getSLoci().size() ; ++sLocIdx )
	recurseFounder( sLocIdx, memberHapState, 0, receiver );
    }



} // ---- end namespace genepi
