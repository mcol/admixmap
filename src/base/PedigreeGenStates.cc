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
/// \file PedigreeGenStates.cc
/// Implementation of the portion of the Pedigree class to generate the
/// hidden-state-space and emission-probabilities (HiddenStateSpace and
/// HiddenStateSpace::State).  See also Pedigree.cc.
//=============================================================================


#include "Pedigree.h"

#include "AncestryVector.h"
#include "HiddenStateSpace.h"
#include "HiddenStateSpace.h"
#include "InheritanceVector.h"



// These control whether debugging-output support is compiled in at all; to get
// the output, it must also be enabled by calling dbgRecursion() and/or
// dbgEmission():
#define DEBUG_EMISSION_PROBS	1
#define DEBUG_RECURSION		1



#if DEBUG_EMISSION_PROBS || DEBUG_RECURSION
    #include <iostream>
#endif



using namespace std;



namespace genepi { // ----



#if 1 // Compile in the enable/disable methods even if they do nothing:
    static bool dRecursion = false;
    static bool dEmission  = false;
    void Pedigree::dbgRecursion( bool nv ) { dRecursion = nv; }
    void Pedigree::dbgEmission ( bool nv ) { dEmission	= nv; }
#endif


#if DEBUG_EMISSION_PROBS
    /// Output the founder-haplotype-state:
    static std::ostream & output_hs( std::ostream & os, const Haplotype * hapState, const Pedigree & ped, bool full )
	{
	os << "FHS{";
	const Pedigree::MemberIdx limit = full ? ped.getNMembers() : ped.getNFounders();
	for ( Pedigree::MemberIdx fIdx = 0 ; fIdx < limit ; ++fIdx )
	    {
	    const Organism &  org = ped.memberAt( fIdx );
	    const Haplotype & hap = hapState[ fIdx ];
	    if ( fIdx == ped.getNFounders() )
		os << "} IHS{";
	    else if ( fIdx != 0 )
		os << ';';
	    os << org.getOrgId() << '(' << hap.paternal() << ',' << hap.maternal() << ')';
	    }
	return os << '}';
	}
#endif



//==============================================================================
// INTERNAL ARRAY OF STATES
//==============================================================================


//-----------------------------------------------------------------------------
// getStateProbs()
//-----------------------------------------------------------------------------

const HiddenStateSpace & Pedigree::getStateProbs( SLocIdxType sLocIdx ) const
    {
    #if AGGRESSIVE_RANGE_CHECK

	if ( stateProbs == 0 )
	    throw runtime_error( estr("Attempt to retrieve probabilities before calculated:"
		    " call genPossibleStatesInternal() first; pedigree ID: ") + getId() );

	if ( sLocIdx >= getSLoci().size() )
	    throw runtime_error( estr("Simple-locus-index ") + sLocIdx +
			" out of range (" + getSLoci().size() + " max)" );

    #endif

    return stateProbs[ sLocIdx ];
    }



//-----------------------------------------------------------------------------
// releaseStateProbs()
//-----------------------------------------------------------------------------

void Pedigree::releaseStateProbs() const
    {
    delete[] stateProbs;
    stateProbs = 0;
    }



//-----------------------------------------------------------------------------
// mendelErrAt() [private]
//-----------------------------------------------------------------------------

bool & Pedigree::mendelErrAt( SLocIdxType t ) const
    {
    if ( mendelErrsByLocus == 0 )
	{
	mendelErrsByLocus = new bool[ getNSLoci() ];
	for ( SLocIdxType idx = getNSLoci() ; idx-- != 0 ; )
	    mendelErrsByLocus[ idx ] = false;
	}

    return mendelErrsByLocus[ t ];
    }



//-----------------------------------------------------------------------------
//
// accumStateInArray() [static, protected]
//
/// This is just here to receive qualifying states (we pass a pointer-to-method
/// into the recursive iterator/checker), store them in an array, and accumulate
/// the emission probabilities.  If we want to hard-code the store-in-array
/// action, this can go away (as well as the receiver method passed into the
/// iterator/checker).
//
//-----------------------------------------------------------------------------

void Pedigree::accumStateInArray( const Pedigree &	    ped		    ,
				  size_t		    sLocIdx	    ,
				  const AncestryVector &    av		    ,
				  const InheritanceVector & iv		    ,
				  const Haplotype *	    founderHapState ,
				  double		    emProb	    )
    {
    // stateProbs is mutable, so we don't need to const_cast<> here:
    ped.stateProbs[sLocIdx].getEProb( av, iv ) += emProb;

    #if DEBUG_EMISSION_PROBS
	if ( dEmission )
	    {
	    std::cout << "Emission prob term for " << av // << '=' << av.to_ulong()
		<< ' ' << iv // << '=' << iv.to_ulong()
		<< ' ';
	    output_hs( std::cout, founderHapState, ped, true )
		<< " is " << emProb
		<< " total so far: " << ped.getStateProbs(sLocIdx).getEProb(av,iv)
		<< '\n';
	    }
    #else
	if ( founderHapState ) {;} // Suppress unused parameter compiler warning
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
///
/// In fact, this could then be a method of that object; and the receiver could
/// also be a virtual method which subclasses could override.
//
//-----------------------------------------------------------------------------

void Pedigree::recurseSib( SLocIdxType		sLocIdx		,
			   Haplotype *		memberHapState	,
			   const AncestryVector&av		,
			   InheritanceVector &	iv		,
			   MemberIdx		memDepth	,
			   double		emProbTerm	,
			   StateReceiver	receiver	) const
    {

    #if DEBUG_RECURSION
	if ( dRecursion )
	    std::cout << "Recurse sib sloc(" << sLocIdx << ") memDepth(" << memDepth
		<< ") " << av << " ep-term(" << emProbTerm << ")\n";
    #endif


    if ( memDepth == getNMembers() )
	{
	(*receiver)( *this, sLocIdx, av, iv, memberHapState, emProbTerm );
	return; // **** RETURN HERE ****
	}

    Haplotype &	     thisMemberHap   = memberHapState[ memDepth ];
    const Member &   member	     = memberAt( memDepth );
    Haplotype &	     patHapState     = memberHapState[ member.getFather()->getPIdx() ];
    Haplotype &	     matHapState     = memberHapState[ member.getMother()->getPIdx() ];
    const Genotype & myUnphasedGType = member.getGType( sLocIdx );

    // Generate every possible (all 4) segregation-indicator combinations for
    // this node, check them against the observed data, and recurse each that is
    // consistent.

    //***************************************************************************
    //*** What about haploid genotypes here???  Only need one loop, but must
    //*** modify the IV to know about haploid organism-loci to reserve only one
    //*** bit.  This will also affect Pedigree::getNMeiosis().
    //***************************************************************************

    // These two odd-looking for-loops each iterate over the two possible
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
		recurseSib( sLocIdx, memberHapState, av, iv, memDepth + 1, emProbTerm, receiver );
		}

	    if ( matSI == InheritanceVector::SI_MATERNAL )
		break;
	    }

	if ( patSI == InheritanceVector::SI_MATERNAL )
	    break;
	}

    }



//-----------------------------------------------------------------------------
//
// recurseFounder() [private]
//
/// Need more documentation here.
///
/// Effectively, we want 2 nested loops here; one loops over all possible
/// haplotypes for this founder, and one loops over all possible ancestries, and
/// inside the inner loop we recurse again to descend to the next member of the
/// pedigree.  The actual code looks a little messier than this because we may
/// have missing or haploid genotypes, so we need several different loops for
/// the haplotypes.
//
//-----------------------------------------------------------------------------

void Pedigree::recurseFounder(	SLocIdxType		sLocIdx		,
				PopIdx			K		,
				const AlleleProbTable & alProbTab	,
				Haplotype *		memberHapState	,
				AncestryVector &	ancestry	,
				MemberIdx		memDepth	,
				double			probProdSoFar	,
				StateReceiver		receiver	) const
    {

    // If we recurse over founder _gametes_ rather than _haplotypes_, we can
    // save some run time by just getting each of these probability factors from
    // the allele-frequency table only _once_, multiplying only once, then
    // re-using during subsequent recursion for the "second gamete" in the
    // haplotype; currently we will re-look-up the same one several times as we
    // get called with each of the different ancestry/haplotype combinations for
    // the entire haplotype.  It also might reduce the number of places where
    // recursive calls are made, making it easier to do this
    // lookup/multiplication (and the memDepth==getNFounders() check) just
    // _prior_ to the recursive call, rather than just _after_ it (as done
    // here), thus avoiding the ugly memDepth==0 and memDepth-1 code here, and
    // allowing the easy elimination of storing the ancestry values in the
    // ancestry vector, then retrieving them again for lookup here (ditto for
    // alleles).

    if ( memDepth == 0 )
	probProdSoFar = 1.0;
    else
	{
	const FounderIdx fDepth = memDepth - 1; // The depth just before this call was made
	const Haplotype & hap = memberHapState[ fDepth ];
	probProdSoFar *= alProbTab.at( hap.getPaternalAllele(), ancestry.at(fDepth,true ) );
	probProdSoFar *= alProbTab.at( hap.getMaternalAllele(), ancestry.at(fDepth,false) );
	}

    if ( memDepth == getNFounders() )
	{
	InheritanceVector iv( *this );
	recurseSib( sLocIdx, memberHapState, ancestry, iv, memDepth, probProdSoFar, receiver );
	return; // **** RETURN HERE ****
	}

    gp_assert_lt( memDepth, getNFounders() );

    Haplotype &		thisMemberHap	= memberHapState[ memDepth ];
    const Member &	member		= memberAt( memDepth );
    const Genotype &	thisUnphasedGT	= member.getGType( sLocIdx );
    const SimpleLocus & locus		= getSLoci()[ sLocIdx ];
    const size_t	nAlleles	= locus.getNumAlleles();

    gp_assert( member.isFounder() );


    const MemberIdx nextMemDepth = memDepth + 1;


    //-------------------------------------------------------------------------
    // Iterate over every possible ancestry for this node:
    //-------------------------------------------------------------------------

    for ( PopIdx popMat = K ; popMat-- != 0 ; )
	{
	ancestry.setAt( memDepth, false, popMat );

	for ( PopIdx popPat = K ; popPat-- != 0 ; )
	    {
	    ancestry.setAt( memDepth, true, popPat );


	    //------------------------------------------------------------------
	    // Now, for each possible ancestry setting at this node, iterate
	    // over every possible haplotype for this node, and recurse each
	    // one of those:
	    //------------------------------------------------------------------

	    // ---- How to deal with haploid? ----

	    if ( thisUnphasedGT.hasOneVal() )
		{
		// Only one observed allele in genotype: haploid or one missing:
		const Genotype::AlleleType obsVal = thisUnphasedGT.getVal1();
		for ( Genotype::AlleleType aIdx = 1; aIdx <= nAlleles; ++aIdx )
		    {
		    thisMemberHap.setVals( aIdx, obsVal );
		    recurseFounder( sLocIdx, K, alProbTab, memberHapState, ancestry, nextMemDepth, probProdSoFar, receiver );
		    if ( aIdx != obsVal )
			{
			thisMemberHap.setVals( obsVal, aIdx );
			recurseFounder( sLocIdx, K, alProbTab, memberHapState, ancestry, nextMemDepth, probProdSoFar, receiver );
			}
		    }
		}
	    else if ( thisUnphasedGT.hasTwoVals() )
		{
		// Two observed alleles in genotype:
		const Genotype::AlleleType a1 = thisUnphasedGT.getVal1();
		const Genotype::AlleleType a2 = thisUnphasedGT.getVal2();
		thisMemberHap.setVals( a1, a2 );
		recurseFounder( sLocIdx, K, alProbTab, memberHapState, ancestry, nextMemDepth, probProdSoFar, receiver );
		if ( a1 != a2 )
		    {
		    thisMemberHap.setVals( a2, a1 );
		    recurseFounder( sLocIdx, K, alProbTab, memberHapState, ancestry, nextMemDepth, probProdSoFar, receiver );
		    }
		}
	    else
		{
		// Untyped (both alleles missing):
		for ( Genotype::AlleleType pIdx = 1; pIdx <= nAlleles; ++pIdx )
		    for ( Genotype::AlleleType mIdx = 1; mIdx <= nAlleles; ++mIdx )
			{
			thisMemberHap.setVals( pIdx, mIdx );
			recurseFounder( sLocIdx, K, alProbTab, memberHapState, ancestry, nextMemDepth, probProdSoFar, receiver );
			}
		}

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

void Pedigree::genPossibleStates( StateReceiver receiver, PopIdx K, const AlleleProbTable & alProbTab, SLocIdxType sLocIdx ) const
    {
    Haplotype	   memberHapState[ getNMembers() ];
    AncestryVector ancestry( *this, K );

    #if DEBUG_EMISSION_PROBS || DEBUG_RECURSION
	if ( dEmission || dRecursion )
	    {
	    std::cout.setf(ios::fixed);
	    std::cout.precision(8);
	    }
    #endif

    recurseFounder( sLocIdx, K, alProbTab, memberHapState, ancestry, 0, 0.0, receiver );
    }



//-----------------------------------------------------------------------------
// genPossibleStates()
//
/// Same as genPossibleStates(StateReceiver,size_t)const but does all simple-loci.
//-----------------------------------------------------------------------------

void Pedigree::genPossibleStates( StateReceiver receiver, PopIdx K, const AlleleProbVect & alProbVect ) const
    {
    Haplotype	   memberHapState[ getNMembers() ];
    AncestryVector ancestry( *this, K );

    #if DEBUG_EMISSION_PROBS || DEBUG_RECURSION
	if ( dEmission || dRecursion )
	    {
	    std::cout.setf(ios::fixed);
	    std::cout.precision(8);
	    }
    #endif

    for ( size_t sLocIdx = 0 ; sLocIdx < getSLoci().size() ; ++sLocIdx )
	recurseFounder( sLocIdx, K, alProbVect.at(sLocIdx), memberHapState, ancestry, 0, 0.0, receiver );
    }



//-----------------------------------------------------------------------------
// genPossibleStatesInternal()
//-----------------------------------------------------------------------------

void Pedigree::genPossibleStatesInternal( PopIdx K, const AlleleProbVect & alProbVect ) const
    {
    if ( stateProbs == 0 )
	{
	const size_t nLoc = getSLoci().size();

	// We want this, but ISO C++ forbids it:
	//stateProbs = new HiddenStateSpace[ nLoc ]( *this, K );

	// So, instead we have this:
	stateProbs = new HiddenStateSpace[ nLoc ];
	for ( size_t sLocIdx = 0 ; sLocIdx < nLoc ; ++sLocIdx )
	    stateProbs[ sLocIdx ] . init( *this, K );
	}

    for ( size_t sLocIdx = 0 ; sLocIdx < getSLoci().size() ; ++sLocIdx )
	stateProbs[ sLocIdx ] . resetEmProbsToZero();

    genPossibleStates( &Pedigree::accumStateInArray, K, alProbVect );


    //------------------------------------------------------------------
    // Check for Mendelian inconsistencies:
    //------------------------------------------------------------------
    for ( size_t sLocIdx = 0 ; sLocIdx < getSLoci().size() ; ++sLocIdx )
	if ( (mendelErrAt(sLocIdx) = HiddenStateSpace::Iterator(stateProbs[sLocIdx]).isFinished()) )
	    {
	    cout << estr("Apparent Mendelian inconsistency in pedigree ") +
		getId() + " at locus #" + sLocIdx + " (" +
		getSLoci()[sLocIdx].getDesc() + ")\n";
	    ++nMendelErrs;
	    }

    }




} // ---- end namespace genepi
