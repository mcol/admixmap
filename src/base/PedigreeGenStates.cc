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
#include "InheritanceVector.h"



// These control whether debugging-output support is compiled in at all; to get
// the output, it must also be enabled by calling dbgRecursion() and/or
// dbgEmission():
#define DEBUG_EMISSION_PROBS	1
#define DEBUG_RECURSION		1

#define DEBUG_INHERITANCE	0


#if DEBUG_EMISSION_PROBS || DEBUG_RECURSION || DEBUG_INHERITANCE
    #include <iostream>
#endif



using namespace std;
typedef genepi::Genotype::AlleleType AlleleType;



namespace genepi { // ----



#if 1 // Compile in the enable/disable methods even if they do nothing:
    static bool dRecursion = false;
    static bool dEmission  = false;
    void Pedigree::dbgRecursion( bool nv ) { dRecursion = nv; }
    void Pedigree::dbgEmission ( bool nv ) { dEmission	= nv; }
#endif


#if DEBUG_EMISSION_PROBS
    /// Output the founder-haplotype-state:
    static std::ostream & output_hs( std::ostream & os, const Haplotype * nFndrHapState, const Pedigree & ped )
	{
	os << "SibHapSt{";
	for ( Pedigree::MemberIdx mIdx = ped.getNFounders() ; mIdx < ped.getNMembers() ; ++mIdx )
	    {
	    const Organism &  org = ped.memberAt( mIdx );
	    const Haplotype & hap = nFndrHapState[ mIdx ];
	    if ( mIdx != ped.getNFounders() )
		os << ';';
	    if ( org.isHaploid() )
		os << org.getOrgId() << '(' << hap.getVal1() << ')';
	    else
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



#if DEBUG_RECURSION

    static std::ostream & operator<<( std::ostream & os, const genepi::Genotype & gt )
	{
	os << '[';
	if ( ! gt.isMissing1() )
	    os << gt.getVal1();
	os << ',';
	if ( ! gt.isMissing2() )
	    os << gt.getVal2();
	os << ']';
	return os;
	}

    static std::ostream & operator<<( std::ostream & os, const genepi::Haplotype & ht )
	{
	return os << static_cast<const genepi::Genotype &>(ht);
	}

#endif



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
				  const Haplotype *	    nonFndrHapState ,
				  double		    emProb	    )
    {

    if ( emProb == 0.0 )
	cerr << "Strangeness: pedigree " << ped.getId() << " at locus " << sLocIdx << " (" <<
	    ped.getSLoci()[sLocIdx].getName() << ") received a zero emission-probability.\n";

    // stateProbs is mutable, so we don't need to const_cast<> here:
    #if 0 // addToEProbAt() tracks the number of non-zero e-probs
	ped.stateProbs[sLocIdx].getEProb( av, iv ) += emProb;
    #else
	ped.stateProbs[sLocIdx].addToEProbAt( av, iv, emProb );
    #endif

    #if DEBUG_EMISSION_PROBS
	if ( dEmission )
	    {
	    std::cout << "Emission prob term for " << av << ' ' << iv << ' ';
	    output_hs( std::cout, nonFndrHapState, ped )
		<< " is " << emProb
		<< " total so far: " << ped.getStateProbs(sLocIdx).getEProb(av,iv)
		<< "; number of non-zero so far: " << ped.getStateProbs(sLocIdx).getNNon0()
		<< '\n';
	    }
    #else
	if ( nonFndrHapState ) {;} // Suppress unused parameter compiler warning
    #endif

    }



//==============================================================================
// GENERATION OF STATES
//==============================================================================


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

void Pedigree::recurseSib( SLocIdxType		  sLocIdx	 ,
			   Genotype::AlleleType * fndrGameteState,
			   Haplotype *		  nonFndrHapState,
			   const AncestryVector & av		 ,
			   InheritanceVector &	  iv		 ,
			   MemberIdx		  memDepth	 ,
			   double		  emProbTerm	 ,
			   StateReceiver	  receiver	 ) const
    {

    #if DEBUG_RECURSION
	if ( dRecursion )
	    std::cout << "Recurse-sib memDepth(" << memDepth
		<< ") " << av << " iv-so-far(" << iv << ") ep-term(" << emProbTerm << ")\n";
    #endif


    if ( memDepth == getNMembers() )
	{
	(*receiver)( *this, sLocIdx, av, iv, nonFndrHapState, emProbTerm );
	return; // **** RETURN HERE ****
	}

    Haplotype &	     thisMemberHap   = nonFndrHapState[ memDepth ];
    const Member &   member	     = memberAt( memDepth );
    const Genotype & myUnphasedGType = member.getGType( sLocIdx );

    // Generate every possible (all 4) segregation-indicator combinations for
    // this node, check them against the observed data, and recurse each that is
    // consistent.

    //***************************************************************************
    //*** What about haploid genotypes here???  Only need one loop, but must
    //*** modify the IV to know about haploid organism-loci to reserve only one
    //*** bit.  This will also affect Pedigree::getNMeiosis().  Alternatively,
    //*** we could use the full complement of 4 possible segregation indicators,
    //*** even though they are not meaningful and we are overstating the number
    //*** of possible states, perhaps the resulting probabilities will not be
    //*** incorrect, i.e. the sum of the probabilities of the 2 states with the
    //*** same value for the paternal SI and (incorrectly) varying maternal SI
    //*** will sum to the correct probability of what would have been the single
    //*** state containing that value of the paternal SI and no maternal SI.
    //***************************************************************************

    // These odd-looking for-loops each iterate over the two possible values of
    // each of the two segregation-indicators.

    for ( InheritanceVector::SegInd patSI = InheritanceVector::SI_PATERNAL ; ;
	    patSI = InheritanceVector::SI_MATERNAL )
	{

	const Member & father = *member.getFather();

	if ( father.isHaploid() )
	    {
	    if ( patSI == InheritanceVector::SI_MATERNAL )
		break;
	    else
		patSI = InheritanceVector::SI_NONE;
	    }

	for ( InheritanceVector::SegInd matSI = InheritanceVector::SI_PATERNAL ; ;
		matSI = InheritanceVector::SI_MATERNAL )
	    {

	    const Member & mother = *member.getMother();

	    if ( mother.isHaploid() )
		{
		if ( matSI == InheritanceVector::SI_MATERNAL )
		    break;
		else
		    matSI = InheritanceVector::SI_NONE;
		}


	    #if DEBUG_RECURSION
		if ( dRecursion )
		    {
		    std::cout
			<< "  propose: paternal-seg-ind=" << patSI
			<< "; maternal-srg-ind=" << matSI << '\n';
		    }
	    #endif

	    AlleleType paternalGamete;
	    if ( patSI == InheritanceVector::SI_NONE )
		{
		gp_assert( father.isHaploid() );
		if ( father.isFounder() )
		    paternalGamete = fndrGameteState[ founderGameteOfFounder(father.getPIdx(),GT_SINGLE) ];
		else
		    {
		    // For the moment, no haploid non-founders; this needs to be fixed for X-chromosomes for pedigrees.
		    throw std::runtime_error( "No haploid non-founders" );
		    }
		}
	    else
		{
		gp_assert( ! father.isHaploid() );
		if ( father.isFounder() )
		    {
		    const FGameteIdx patGameteIdx = founderGameteOfFounder( father.getPIdx(),
						    (patSI == InheritanceVector::SI_PATERNAL) ? GT_PATERNAL : GT_MATERNAL );
		    paternalGamete = fndrGameteState[ patGameteIdx ];
		    #if DEBUG_INHERITANCE
			std::cout << "    father.pedIdx=" << father.getPIdx()
			    << ", patGameteIdx=" << patGameteIdx
			    << " -> paternalGamete=" << paternalGamete << '\n';
		    #endif
		    }
		else
		    {
		    Haplotype & paternalHap = nonFndrHapState[ father.getPIdx() ];
		    paternalGamete =
			(patSI == InheritanceVector::SI_PATERNAL) ?
			    paternalHap.getPaternalAllele()	  :
			    paternalHap.getMaternalAllele()	  ;
		    }
		}


	    AlleleType maternalGamete;
	    if ( matSI == InheritanceVector::SI_NONE )
		{
		gp_assert( mother.isHaploid() );
		if ( mother.isFounder() )
		    maternalGamete = fndrGameteState[ founderGameteOfFounder(mother.getPIdx(),GT_SINGLE) ];
		else
		    {
		    // For the moment, no haploid non-founders; this needs to be fixed for X-chromosomes for pedigrees.
		    throw std::runtime_error( "No haploid non-founders" );
		    }
		}
	    else
		{
		gp_assert( ! mother.isHaploid() );
		if ( mother.isFounder() )
		    {
		    const FGameteIdx matGameteIdx = founderGameteOfFounder( mother.getPIdx(),
						    (matSI == InheritanceVector::SI_PATERNAL) ? GT_PATERNAL : GT_MATERNAL );
		    maternalGamete = fndrGameteState[ matGameteIdx ];
		    #if DEBUG_INHERITANCE
			std::cout << "    mother.pedIdx=" << mother.getPIdx()
				<< ", matGameteIdx=" << matGameteIdx
				<< " -> maternalGamete=" << maternalGamete << '\n';
		    #endif
		    }
		else
		    {
		    Haplotype & maternalHap = nonFndrHapState[ mother.getPIdx() ];
		    maternalGamete =
			(matSI == InheritanceVector::SI_PATERNAL) ?
			    maternalHap.getPaternalAllele()	  :
			    maternalHap.getMaternalAllele()	  ;
		    }
		}


	    const InheritanceVector::Bits ivSeg( patSI, matSI );
	    thisMemberHap.setVals( paternalGamete, maternalGamete );


	    // One possible added optimization: if two different IVs produce the
	    // same haplotype here, we don't need to recurse each separately,
	    // but we may need a way to express the set of possible resulting
	    // IVs as an equivalence class of states, so for the moment we'll
	    // generate all of them separately.

	    #if DEBUG_RECURSION
		if ( dRecursion )
		    std::cout
			<< "  check consistency of proposed phased nonf: "
			<< thisMemberHap
			<< " versus actual unphased: " << myUnphasedGType
			<< ": " << (myUnphasedGType.consistent(thisMemberHap)?"yes":"no")
			<< '\n';
	    #endif

	    if ( myUnphasedGType.consistent( thisMemberHap ) )
		{
		// Store the bits in the IV so that if/when we reach the bottom
		// of the recursion, we'll have accumulated the whole qualifying
		// IV for emission:
		iv.setMember( memDepth, ivSeg );
		recurseSib( sLocIdx, fndrGameteState, nonFndrHapState, av, iv, memDepth + 1, emProbTerm, receiver );
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

void Pedigree::recurseFndrGamete( SLocIdxType		  sLocIdx	  ,
				  PopIdx		  K		  ,
				  const AlleleProbTable & alProbTab	  ,
				  AlleleType *		  fndrGameteState ,
				  AncestryVector &	  ancestry	  ,
				  FGameteIdx		  fgDepth	  ,
				  double		  probProdSoFar   ,
				  StateReceiver		  receiver	  ) const
    {

    #if DEBUG_RECURSION
	if ( dRecursion )
	    {
	    std::cout << "Recurse-fndr fgDepth(" << fgDepth << ':';
	    for ( FGameteIdx fgi = 0 ; fgi < fgDepth ; ++fgi )
		{
		if ( fgi != 0 )
		    std::cout << ',';
		std::cout << fndrGameteState[fgi];
		}
	    std::cout << ") ";
	    output(std::cout,ancestry,fgDepth) << " prob-prod-so-far(" << probProdSoFar << ")\n";
	    }
    #endif


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

    if ( fgDepth == 0 )
	probProdSoFar = 1.0;
    else
	{
	const FGameteIdx pDepth = fgDepth - 1; // The depth just before this call was made
	const AlleleType & al = fndrGameteState[ pDepth ];
	probProdSoFar *= alProbTab.at( al, ancestry.at(pDepth) );
	}


    // If recursed to the end of the founders, start recursing the non-founders:
    if ( fgDepth == getNFounderGametes() )
	{

	InheritanceVector iv( *this );

	// We currently index the non-founder-haplotype-state stack on
	// member-index rather than non-founder-index, so we must allocate more
	// than we need.
	#if 0
	    Haplotype nonFndrHapState[ getNNonFndrs() ];
	#else
	    Haplotype nonFndrHapState[ getNMembers() ];
	#endif

	// Start memDepth at nFounders, which is the member-index of the first non-founder.
	recurseSib( sLocIdx, fndrGameteState, nonFndrHapState, ancestry, iv, getNFounders(), probProdSoFar, receiver );

	return; // **** RETURN HERE ****

	}

    gp_assert_lt( fgDepth, getNFounderGametes() );

    const SimpleLocus & locus		  = getSLoci()[ sLocIdx ];
    const size_t	nAlleles	  = locus.getNumAlleles();
    const bool		isX		  = locus.isXChrom();
    AlleleType &	thisAllele	  = fndrGameteState[ fgDepth ];
    GameteType		whichOne;
    const FounderIdx	fIdx		  = founderOfGameteIdx( fgDepth, whichOne, isX );
    const Member &	founder		  = founderAt( fIdx );
    const bool		haplotypeComplete = (whichOne != GT_PATERNAL); // NOTE *X91* in Pedigree.h
    const Genotype &	thisUnphasedGT	  = founder.getGType( sLocIdx );

    gp_assert( founder.isFounder() );
    gp_assert( (whichOne == GT_SINGLE) || (! founder.isHaploid()) );


    //-------------------------------------------------------------------------
    // Iterate over every possible ancestry for this node:
    //-------------------------------------------------------------------------

    for ( PopIdx pop = K ; pop-- != 0 ; )
	{

	ancestry.setAt( fgDepth, pop );


	//------------------------------------------------------------------
	// Now, for each possible ancestry setting at this node, iterate over
	// every possible allele for this founder-gamete, and recurse each one
	// of those; yet eliminating any phased-haplotype combinations that are
	// not consistent with observed genotype data for this founder.
	//------------------------------------------------------------------

	for ( Genotype::AlleleType al = 1; al <= nAlleles; ++al )
	    if ( thisUnphasedGT.consistent( al ) )
		{
		bool consistent = (! haplotypeComplete) || founder.isHaploid();
		if ( ! consistent )
		    {
		    gp_assert( fgDepth > 0 );
		    const FGameteIdx pDepth = fgDepth - 1; // The depth just before this call was made
		    consistent = thisUnphasedGT.consistent( Haplotype(fndrGameteState[pDepth],al) );
		    #if DEBUG_RECURSION
			if ( dRecursion )
			    std::cout
				<< "  check consistency of proposed phased fndr: "
				<< Haplotype( fndrGameteState[pDepth], al )
				<< " versus actual unphased: " << thisUnphasedGT
				<< ": " << (consistent?"yes":"no") << '\n';
		    #endif
		    }
		if ( consistent )
		    {
		    thisAllele = al;
		    recurseFndrGamete( sLocIdx, K, alProbTab, fndrGameteState, ancestry,
				    fgDepth + 1, probProdSoFar, receiver );
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
    AlleleType	   fndrGameteState[ getNFounderGametes() ];
    AncestryVector ancestry( *this, K );

    #if DEBUG_EMISSION_PROBS || DEBUG_RECURSION
	if ( dEmission || dRecursion )
	    {
	    std::cout.setf(ios::fixed);
	    std::cout.precision(8);
	    std::cout << "\n\n====== LOCUS #" << sLocIdx << ' '
		      << getSLoci()[sLocIdx].getName() << " =======\n\n";
	    }
    #endif

    recurseFndrGamete( sLocIdx, K, alProbTab, fndrGameteState, ancestry, 0, 0.0, receiver );
    }



//-----------------------------------------------------------------------------
// genPossibleStates()
//
/// Same as genPossibleStates(StateReceiver,size_t)const but does all simple-loci.
//-----------------------------------------------------------------------------

void Pedigree::genPossibleStates( StateReceiver receiver, PopIdx K, const AlleleProbVect & alProbVect ) const
    {
    AlleleType	   fndrGameteState[ getNFounderGametes() ];
    AncestryVector ancestry( *this, K );

    #if DEBUG_EMISSION_PROBS || DEBUG_RECURSION
	if ( dEmission || dRecursion )
	    {
	    std::cout.setf(ios::fixed);
	    std::cout.precision(8);
	    }
    #endif

    for ( size_t sLocIdx = 0 ; sLocIdx < getSLoci().size() ; ++sLocIdx )
	{
	#if DEBUG_EMISSION_PROBS || DEBUG_RECURSION
	    if ( dEmission || dRecursion )
		std::cout << "\n\n====== LOCUS #" << sLocIdx << ' '
			  << getSLoci()[sLocIdx].getName() << " =======\n\n";
	#endif
	recurseFndrGamete( sLocIdx, K, alProbVect.at(sLocIdx), fndrGameteState, ancestry, 0, 0.0, receiver );
	}
    }



//-----------------------------------------------------------------------------
// genPossibleStatesInternal()
//-----------------------------------------------------------------------------

void Pedigree::genPossibleStatesInternal( PopIdx K, const AlleleProbVect & alProbVect ) const
    {
    if ( stateProbs == 0 )
	{
	const size_t nLoc = getSLoci().size();

	#if 0
	    // We want this, but ISO C++ forbids it:
	    stateProbs = new HiddenStateSpace[ nLoc ]( *this, K );
	#else
	    // So, instead we have this:
	    stateProbs = new HiddenStateSpace[ nLoc ];
	    for ( size_t sLocIdx = 0 ; sLocIdx < nLoc ; ++sLocIdx )
		stateProbs[ sLocIdx ] . init( *this, K );
	#endif
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
		getId() + " at locus #" + sLocIdx + " [" +
		getSLoci()[sLocIdx].getDesc() + "]\n";
	    ++nMendelErrs;
	    }

    }




} // ---- end namespace genepi
