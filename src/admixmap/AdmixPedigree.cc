//=============================================================================
//
// Copyright (c) 2002-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
// Portions Copyright (C) 2009 David D. Favro
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 as published by the Free
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
/// \file PedigreeAdmix.cc
/// Implementation of the admixmap-specific methods of the Pedigree class.
///
/// Additional methods to run the admixmap model.  These should be moved
/// into a derived class, not here in Pedigree.
//
//=============================================================================


#include "Pedigree.h"


#include <cmath>    // exp(), log()
#include <limits>   // numeric_limits<>

#include <bclib/misc.h>
#include <bclib/rand.h>

#include "config.h"			// AGGRESSIVE_RANGE_CHECK
#include "TransProbCache.h"
#include "HiddenMarkovModel.new.h"
#include "AdmixIndividualCollection.h"	// PARALLELIZE_PEDIGREE_LOOP, SAMPLE_THETA_CALL_CRITICAL


#define ALL_FOUNDERS_ARE_ADMIXED    1
#if ALL_FOUNDERS_ARE_ADMIXED
    /// For pedigrees, we assume that all founders are admixed.
    #define IS_ADMIXED()  true
#else
    /// The way it was for individuals.
    #define IS_ADMIXED()  options.isAdmixed( this_gamete )
#endif


#define SEPARATE_HMM_FOR_EACH_CHROM 0

// This is very broken, don't use it unless you fix it:
#define USE_AO_CACHE 0

// Print the details of the M-calculation
#define DEBUG_CALC_M 0


#if 0
    #define DEBUG_TH_PROP(X) X
#else
    #define DEBUG_TH_PROP(X)
#endif

#if AGGRESSIVE_RANGE_CHECK
    #define AGGRESSIVE_ONLY(X)	X
#else
    #define AGGRESSIVE_ONLY(X)
#endif



//------------------------------------------------------------------------
// NOTE X1:
// These were defined as forward references in the header file:
#include "AdmixOptions.h"
#include "CopyNumberAssocTest.h"
#include "bclib/DataMatrix.h"
#include "Chromosome.h"
#include "AdmixOptions.h"
#include "Individual.h" // For Individual::Loci
#include "AffectedsOnlyTest.h"
//------------------------------------------------------------------------



using bclib::pvector;
typedef pvector<double> pvectord;

#if ! PED_HAS_OWN_PRNG
    using bclib::Rand;
#endif


static const double SOFTMAX_0_FLAG = std::numeric_limits<double>::infinity();



// Stages of implementing PedBase methods, ported over from Individual/AdmixedIndividual.
#define NEEDED_ONLY_FOR_CONJUGATE_UPDATE		    0
#define NOT_NEEDED_FOR_FIRST_VERSION			    0
#define USE_SINGLE_POPULATION_SPECIAL_CASE		    0
#define NON_GLOBAL_RHO_WORKS				    0
#define NOT_NEEDED_FOR_FIXEDALLELEFREQ_EQ_1		    0
#define NOT_NEEDED_UNLESS_CONJUGATE_UPDATE		    0
#define SUPPORT_ASSOCIATION_TESTS			    0
#define TRACK_PEDIGREE_MISSING				    0
#define PEDIGREES_HAVE_PLOIDINESS			    0
#define SUM_PROBS_IS_IMPLEMENTED			    0


#if ! NOT_NEEDED_FOR_FIRST_VERSION
    #define IAmUnderTest    0
#endif



namespace genepi { // ----



// Documented in header.
const AdmixOptions * Pedigree::optionsRef = 0;



//-----------------------------------------------------------------------------
// getGenome() [static]
/// This should probably be replaced with a per-instance pointer to a context object.
//-----------------------------------------------------------------------------

inline static Genome & getGenome()
    {
    return Individual::getGenome();
    }



//-----------------------------------------------------------------------------
// getK() [perhaps could be static but currently isn't]
// This hack implementation should be replaced by a per-pedigree reference to a
// context object.
//-----------------------------------------------------------------------------

PopIdx Pedigree::K = 0;

void Pedigree::setK( PopIdx _K )
    {
    if ( K != 0 )
	throw std::runtime_error( genepi::estr("Set K to ") + _K
				    + " after setting to " + K );

    K = _K;
    }



//-----------------------------------------------------------------------------
// getNumChromosomes() [perhaps could be static but currently isn't]
// Context access.
//-----------------------------------------------------------------------------

ChromIdxType Pedigree::getNumChromosomes() const
    {
    return getSLoci().getNChromosomes();
    }



//-----------------------------------------------------------------------------
// LLBuf
//-----------------------------------------------------------------------------

Pedigree::LLBuf::LLBuf( Pedigree & _ped ) :
	ped	      ( _ped  ) ,
	accValIsValid ( false ) ,
	propActive    ( false )
    {
    }


Pedigree::LLBuf::LLBuf( Pedigree & _ped, const LLBuf & rhs ) :
	ped	      ( _ped ) ,
	accVal	      ( rhs.accVal ) ,
	accValIsValid ( rhs.accValIsValid ) ,
	propActive    ( false )
    {
    gp_assert( ! rhs.propActive );
    }


Pedigree::LLBuf & Pedigree::LLBuf::operator=( const LLBuf & rhs )
    {
    accVal	  = rhs.accVal	      ;
    accValIsValid = rhs.accValIsValid ;
    propActive	  = false	      ;

    gp_assert( ! rhs.propActive );

    return *this;
    }


/// "Begin"
void Pedigree::LLBuf::startProposal()
    {
    gp_assert( ! propActive );
    propActive = true;
    propIsValid = false;
    }


/// "Commit"
/// It would be very strange if a proposal was accepted without first computing
/// the LL, but just in case, we'll handle it gracefully, so long as this method
/// is called prior to the parameters (rho, theta) being accepted.
void Pedigree::LLBuf::acceptProposal()
    {
    gp_assert( propActive );
    if ( propIsValid )
	accVal = propVal;
    else
	accVal = getHMM().getLogLikelihood();
    propActive = false;
    }


/// "Rollback"
void Pedigree::LLBuf::rejectProposal()
    {
    gp_assert( propActive );
    propActive = false;
    }


/// Get the log-likelihood of the last accepted parameters (rho,theta)
double Pedigree::LLBuf::getAcceptedVal()
    {
    if ( ! accValIsValid )
	{
	gp_assert( ! propActive ); // Maybe not, could the first time be during a proposal?
	accVal = getHMM().getLogLikelihood();
	accValIsValid = true;
	}
    return accVal;
    }


/// Get the log-likelihood of the proposed parameters (rho,theta)
double Pedigree::LLBuf::getProposedVal()
    {
    gp_assert( propActive );
    if ( ! propIsValid )
	{
	propVal = getHMM().getLogLikelihood();
	propIsValid = true;
	}
    return propVal;
    }



/// Invalidate the cache (proposed if a proposal is active, otherwise
/// last-accepted): inform it that the underlying parameters (rho or theta)
/// changed.  Ordinarily, this need not be called since the proposed cache is
/// already automatically invalid when a new proposal is started, and once
/// computed, the value is either accepted or rejected without further changes
/// to the parameters.

void Pedigree::LLBuf::parmsChanged()
    {
    if ( propActive )
	propIsValid = false;
    else
	accValIsValid = false;
    }



//-----------------------------------------------------------------------------
// Overridden virtuals from PedBase
//-----------------------------------------------------------------------------

unsigned int Pedigree::getMyNumber() const
    {
    return myNumber;
    }


unsigned int Pedigree::getIndex() const
    {
    return myNumber - 1;
    }


unsigned int Pedigree::getNumObs() const
    {
    return getNFounders();
    }



//-----------------------------------------------------------------------------
// setMyNumber()
//-----------------------------------------------------------------------------

void Pedigree::setMyNumber( unsigned int nv )
    {
    myNumber = nv;
    #if PED_HAS_OWN_PRNG
	//fprintf( stderr, "Seed %u (%s) -> %ld\n", nv, getId().c_str(), getOptions().getSeed()+nv ); //DEBUG
	rng.seed( getOptions().getSeed() + nv );
    #endif
    }



//-----------------------------------------------------------------------------
// setOptions(), very bad.
//-----------------------------------------------------------------------------

void Pedigree::setOptions( const AdmixOptions & opts )
    {
    optionsRef = &opts;
    }


//-----------------------------------------------------------------------------
//
/// These intitialisations should be in the constructor, when we have a derived
/// class.  See NOTE *2* in Pedigree.h.
///
/// We can't initialise some of the members prior to knowing, amongst other
/// things, the number of founders.  They should of course be in a separate
/// subclass, in which case they wouldn't be allocated until much later anyhow.
//
//-----------------------------------------------------------------------------

void Pedigree::InitialiseAdmixedStuff( const AdmixOptions & options )
    {

    // See getOptions() in Pedigree.h for an explanation of this.
    if ( optionsRef == 0 )
	optionsRef = &options;
    else
	gp_assert( optionsRef == &options );


    const PopIdx K = getK();
    for ( ProposalRingBuffer<ThetaType>::Iterator it = thetas.begin() ; it != thetas.end() ; ++it )
	{
	it->resize( getNTheta() );
	for ( ThetaIdx fIdx = getNTheta() ; fIdx-- != 0 ; )
	    (*it)[ fIdx ].resize( K );
	}


    pvectord Kzero;
    Kzero.resize( getK(), 0.0 );
    SumSoftmaxTheta.resize( getNTheta(), Kzero );


    // Set the initial values for Theta:
    SetUniformAdmixtureProps();


    ThetaTuner.SetParameters( step, 0.0001, 10.0, 0.44 );


    #if ! NON_GLOBAL_RHO_WORKS
	if ( ! options.isGlobalRho() )
	    throw std::runtime_error( "Only global-rho currently is supported." );
    #endif

    #if ! NOT_NEEDED_UNLESS_CONJUGATE_UPDATE
	if ( ! options.getNoConjugateUpdate() )
	    throw std::runtime_error( "Pedigrees (as yet) only support the no-conjugate-update option."
					" Re-run with no-conjugate-update turned on." );
    #endif

    #if ! NOT_NEEDED_FOR_FIXEDALLELEFREQ_EQ_1
	if ( ! options.getFixedAlleleFreqs() )
	    throw std::runtime_error( "Sorry, pedigrees are not (yet) compatible with fixedallelefreqs=0."
					" Re-run with fixedallelefreqs=1." );
    #endif

    #if ! SUPPORT_ASSOCIATION_TESTS
	if ( options.hasAnyAssociationTests() )
	    throw std::runtime_error( "Sorry, pedigrees are not (yet) compatible with association tests."
					" Re-run with some options turned off." );
    #endif


    if ( getSLoci().size() != getSLoci().getNComposite() )
	throw std::runtime_error( "Sorry, pedigrees are not (yet) compatible with composite loci." );


    // This was AdmixedIndividual::InitialiseSumIntensities():
    double init;
    if ( !options.isGlobalRho() && options.getIndAdmixHierIndicator())
	{//model with individual- or gamete-specific sumintensities

	//set prior mean as initial value for rho
	if ( options.getRhobetaShape() > 1 )
	    init = options.getRhoalpha() * options.getRhobetaRate() / (options.getRhobetaShape() - 1 );
	else
	    init = options.getRhoalpha() * options.getRhobetaRate() / options.getRhobetaShape() ;//conditional prior mean
	}
    else // no hierarchical model or globalrho prior
	{
	init = options.getRhoalpha() / options.getRhobeta();
	}


    for ( RhoType * it = rhos.begin() ; it != rhos.end() ; ++it )
	it->resize( options.isGlobalRho() ? 1 : NumGametes );

    // Initialise rho to 0 for unadmixed gametes
    for ( size_t idx = getCurRho().size() ; idx-- != 0 ; )
	getCurRho()[idx] = options.isAdmixed(idx) ? init : 0.0;

    sumlogrho.assign( getCurRho().size(), 0.0 );

    }



//-----------------------------------------------------------------------------
/// Overridden virtual from PedBase, here only for compatibility with the
/// individual objects.  The arguments are all ignored.  Returns the
/// log-likelihood of the current parameters (rho,theta) -- i.e. proposed if a
/// proposal is active, last accepted value otherwise.
//-----------------------------------------------------------------------------

double Pedigree::getLogLikelihood( const Options & , bool /*forceUpdate*/ , bool /*store*/ )
    {
    double rv;

    #if USE_SINGLE_POPULATION_SPECIAL_CASE
      if ( getK() == 1 )
	  rv = getLogLikelihoodOnePop();
      else
    #endif
	rv = llCache.getCurVal();

    return rv;
    }



//-----------------------------------------------------------------------------
// calcNInheritedByAffected() [private]
//
// Helper method for affected-only test computations [getNInheritedByAffected()].
//
// METHOD-NOTE *1*: We could just allocate ancStack of size
//	[getNMembers()-fIdx], then subtract the offset from member-indexes when
//	looking up, but this is easier and not too large.  We need to initialise
//	the inheritance of all founders other than the one in question (fIdx) to
//	"false", however, even if the ancestry-vector in question shows "k" at
//	that gamete, because we are only interested in gametes inherited from
//	_this_ founder with ancestry k, so initialising those to false (which is
//	done in the default constructor) frees us up to look up the parents of
//	any descendant nodes and increment the return-value if they inherit from
//	that population and are affected, without at that time checking from
//	which founder-gamete they inherited the gamete in question.
//
//-----------------------------------------------------------------------------

inline short Pedigree::calcNInheritedByAffected( PopIdx k, FounderIdx fIdx,
				const AncestryVector & av, const InheritanceVector & iv ) const
    {

    int rv = 0;

    #if DEBUG_CALC_M
	cerr << "calcNInherited[ped=" << getId() << "](k=" << k << ",fIdx=" << fIdx
	    << '(' << founderAt(fIdx).getOrgId() << ")," << av << ',' << iv << "):\n";
    #endif

    struct Ancestry
	{
	bool pAncestry;
	bool mAncestry;

	Ancestry( bool p, bool m ) : pAncestry(p), mAncestry(m) {}
	Ancestry() :
	    pAncestry(false),	// Initialise to false:
	    mAncestry(false) {} //  see METHOD-NOTE *1*

	bool hasAncOfType( const InheritanceVector::SegInd & si ) const
	    {
	    return (si == InheritanceVector::SI_PATERNAL) ? pAncestry : mAncestry;
	    }
	};


    Ancestry ancStack[ getNMembers() ]; // See METHOD-NOTE *1*


    if ( founderAt(fIdx).isHaploid() )
	{
	const PopIdx ancestry = av.at( fIdx, GT_SINGLE );

	// NOTE *X9*: Use the paternal slot to hold the ancestry in the case of
	// a single-gamete founder.  This will later be recognized when tracing
	// back inheritance of ancestry.
	ancStack[fIdx].pAncestry = (ancestry == k);
	}
    else
	{
	ancStack[fIdx].pAncestry = (av.at(fIdx,GT_PATERNAL) == k);
	ancStack[fIdx].mAncestry = (av.at(fIdx,GT_MATERNAL) == k);
	}

    #if 0

	const Ancestry & parentAncestry = ancStack[ fIdx ];

	for ( AdmixPedigree::ChConstIter it = childrenBegin() ; it != childrenEnd() ; ++it )
	    {
	    const Organism & child = *it;
	    const MemberIdx cIdx = child.getPIdx();
	    const InheritanceVector::SegInd si = isFemale() ? iv.maternal(cIdx) : iv.paternal(cIdx);
	    inheritsPop = (si == SI_PATERNAL) ? parentAncestry.pAncestry : parentAncestry.mAncestry;
	    if ( inheritsPop )
		{
		if ( isFemale() )
		    ancStack[cIdx].mAncestry = true;
		else
		    ancStack[cIdx].pAncestry = true;
		if ( child.isAffected() )
		    {
		    ++rv;
		    #if DEBUG_CALC_M
			cerr << "   path to affected-child " << child.getOrgId() << "; result so far: " << rv << '\n';
		    #endif
		    }
		}
	    Then recurse here...
	    }

    #else

	// Or do something like this...

	for ( MemberIdx cIdx = fIdx + 1 ; cIdx < getNMembers() ; ++cIdx )
	    {

	    const Organism & child	= memberAt( cIdx );
	    Ancestry &	     cAncestry	= ancStack[ cIdx ];

	    const Organism * const father = child.getFather();
	    if ( father != 0 )
		{
		const Ancestry & pAncestry = ancStack[ father->getPIdx() ];
		if ( father->isHaploid()				?
			    pAncestry.pAncestry				:
			    pAncestry.hasAncOfType( iv.paternal(cIdx) ) )
		    {
		    cAncestry.pAncestry = true;
		    if ( child.getOutcome() == Organism::OUTCOME_AFFECTED ) // missing counted as unaffected
			{
			++rv;
			#if DEBUG_CALC_M
			    cerr << "   pat-path to affected-child " << child.getOrgId() << "; result so far: " << rv << '\n';
			#endif
			}
		    }
		}


	    // If we have already found a transmitted gamete on the father's
	    // side above, we would only need to continue to look for one on the
	    // mother's side in the case of incest in a multi-generational
	    // pedigree, which might result in a cycle in the pedigree (treated
	    // as a non-directed graph).  That is, we could "continue" to the
	    // bottom of the loop after having incremented rv above; but we
	    // continue checking here anyhow, since such things do sometimes
	    // happen.


	    const Organism * const mother = child.getMother();
	    if ( mother != 0 )
		{
		const Ancestry & mAncestry = ancStack[ mother->getPIdx() ];
		if ( mother->isHaploid()				?
			    mAncestry.pAncestry				:
			    mAncestry.hasAncOfType( iv.maternal(cIdx) ) )
		    {
		    cAncestry.mAncestry = true;
		    if ( child.getOutcome() == Organism::OUTCOME_AFFECTED ) // missing counted as unaffected
			{
			++rv;
			#if DEBUG_CALC_M
			    cerr << "   mat-path to affected-child " << child.getOrgId() << "; result so far: " << rv << '\n';
			#endif
			}
		    }
		}

	    }

    #endif

    #if DEBUG_CALC_M
	cerr << "   final result: " << rv << '\n';
    #endif

    return rv;

    }



//-----------------------------------------------------------------------------
// getNInheritedByAffected() [private]
//
// Helper method for affected-only test computations [accumAOScore()].
//
/// Given a particular hidden-state (i.e. an AncestryVector and
/// InheritanceVector for the pedigree), compute the number of affected
/// offspring of a given founder (@a fIdx) who inherit ancestry in population
/// @a k from founder @a fIdx.
///
/// Since our current assumption is that a two-gamete founder's ancestry is the
/// same for both gametes, XXX.
//-----------------------------------------------------------------------------

inline int Pedigree::getNInheritedByAffected( PopIdx k, FounderIdx fIdx, const AncestryVector & av, const InheritanceVector & iv ) const
    {

    // Is this optimization worth having?  Does it come up in real datasets?
    #if 0
	if ( getNAffected() == 0 )
	    return 0; // *** RETURN HERE ***
    #endif


    #if USE_AO_CACHE
	if ( aoCache == 0 )
	    {

	    aoCache = new TwoDimArray<short,PopIdx,FounderIdx>( getK(), getNFounders() );

	    for ( PopIdx i_k = getK() ; i_k-- != 0 ; )
		for ( FounderIdx i_f = getNFounders() ; i_f-- != 0 ; )
		    (*aoCache).get(i_k,i_f) = calcNInheritedByAffected( i_k, i_f, av, iv );

	    }


	return (*aoCache).get( k, fIdx );

    #else // USE_AO_CACHE

	return calcNInheritedByAffected( k, fIdx, av, iv );

    #endif

    }



//-----------------------------------------------------------------------------
// accumAAScore()
//-----------------------------------------------------------------------------


// Static helper:
static inline double aa_score( int nFromK, double nAffOver2, double negNAffOver2, double nAffOver4, double mu )
    {
    if ( nFromK == 0 )
	return (mu * negNAffOver2);
    else if ( nFromK == 1 )
	return nAffOver4 + (negNAffOver2 * mu);
    else
	{
	AGGRESSIVE_ONLY( if ( nFromK != 2 ) throw std::logic_error( "invalid nFromK" ); )
	return nAffOver2 * (1.0 - mu);
	}
    }


// The method:
void Pedigree::accumAOScore( AffectedsOnlyTest & aoTest ) const
    {

    // Is this optimization worth having?  Does it come up in real datasets?
    // Note that it affects the unrelated-individual special-case, so test that
    // for is-affected if this is removed.
    if ( getNAffected() == 0 )
	return; // *** RETURN HERE ***


    const MemberIdx nAff = getNAffNonFndr();

    // Precompute invariants outside the loop:
    const double nAffOver2	= double(nAff) / 2;
    const double negNAffOver2	= - nAffOver2;
    const double nAffOver4	= double(nAff) / 4;
    const double nAffOver64	= double(nAff) / 64;
    const double nAffSqOver8	= double(nAff * nAff) / 8;

    const SimpleLocusArray & loci = getSLoci();

    // If K==2, only evaluate for 1 population, otherwise for all of them.
    // Existing code uses k==1 (not 0), yet treats it within the AffectedsOnlyTest object as "k==0", so:
    const PopIdx k_begin = (getK() == 2) ? 1 : 0;

    for ( SLocIdxType t = 0 ; t < loci.size() ; ++t )
	{

	#if DEBUG_CALC_M
	    cerr << "\n=== AO-debug: locus " << t << " (" << loci[t].getName() << ")\n";
	#endif

	// If K==2, only evaluate for 1 population, otherwise for all of them.
	// Existing code uses k==1 (not 0), so:
	for ( PopIdx k = k_begin ; k < getK() ; ++k )
	    {

	    #define DEBUG_AOTEST 0
	    #if DEBUG_AOTEST
		fprintf( stderr, "ped:%s(%d) t:%zd k:%zd\n", getId().c_str(), getMyNumber(), t, k );
	    #endif

	    double scoreAvg   = 0.0;
	    double scoreSqAvg = 0.0;
	    double infoAvg    = 0.0;

	    const HiddenStateSpace & hss	    = getStateProbs( t );
	    const cvector<double> &  condStateProbs = getHMM().getCondStateProbsAtLocus( t );

	    for ( HiddenStateSpace::Iterator stIt( hss ) ; stIt ; ++stIt )
		{

		double stScore;
		double stInfo;

		// Special case for a single unrelated individual: this can be
		// removed should we model it as two single-gamete parents.
		if ( getNMembers() == 1 )
		    {

		    const double mu = getCurTheta()[ 0 ][ k ]; // Could move outside the inner HSS loop.

		    stScore = -mu;

		    //gp_assert( stId->getIV() == InheritanceVector::null_IV() );
		    const AncestryVector & av = stIt.getAV();

		    // Assert: the sole member (unrelated individual) is not a
		    // single-gamete-founder, but rather has two gametes with
		    // observed genotyped data.  This is not necessary, but it's
		    // not clear (to me) why we would have in the input dataset
		    // an unrelated individual with no genotyped data.  This can
		    // be removed if the score computation below is adjusted
		    // accordingly.
		    gp_assert( av.size() == 2 );

		    if ( av.at(0) == k ) // Paternal ancestry
			stScore += 0.5;
		    if ( av.at(1) == k ) // Maternal ancestry
			stScore += 0.5;

		    stInfo = 0.5 * (1 - mu) * mu;

		    }

		else
		    {

		    stScore = 0.0;
		    stInfo  = 0.0;

		    for ( FounderIdx fIdx = getNFounders() ; fIdx-- != 0 ; )
			{

			const AncestryVector & av = stIt.getAV();

			const Member & founder = founderAt( fIdx );

			#if DEBUG_CALC_M
			  cerr << "F-idx " << fIdx << " (ID " << founder.getOrgId() <<
			    "): is " << (founder.isHaploid()?"":"not ") << "haploid, is "
			    << (av.isHetrozygousForPop(fIdx,k)?"":"not ") << "heterozygous in "
			    << av << " WRT pop " << k << ".\n";
			#endif

			const double mu = getCurTheta()[ fIdx ][ k ];

			//-------------------------------------------------------
			// Special case for founders modeled as a single-gamete:
			// it may be possible to use the 2-gamete computation
			// for this, if we can verify that it simplified
			// correctly to the same computation here.
			//
			// If not, consider passing the is-haploid flag into
			// AncestryVector's nCopiesFromKAtFounder() and
			// isHetrozygousForPop() to avoid their looking it up (or
			// just having them assume that it is diploid if we only
			// use those methods for that case).
			//
			// Alternatively, just get both ancestries once here,
			// keep them, and do the isHetro and nCopies
			// computations here.
			//-------------------------------------------------------

			if ( founder.isHaploid() )
			    {
			    stScore += 0.5 * (((av.at(fIdx,GT_SINGLE) == k) ? 1 : 0) - mu);
			    stInfo += 0.25 * (1 - mu) * mu;
			    }
			else
			    {
			    // AncestryVector::isHetrozygousForPop() will return
			    // false for haploid founders (i.e. those that we model
			    // as a single gamete), so we need no special-case for them here.
			    if ( av.isHetrozygousForPop( fIdx, k ) )
				{
				const int m = getNInheritedByAffected( k, fIdx, av, stIt.getIV() );
				stScore += 0.25 * (m - nAffOver2);
				stInfo += nAffOver64;
				#if DEBUG_AOTEST
				    fprintf( stderr, " hetro-m:%d; score:%.12lf; info:%.12lf", m, stScore, stInfo );
				#endif
				}

			    const int nFromK = av.nCopiesFromKAtFounder( fIdx, k );

			    #if DEBUG_AOTEST
				fprintf( stderr, "	fIdx:%zd st:%zd nFromK:%d hetro:%s nAff:%zd",
					fIdx, stIt.getNon0Index(), nFromK,
					av.isHetrozygousForPop(fIdx,k) ? "yes" : "no", nAff );
			    #endif

			    stScore += aa_score( nFromK, nAffOver2, negNAffOver2, nAffOver4, mu );
			    stInfo  += nAffSqOver8 * (1 - mu) * mu;

			    #if DEBUG_AOTEST
				fprintf( stderr, " mu:%.12lf; stScore:%.12lf; stInfo:%.12lf\n", mu, stScore, stInfo );
			    #endif
			    }
			}

		    }


		const double stWeight = condStateProbs[ stIt.getNon0Index() ];
		const double wtScore = stScore * stWeight;
		scoreAvg += wtScore;
		scoreSqAvg += wtScore * stScore;
		infoAvg += stInfo * stWeight;

		#if DEBUG_AOTEST
		    fprintf( stderr, "-> stWeight=%.12lf wtScore=%.12lf scoreAvg=%.12lf\n", stWeight, wtScore, scoreAvg );
		#endif

		}

	    const double condVarianceOfU = scoreSqAvg - (scoreAvg * scoreAvg);

	    // In practice, the aoTest object that is passed in is effectively a
	    // shared global object.  We should probably find a way to address
	    // this more cleanly and efficiently.  We can disable this critical
	    // section by marking as critical-section the entire call to
	    // el.UpdateScores() inside the parallelized for-loop in
	    // AdmixIndividualCollection::HMMUpdates().
	    #if defined(_OPENMP) && PARALLELIZE_PEDIGREE_LOOP
		#pragma omp critical
	    #endif
		{ // BEGIN SCOPE BLOCK
		aoTest.getAffectedsScore    ( t, k-k_begin ) += scoreAvg	;
		aoTest.getAffectedsVarScore ( t, k-k_begin ) += condVarianceOfU ;
		aoTest.getAffectedsInfo	    ( t, k-k_begin ) += infoAvg		;
		} // END SCOPE BLOCK

	    #if DEBUG_AOTEST
		fprintf( stderr, "==> cum-score: %.12lf cum-var=%.12lf cum-info=%.12lf\n\n",
		    aoTest.getAffectedsScore(t,k-k_begin)	 ,
		    aoTest.getAffectedsVarScore(t,k-k_begin) ,
		    aoTest.getAffectedsInfo(t,k-k_begin)	 );
	    #endif
	    }

	}

    }



//-----------------------------------------------------------------------------
// HMMIsBad()
/// Overridden from PedBase, but null implementation because only applies to
/// one-HMM-per-chromosome.
//-----------------------------------------------------------------------------

void Pedigree::HMMIsBad( bool /*loglikisbad*/ )
    {
    }



//-----------------------------------------------------------------------------
// getTPC(), getHMM()
//-----------------------------------------------------------------------------

TransProbCache & Pedigree::getTPC() const
    {
    #if NON_GLOBAL_RHO_WORKS
	#error Which rho to use here?
	// Issue yet to resolve: rho indexing.  On what is it indexed for
	// pedigrees?  Which element to pass into TPC constructor?  For the
	// moment, we are restrictng pedigrees to globalrho=1, so it's not
	// relevant.
    #else
	if ( tpCache == 0 )
	    tpCache = new TransProbCache( *this, getGlobalRhoVal(), getCurTheta() );
    #endif
    return *tpCache;
    }


/// Allocate the HiddenMarkovModel for the Pedigree object (only if have one HMM
/// for each pedigree rather than one per chromosome).

HiddenMarkovModel & Pedigree::getHMM() const
    {
    if ( hmm == 0 )
	hmm = new HiddenMarkovModel( *this, getTPC(), &getCurTheta() );
    return *hmm;
    }

void Pedigree::freeHMM() const
    {
    delete tpCache;
    tpCache = 0;
    delete hmm;
    hmm = 0;
    }



//-----------------------------------------------------------------------------
//
// SampleTheta()
//
/// Samples admixture proportions.
///
/// @param RW true for a random-walk proposal, false for a conjugate proposal
///	    NOTE *1*: For Pedigree however, we are ignoring @a RW (treating it
///	    as true for every call) so as to avoid conjugate-update, thus
///	    avoiding calling ResetSufficientStats().  There is also a runtime
///	    option, no-conjugate-update, which should be specified every time
///	    that pedigrees are used.
//
//-----------------------------------------------------------------------------

void Pedigree::SampleTheta(
		 int				    iteration		,
		 double *			    SumLogTheta		, ///< Indexed on 0 <= k < K
		 const bclib::DataMatrix *	    Outcome		,
		 const DataType *		    OutcomeType		,
		 const std::vector<double> &	    lambda		,
		 int				    NumCovariates	,
		 bclib::DataMatrix *		    Covariates		,
		 const std::vector<const double*> & beta		,
		 const PopThetaType &		    poptheta		,
		 const AdmixOptions &		    options		,
		 const AlphaType &		    alpha		,
		 double				    /*DInvLink*/	,
		 double				    /*dispersion*/	,
		 CopyNumberAssocTest &		    /*ancestryAssocTest*/,
		 bool				    /*RW*/		,
		 bool				    anneal		)
    {

    // NOTE *1*: We ignore the value of RW, use true every time (see comments above):
    #define RW true

    const PopIdx K = getK();

    double logpratio = 0.0;
    try
	{
	if ( RW )
	    {
	    NumberOfUpdates++;
	    logpratio += ProposeThetaWithRandomWalk( alpha );
	    }
	  #if NEEDED_ONLY_FOR_CONJUGATE_UPDATE
	    else
		ProposeTheta( options, alpha, SumLocusAncestry, SumLocusAncestry_X );
	    #error this must be converted
	  #endif
	}
    catch ( string & s )
	{
	throw std::runtime_error(
	    "Error generating proposed pedigree admixture proportions:\n  " + s );
	}

    // Calculate Metropolis acceptance probability ratio for proposal theta
    if ( (! options.getTestForAdmixtureAssociation()) && (getMyNumber() < Outcome->nCols()) )
	{
	const int NumOutcomes = Outcome->nCols();
	for ( int k = 0; k < NumOutcomes; k++ )
	    {
	    const RegressionType RegType = (OutcomeType[k] == Binary) ? Logistic : Linear;
	    logpratio += LogAcceptanceRatioForRegressionModel( RegType, options.isRandomMatingModel(), K,
						NumCovariates, Covariates, beta[k],
						Outcome->get( getIndex(), k ), poptheta, lambda[k] );
	    }
	}

    // Accept or reject proposed value - if conjugate update and no regression
    // model, proposal will be accepted because logpratio = 0
    Accept_Reject_Theta( logpratio, options.isRandomMatingModel(), RW );

    // Update the value of admixture proportions used in the regression model, but
    // only if pedigree is a single individual:
    if ( (options.getNumberOfOutcomes() > 0) && (getNMembers() == 1) )
	updateAdmixtureForRegression( NumCovariates, poptheta, options.isRandomMatingModel(), Covariates );


    // See NOTE *1* in list-of-things-dont-agree here.
    if ( (! anneal) && (iteration > options.getBurnIn()) )
	{
	// Accumulate sums in softmax basis for calculation of posterior means

	double a[ K ]; // Or could use: pvectord a( K );

	for ( ThetaIdx tIdx = 0 ; tIdx < getNTheta() ; ++tIdx )
	    {
	    getCurTheta()[tIdx].inv_softmax_gt0( a, 0.0 ); // Equivalent: getCurTheta()[tIdx].inv_softmax( a, gt_0 );
	    SumSoftmaxTheta[tIdx] += a;
	    }
	}


    if ( ! IAmUnderTest )
	{
	// See NOTE *1* in list-of-things-dont-agree here.
	for ( ThetaIdx tIdx = 0 ; tIdx < getNTheta() ; ++tIdx )
	    {
	    //SumLogTheta.addAfterFunction( getCurTheta()[tIdx], log );
	    const ThetaElType & th = getCurTheta()[ tIdx ];

	    // The SumLogTheta object that is passed in is effectively a shared
	    // global object which gets modified here, so we define a critical
	    // section.  We should probably find a way to address this more
	    // cleanly and efficiently, perhaps passing the partial-sum of
	    // log(theta) for this iteration of pedigree out, and let the
	    // parallelized loop accumulate the final values, using OpenMP
	    // reduction().  Alternately, use some combination of
	    // copyin/copyout/firstprivate/lastprivate.  If the entire call to
	    // SampleTheta() is inside a critical section (every place that it
	    // is called from a parallel section), this does not need to be
	    // protected here.
	    #if defined(_OPENMP) && PARALLELIZE_PEDIGREE_LOOP && (! SAMPLE_THETA_CALL_CRITICAL)
		#pragma omp critical
	    #endif
		{ // BEGIN SCOPE BLOCK
		for ( size_t k = th.size() ; k-- != 0 ; )
		    SumLogTheta[k] += log( th[k] );
		} // END SCOPE BLOCK
	    }

	#if NOT_NEEDED_FOR_FIRST_VERSION
	    //increment B using new Admixture Props
	    //Xcov is a vector of admixture props as covariates as in UpdateScoreForAncestry
	    if(iteration >= options.getBurnIn() && options.getTestForLinkageWithAncestry()){
	      double* admixtureCovars = new double[K-1];
	      for(int t = 0; t < K-1; ++t)admixtureCovars[t] = Covariates->get(getIndex(), Covariates->nCols()-K+1+t);
	       ancestryAssocTest.UpdateB(DInvLink, dispersion, admixtureCovars);
	      delete[] admixtureCovars;
	    }
	#endif
	}


    #undef RW // NOTE *1*

    }



//-----------------------------------------------------------------------------
// LogAcceptanceRatioForRegressionModel()
//
/// Returns log of ratio of likelihoods of new and old values of population
/// admixture in regression models.  individual admixture theta is
/// standardized about the mean poptheta calculated during burn-in.
//-----------------------------------------------------------------------------

double Pedigree::LogAcceptanceRatioForRegressionModel( RegressionType /*RegType*/, bool /*RandomMatingModel*/,
							PopIdx /*K*/, int /*NumCovariates*/,
							const bclib::DataMatrix * /*Covariates*/, const double * /*beta*/,
							const double /*Outcome*/, const PopThetaType & /*poptheta*/,
							double /*lambda*/ ) const
    {

    #if NOT_NEEDED_FOR_FIRST_VERSION

	double logprobratio = 0.0;
	double XBeta	= 0.0;
	double currentXBeta = 0.0;

	vector<double> avgtheta( K );
	avgtheta[0] = 0.0;

	vector<double> currentavgtheta( K );
	currentavgtheta[0] = 0.0;

	// Resume here...
	if( RandomMatingModel && NumGametes==2)
	  for(int k = 1;k < NumHiddenStates; ++k){
	    avgtheta[k] = (ThetaProposal[k] + ThetaProposal[k+NumHiddenStates]) / 2.0 - poptheta[k];
	    currentavgtheta[k] = (Theta[k] + Theta[k + NumHiddenStates ])/ 2.0 - poptheta[k];
	  }
	else
	  for(int k = 1;k < NumHiddenStates; ++k){
	    avgtheta[k] = ThetaProposal[k]  - poptheta[k];
	    currentavgtheta[k] = Theta[k]  - poptheta[k];
	  }

	for( int jj = 0; jj < NumCovariates - NumHiddenStates + 1; jj++ ){
	  XBeta += Covariates->get( getIndex(), jj ) * beta[jj];
	  currentXBeta += Covariates->get( getIndex(), jj ) * beta[jj];
	}
	for( int k = 1; k < NumHiddenStates; k++ ){
	  XBeta += avgtheta[ k ] * beta[NumCovariates - NumHiddenStates + k ];
	  currentXBeta += currentavgtheta[ k ] * beta[NumCovariates - NumHiddenStates + k];
	}
	if(RegType == Linear) {
	  logprobratio = 0.5 * lambda * (XBeta - currentXBeta) * (Outcome + Outcome - currentXBeta - XBeta);
	  // ( currentXBeta - Outcome )^2 - ( XBeta - Outcome )^2 factorized
	}
	else if(RegType == Logistic) {
	  if( Outcome == 1 ) {
	    logprobratio =  log( ( 1.0 + exp( -currentXBeta ) ) / ( 1.0 + exp( -XBeta ) ) );
	  } else {
	    logprobratio =  log( ( 1.0 + exp( currentXBeta ) ) / ( 1.0 + exp( XBeta ) ) );
	  }
	}
	return( logprobratio );
    #else
	throw std::runtime_error( "Not implemented" );
    #endif
    }



//-----------------------------------------------------------------------------
// updateAdmixtureForRegression() [private]
/// Update individual admixture values (mean of both gametes) used in the
/// regression model.  Will use the proposed theta if there is an active
/// proposal, or the last accepted values otherwise.
//-----------------------------------------------------------------------------

void Pedigree::updateAdmixtureForRegression( int NumCovariates,
                                               const PopThetaType & poptheta, bool /*RandomMatingModel*/,
					       bclib::DataMatrix * Covariates )
    {
    double avgtheta[ getNTheta() ];


    for ( PopIdx k = getK() ; k-- != 0 ; )
	{
	double sum = 0.0;
	for ( ThetaIdx tIdx = getNTheta() ; tIdx-- != 0 ; )
	    sum += getCurTheta()[ tIdx ][ k ];
	avgtheta[ k ] = sum / getNTheta();
	}


    const PopIdx K = getK();
    const size_t NumberOfInputCovariates = NumCovariates - K;


    // The Covariates object that is passed in is effectively a shared global
    // object which gets modified here, so we define a critical section.  We
    // should probably find a way to address this more cleanly and efficiently,
    // perhaps passing the partial-sum of log(theta) for this iteration of
    // pedigree out, and let the parallelized loop accumulate the final values,
    // using OpenMP reduction().  Alternately, use some combination of
    // copyin/copyout/firstprivate/lastprivate.  If the entire call to
    // SampleTheta() is inside a critical section (every place that it is called
    // from a parallel section), this does not need to be protected here.
    #if defined(_OPENMP) && PARALLELIZE_PEDIGREE_LOOP && (! SAMPLE_THETA_CALL_CRITICAL)
	#pragma omp critical
    #endif
	{ // BEGIN SCOPE BLOCK
	// Seems strange, why does this loop start at 1?
	for( PopIdx k = 1 , covIdx = NumberOfInputCovariates ; k < K ; k++, covIdx++ )
	    Covariates->set( getIndex(), covIdx, avgtheta[ k ] - poptheta[ k ] );
	} // END SCOPE BLOCK
    }



//-----------------------------------------------------------------------------
// ProposeThetaWithRandomWalk() [private]
//-----------------------------------------------------------------------------

double Pedigree::ProposeThetaWithRandomWalk( const AlphaType & alpha )
    {
    double LogLikelihoodRatio = 0.0;
    double LogPriorRatio      = 0.0;


    //-------------------------------------------------------------------------
    // Generate proposals
    //-------------------------------------------------------------------------

    const PopIdx K = getK();

    // NOTE *M1*: Maybe this should be be allocated once for the class?
    pvectord a( K );

    ThetaType & theta_proposal = startThetaProposal();
    ThetaType & theta	       = thetas.getAccepted();

    for ( ThetaIdx tIdx = 0 ; tIdx < getNTheta() ; ++tIdx )
	{

	if ( IS_ADMIXED() )
	    {
	    const ThetaElType & th	= theta[ tIdx ];
	    ThetaElType &	th_prop = theta_proposal[ tIdx ];

	    th.inv_softmax_gt0( a, SOFTMAX_0_FLAG );

	    DEBUG_TH_PROP(
		fprintf( stderr, "DBG-TH-PR-1: %d %s %zu", getMyNumber(), getId().c_str(), tIdx );
		for ( PopIdx k = 0 ; k < K ; ++k )
		    fprintf( stderr, " %.9lf", th[k] );
		putc( '\n', stderr );
		fprintf( stderr, "DBG-TH-PR-2: %d %s %zu", getMyNumber(), getId().c_str(), tIdx );
		for ( PopIdx k = 0 ; k < K ; ++k )
		    fprintf( stderr, " %.9lf", a[k] );
		putc( '\n', stderr );

		fprintf( stderr, "PRNG POISON TEST: %.12lf\n", RNG_UNIFORM() );
	    ) // end DEBUG_TH_PROP()

	    // Random walk step - on all elements of array a
	    for ( PopIdx k = 0 ; k < K ; ++k )
		if ( a[k] != SOFTMAX_0_FLAG ) // Equivalent: if ( th[k] > 0.0 )
		    a[k] = RNG_NORMAL( a[k], step );

	    DEBUG_TH_PROP(
		fprintf( stderr, "DBG-TH-PR-3: %d %s %zu", getMyNumber(), getId().c_str(), tIdx );
		for ( PopIdx k = 0 ; k < K ; ++k )
		    fprintf( stderr, " %.9lf", a[k] );
		putc( '\n', stderr );
	    ) // end DEBUG_TH_PROP()

	    // Reverse transformation from numbers on real line to proportions
	    a.softmax( th_prop, pvectord::not_equal_to(SOFTMAX_0_FLAG), 0.0 );

	    DEBUG_TH_PROP(
		fprintf( stderr, "DBG-TH-PR-4: %d %s %zu", getMyNumber(), getId().c_str(), tIdx );
		for ( PopIdx k = 0 ; k < K ; ++k )
		    fprintf( stderr, " %.9lf", th_prop[k] );
		putc( '\n', stderr );
	    ) // end DEBUG_TH_PROP()

	    // Compute contribution of this founder to log prior ratio
	    // Prior densities must be evaluated in softmax basis
	    for ( PopIdx k = 0; k < K; ++k )
		if ( th[k] > 0.0 ) // Equivalent: if (a[k] != SOFTMAX_0_FLAG)
		    LogPriorRatio += alpha[0][k] * (log(th_prop[k]) - log(th[k])); // Equivalent: log(th_prop[k]/th[k])
		else
		    th_prop[k] = th[k];
	    }

	else // IS_ADMIXED()
	    {
	    #if ALL_FOUNDERS_ARE_ADMIXED
		throw std::logic_error( "Logic error. All founders are admixed." );
	    #else
		Theta[tIdx] = ThetaProposal[tIdx];
	    #endif
	    }

	} // End loop over Theta's elements (i.e. founders)

    DEBUG_TH_PROP( dumpTheta( "New proposed theta" ); )

    // Get log likelihood at current parameter values - do not force update, store result of update
    LogLikelihoodRatio -= llCache.getAcceptedVal();

    // Get log likelihood at proposal theta and current rho and accumulate
    LogLikelihoodRatio += llCache.getProposedVal();

    DEBUG_TH_PROP( fprintf( stderr, "Accept-LLR: %.9lf  Propose-LLR: %.9lf  Ratio: %.9lf  PriorRat: %.9lf\n",
	llCache.getAcceptedVal(), llCache.getProposedVal(), LogLikelihoodRatio, LogPriorRatio ); )

    return LogLikelihoodRatio + LogPriorRatio; // Log ratio of full conditionals
    }



//-----------------------------------------------------------------------------
// Accept_Reject_Theta() [private]
//
/// Metropolis update for admixture proportions theta, taking log of acceptance
/// probability ratio as argument uses log ratio because this is easier than
/// ratio to calculate for linear regression model if conjugate update and no
/// regression model, logpratio remains set to 0, so all proposals are accepted
//
//-----------------------------------------------------------------------------

void Pedigree::Accept_Reject_Theta( double logpratio, bool /*isRandomMatingModel*/, bool isRandomWalk )
    {

    const ThetaType & theta	    = thetas.getAccepted();
    const ThetaType & thetaProposal = thetas.getProposed();


    // Loop over populations: if any element of proposed Dirichlet parameter
    // vector is too small, reject update without test step
    bool should_do_test = true;
    for ( ThetaIdx tIdx = getNTheta() ; should_do_test && (tIdx-- != 0) ; )
	{
	const ThetaElType & th	    = theta	   [ tIdx ];
	const ThetaElType & th_prop = thetaProposal[ tIdx ];

	for ( PopIdx k = getK() ; should_do_test && (k-- != 0) ; )
	    if ( (th[k] > 0.0) && (th_prop[k] < 0.0001) )
		should_do_test = false;
	}


    bool accept = false;
    double AccProb = 1.0;

    if ( should_do_test )
	{
	// Generic Metropolis step
	if ( logpratio < 0 )
	    {
	    AccProb = exp( logpratio );
	    if ( RNG_UNIFORM() < AccProb )
		accept = true;
	    }
	else
	    accept = true;
	}


    DEBUG_TH_PROP( fprintf( stderr, "Accept: %s\n", accept ? "yes" : "no" ); )

    if ( accept )
	{
	// Set proposed values as new values
	acceptThetaProposal();

	gp_assert( isRandomWalk );

	// Copied comment:
	// Conjugate update of parameters invalidates both HMM forward probs and stored loglikelihood.
	}
    else
	{
	rejectThetaProposal();

	// Copied comment:
	// If RW proposal is rejected, loglikelihood.HMMisOK is already set to
	// false, and stored log-likelihood is still valid.
	}


    if ( isRandomWalk )
	{
	// Update step size in tuner object every w updates
	if ( (NumberOfUpdates % w) == 0 )
	    step = ThetaTuner.UpdateStepSize( AccProb );
	}

    }



//-----------------------------------------------------------------------------
// storeLogLikelihood()
//
/// Overridden from PedBase; deprecated in the case of pedigrees, which
/// encapsulate the proposal-caching logic.  We need to preserve it since it is
/// called externally from PopAdmix.
//-----------------------------------------------------------------------------

void Pedigree::storeLogLikelihood( bool /*setHMMAsOK*/ )
    {
    llCache.acceptProposal();

    if ( thetas.proposalActive() )
	thetas.acceptProposal();

    if ( rhos.proposalActive() )
	rhos.acceptProposal();
    }



//--------------------------------------------------------------------------
// Public rho-proposal methods, overridden from PedBase, called from PopAdmix,
// ignored for individuals, which store rho in the Loci.
//--------------------------------------------------------------------------

#if AGGRESSIVE_RANGE_CHECK
    double & Pedigree::getGlobalRhoVal()
	{
	gp_assert( getOptions().isGlobalRho() );
	return getCurRho()[0];
	}
#endif


// setRho()
// Virtual, overridden from PedBase.  Only makes sense if each pedigree keeps its own HMM and TPC.
void Pedigree::setRho( double nv )
    {
    getGlobalRhoVal() = nv;
    getTPC().setRho( nv );
    getHMM().transProbsChanged();
    llCache.parmsChanged();
    }

void Pedigree::startRhoProposal()
    {
    llCache.startProposal();
    rhos.startNewProposal();
    }

void Pedigree::acceptRhoProposal()
    {
    llCache.acceptProposal();
    rhos.acceptProposal();
    }

void Pedigree::rejectRhoProposal()
    {
    llCache.rejectProposal();
    rhos.rejectProposal();
    getTPC().setRho( getGlobalRhoVal() );
    getHMM().transProbsChanged();
    }



//--------------------------------------------------------------------
// Combined Theta-proposal methods
//--------------------------------------------------------------------

/// Allocates new theta on the ring-buffer for proposal and starts a new
/// logLikelihood proposal.
Pedigree::ThetaType & Pedigree::startThetaProposal()
    {
    llCache.startProposal();
    ThetaType & rv = thetas.startNewProposal();
    getHMM().setTheta( &rv );
    return rv;
    }

/// Reject Theta Proposal.
void Pedigree::rejectThetaProposal()
    {
    llCache.rejectProposal();
    thetas.rejectProposal();
    getHMM().setTheta( &getCurTheta() );
    }

/// Accept Theta Proposal.
void Pedigree::acceptThetaProposal()
    {
    llCache.acceptProposal();
    thetas.acceptProposal();
    }



//*****************************************************************************
//***************************  Log-Likelihoods  *******************************
//*****************************************************************************

//-----------------------------------------------------------------------------
// getLogLikelihoodAtPosteriorMeans()
//-----------------------------------------------------------------------------

double Pedigree::getLogLikelihoodAtPosteriorMeans( const Options & options )
    {

    double rv;

    // Should set allele freqs also to posterior means, and recalculate prob
    // genotypes at these freqs before calling getloglikelihood

    #if USE_SINGLE_POPULATION_SPECIAL_CASE
      if ( getK() == 1 )
	rv = getLogLikelihoodOnePop();
      else
    #endif
	{
	ThetaType & thetabar = startThetaProposal();

	// AdmixedIndividual's updateHMMInputs() only changes rho if not global-rho -- WHY???
	#if NON_GLOBAL_RHO_WORKS
	    if ( ! options.isGlobalRho() )
		{
		RhoType & rhobar = rhos.startNewProposal();
		gp_assert_eq( sumlogrho.size(), rhobar.size() );
		for ( unsigned int i = sumlogrho.size() ; i-- != 0 ; )
		    rhobar[i] = exp( sumlogrho[i] / double(options.getTotalSamples() - options.getBurnIn()) );
		#warning is/should this already done by rhos.rejectProposal()?
		getTPC().setRho( getGlobalRhoVal() );
		}
	#endif

	const double scale_factor = options.getTotalSamples() - options.getBurnIn();

	ThetaType scaledSSMT( SumSoftmaxTheta );
	scaledSSMT /= scale_factor;

	// Apply softmax transformation to obtain thetabar
	for ( ThetaIdx tIdx = 0 ; tIdx < getNTheta() ; ++tIdx )
	    scaledSSMT[tIdx].softmax_gt0( thetabar[tIdx], 0.0 ); // Equivalent: scaledSSMT[tIdx].softmax( thetabar[tIdx], gt_0 );

	#if 0
	    llCache.parmsChanged();	   // No need since starting the proposal automatically invalidates
	    getHMM().setTheta( thetabar ); // No need since we used a cache-proposal
	    getHMM().thetaChanged();	   // No need since we used a cache-proposal
	#endif

	rv = llCache.getProposedVal(); // getHMM().getLogLikelihood(), but caches

	#if NON_GLOBAL_RHO_WORKS
	    if ( ! options.isGlobalRho() )
		{
		rhos.rejectProposal();
		#warning is/should this already done by rhos.rejectProposal()?
		getTPC().setRho( getGlobalRhoVal() );
		}
	#endif

	rejectThetaProposal(); // We were never going to accept it anyhow.
	}

    return rv;
    }



//-----------------------------------------------------------------------------
// getLogLikelihoodOnePop()
//
// Convenient for a single population as no arguments required
// Not as yet converted to pedigrees.  The "GetGenotypeProbs()" is the emission
// probabilities for individuals.
//-----------------------------------------------------------------------------

#if USE_SINGLE_POPULATION_SPECIAL_CASE

    #error This needs, at the very least, to have a look-over, and probably more than that.

    double Pedigree::getLogLikelihoodOnePop() const
	{
	double rv = 0.0;

	// Loop over composite loci:
	for ( unsigned int j = 0; j < getGenome().GetNumberOfCompositeLoci(); j++ )
	    {
	    if ( ! GenotypeIsMissing(j) )
		{
		double Prob; // One pop so 1x1 array
		(*Loci)(j)->GetGenotypeProbs( &Prob, getPossibleHapPairs(j), false );
		rv += log( Prob );
		}
	    }

	return rv;
	}

#endif



//-----------------------------------------------------------------------------
// SetUniformAdmixtureProps() [protected]
/// Set the initial values of theta to a uniform distribution.  Operates on the
/// "current" theta.
//-----------------------------------------------------------------------------

void Pedigree::SetUniformAdmixtureProps()
    {
    const size_t K	    = getK();
    const double one_over_K = 1.0 / K;

    for ( ThetaIdx fIdx = getNTheta() ; fIdx-- != 0 ; )
	{
	ThetaElType & thetaEl = getCurTheta()[ fIdx ];
	thetaEl.resize( K );
	for ( PopIdx k = K ; k-- != 0 ; )
	    thetaEl[ k ] = one_over_K;
	}
    }



void Pedigree::SampleHiddenStates( const Options & /*options*/ )
    {

    #if NOT_NEEDED_FOR_FIXEDALLELEFREQ_EQ_1

	for ( ChromIdxType j = 0; j < numChromosomes; j++ )
	    {

	    const Chromosome & C = Loci->getChromosome( j );

	    // Update of forward probs here is unnecessary if SampleTheta was called
	    // and proposal was accepted.  Update Forward/Backward probs in HMM
	    if ( ! logLikelihood.HMMisOK )
		UpdateHMMInputs(j, options, Theta, _rho);

	    // Sampling locus ancestry can use current values of forward probability vectors alpha in HMM
	    C->HMM->SampleHiddenStates( LocusAncestry[j], (!isHaploid && (!Loci->isXChromosome(j) || SexIsFemale)) );

	    } // end chromosome loop

	logLikelihood.HMMisOK = true;

    #endif

    }


void Pedigree::SampleJumpIndicators( bool /*sampleArrivals*/ )
    {
    #if NOT_NEEDED_UNLESS_CONJUGATE_UPDATE
	for ( ChromIdxType chromIdx = 0 ; chromIdx < getNumChromosomes() ; chromIdx++ )
	    {
	    const Chromosome & chrom = getChromosome( chromIdx );

	    // Don't need to sample jump indicators if globalrho and no conjugate
	    // update of admixture this iteration sample number of arrivals, update
	    // SumNumArrivals and SumLocusAncestry

	    if ( ! chrom.isXChromosome() )
		chrom.SampleJumpIndicators( LocusAncestry[chromIdx], gametes[chromIdx],
					    SumLocusAncestry  , SumNumArrivals, sampleArrivals );
	    else
		chrom.SampleJumpIndicators( LocusAncestry[chromIdx], gametes[chromIdx],
					    SumLocusAncestry_X, SumNumArrivals, sampleArrivals );
	    } //end chromosome loop
    #endif
    }



//-----------------------------------------------------------------------------
// Overridden from PedBase
//-----------------------------------------------------------------------------

void Pedigree::SampleMissingOutcomes(bclib::DataMatrix *Outcome, const vector<bclib::Regression*>& R)
    {
    // This code directly copied from Individual::SampleMissingOutcomes() --
    // could be moved up to base class.

    int NumOutcomes = Outcome->nCols();
    // sample missing values of outcome variable
    for ( int k = 0; k < NumOutcomes; k++ )
	if ( Outcome->isMissing( getIndex(), k ) )
	    {
	    if ( R[k]->getRegressionType() == Linear )
		Outcome->set( getIndex(), k, RNG_NORMAL( R[k]->getExpectedOutcome(getIndex()), 1 / sqrt( R[k]->getlambda() ) ) );
	    else
		Outcome->set( getIndex(), k, (RNG_UNIFORM() * R[k]->getExpectedOutcome(getIndex()) < 1) ? 1 : 0 );
	    }
    }


bool Pedigree::GenotypeIsMissing( unsigned int clIdx ) const
    {
    #if TRACK_PEDIGREE_MISSING
	throw std::runtime_error( "Pedigrees don't yet track missing genotypes." );
    #else
	if ( clIdx ) {;} // Suppress compiler warning
	return false;
    #endif
    }


bool Pedigree::simpleGenotypeIsMissing( SLocIdxType slIdx ) const
    {
    #if TRACK_PEDIGREE_MISSING
	throw std::runtime_error( "Pedigrees don't yet track missing genotypes." );
    #else
	if ( slIdx ) {;} // Suppress compiler warning
	return false;
    #endif
    }


#if ! PEDIGREES_HAVE_PLOIDINESS

    bool Pedigree::isHaploidatLocus( unsigned int sLocIdx ) const
	{
	if ( getNMembers() == 1 )
	    return memberAt(0).getGType( sLocIdx ).isHaploid();
	else
	    throw std::runtime_error( "Pedigrees don't yet know what haploid is." );
	}

    bool Pedigree::isHaploidIndividual() const
	{
	if ( getNMembers() == 1 )
	    return memberAt(0).isHaploid();
	else
	    throw std::runtime_error( "Pedigrees don't yet know what haploid is." );
	}

#endif



//=============================================================================
// Methods from AdmixedIndividual (overridden from PedBase):
//=============================================================================

//-----------------------------------------------------------------------------
//
/// We don't need this for pedigrees because we calculate the
/// emission-probabilities (AKA genotype-probabilities) via a different
/// mechanism (Pedigree::genPossibleStates() et al, which gets called by
/// InputAdmixData::finishConstructing()).
///
/// @param chibIndicator Passed to CompositeLocus object.  If set to true,
///		CompositeLocus will use HapPairProbsMAP instead of HapPairProbs
///		when allelefreqs are not fixed.
//
//-----------------------------------------------------------------------------

void Pedigree::SetGenotypeProbs( /*ChromIdxType*/ int /*chrmIdx*/	    ,
				 int		      /*locusWithinChrom*/  ,
				 unsigned int	      /*absCompLocIdx*/	    ,
				 bool		      /*chibIndicator*/	    )
    {

    #if 1

	#if 0 // DEBUG
	    fprintf( stderr, "Warning: called unimplemented "
		"Pedigree::SetGenotypeProbs(this=%s,chrmIdx=%d,locusWithinChrom=%d,absCompLocIdx=%u,chibIndicator=%s)\n",
		getId().c_str(), chrmIdx, locusWithinChrom, absCompLocIdx, chibIndicator ? "true" : "false" );
	#endif

    #else


	#if THE_WAY_IT_WAS_BUT_CLEANED_UP_A_LITTLE

	    // === Code copied from AdmixedIndividual, may or may not be appropriate here ===

	    // DDF: Let's encapsulate isHaploid() in a method [of Pedigree?]
	    //	    Parameterized on what? simple-locus-index?
	    //	    (chromosome-index,composite-locus-within-chromosome-index)?

	    if( ! GenotypesMissing[chrmIdx][locusWithinChrom] )
		{
		if ( !isHaploid && ((chrmIdx != (int)X_posn) || SexIsFemale) )
		    { //diploid genotype
		    (*Loci)(locus)->GetGenotypeProbs(   GenotypeProbs[chrmIdx] + locusWithinChrom*NumHiddenStates*NumHiddenStates,
							PossibleHapPairs[locus],
							chibindicator);
		    }
		else
		    { //haploid genotype
		    (*Loci)(locus)->GetHaploidGenotypeProbs(
							GenotypeProbs[chrmIdx] + locusWithinChrom*NumHiddenStates,
							PossibleHapPairs[locus],
							chibindicator);
		    }
		}
	    else
		{
		if ( !isHaploid && (chrmIdx!=(int)X_posn || SexIsFemale))  //diploid genotype
		    for ( int k = 0; k < NumHiddenStates*NumHiddenStates; ++k )
			GenotypeProbs[chrmIdx][locusWithinChrom*NumHiddenStates*NumHiddenStates + k] = 1.0;
		else //haploid genotype
		    for ( int k = 0; k < NumHiddenStates; ++k )
			GenotypeProbs[chrmIdx][locusWithinChrom*NumHiddenStates + k] = 1.0;
		}

	#endif

    #endif

    }



void Pedigree::SampleHapPair( unsigned /*j*/, unsigned /*jj*/, unsigned /*locus*/, AlleleFreqs * /*A*/, bool /*skipMissingGenotypes*/, bool /*annealthermo*/, bool /*UpdateCounts*/ )
    {
    #if SUPPORT_ASSOCIATION_TESTS
	PedBase::SampleHapPair( j, jj, locus, A, skipMissingGenotypes, annealthermo, UpdateCounts );
    #endif
    }



#if PEDBASE_DEBUG_METHODS
    void Pedigree::dumpTheta( const char * prefix ) const
	{
	for ( FounderIdx fIdx = 0 ; fIdx < getNFounders() ; ++fIdx )
	    {
	    fprintf( stderr, "%s(%d) theta[%zd]: ", prefix, getMyNumber(), fIdx );
	    for ( size_t x = 0 ; x < thetas.getAccepted()[fIdx].size() ; ++x )
		fprintf( stderr, " %18.15lf", thetas.getAccepted()[fIdx][x] );
	    fprintf( stderr, "\n" );
	    if ( thetas.proposalActive() )
		{
		fprintf( stderr, "	theta-prop: " );
		for ( size_t x = 0 ; x < thetas.getProposed()[fIdx].size() ; ++x )
		    fprintf( stderr, " %18.15lf", thetas.getProposed()[fIdx][x] );
		fprintf( stderr, "\n" );
		}
	    }
	}
#endif



//-----------------------------------------------------------------------------
//
// drawInitialAdmixtureProps()
//
/// Draw initial values for admixture proportions theta from Dirichlet prior.
/// Operates on the "current" theta.
//
//-----------------------------------------------------------------------------

void Pedigree::drawInitialAdmixtureProps( const AlphaType & alpha )
    {
    const PopIdx K = getK();

    for ( FounderIdx fIdx = 0 ; fIdx < getNFounders() ; ++fIdx )
	{
	// Draw theta from Dirichlet with dirichlet parameters alpha:
	bclib::Rand::gendirichlet( K, alpha[0].getVector_unsafe(), getCurTheta()[fIdx] );
	}

    }



void Pedigree::ResetSufficientStats()
    {
    #if NEEDED_ONLY_FOR_CONJUGATE_UPDATE
	//************** Updating (Public) **********************************************************
	if(NumHiddenStates>1) {
	    // ** reset SumLocusAncestry to zero
	    for(int j = 0; j < NumHiddenStates *2; ++j) {
		SumLocusAncestry[j] = 0;
		SumLocusAncestry_X[j] = 0;
	    }

	    //SumNumArrivals is the number of arrivals between each pair of adjacent loci
	    fill(SumNumArrivals.begin(), SumNumArrivals.end(), 0);
	}
    #else
	#if 0
	    if ( ! options.getNoConjugateUpdate() )
		fprintf( stderr ,
			 "Warning: called unimplemented Pedigree::ResetSufficientStats(this=%s)\n" ,
			 getId().c_str() );
	#endif
    #endif
    }



//-----------------------------------------------------------------------------
// UpdateScores()
//-----------------------------------------------------------------------------

void Pedigree::UpdateScores( const AdmixOptions & options, bclib::DataMatrix *, bclib::DataMatrix *,
			const vector<bclib::Regression*> &, AffectedsOnlyTest & aoTest,
			CopyNumberAssocTest & )
    {
    if ( options.getTestForAffectedsOnly() )
	accumAOScore( aoTest );
    }



//-----------------------------------------------------------------------------
// getPosteriorMeans() [private helper for WritePosteriorMeans()
// DDF: code is adapted from code copied from AdmixedIndividual::getPosteriorMeans().
//-----------------------------------------------------------------------------

void Pedigree::getPosteriorMeans( ThetaType & thetaBar, RhoType & rhoMean, unsigned int samples ) const
    {

    const double dSamples = samples;

    rhoMean.resize( sumlogrho.size() );
    for ( size_t i = sumlogrho.size() ; i-- != 0 ; )
	rhoMean[i] = exp( sumlogrho[i] / dSamples );

    if ( thetaBar.size() != getNTheta() )
	thetaBar.resize( getNTheta() );


    // NOTE *1*: copy each element of ssmt rather than keep a reference, so the
    // scaled values can just be thrown away rather than scaled back when
    // finished.  N is not large.  This could perhaps be done more efficiently
    // on the stack.
    ThetaElType ssmtEl;


    // Apply softmax transformation to obtain thetabar
    for ( FounderIdx fIdx = 0 ; fIdx < getNTheta(); fIdx++ )
	{

	// Copy ssmt rather than keep a reference, so the scaled values can just
	// be thrown away rather than scaled back when finished.  N is not
	// large.  This could perhaps be done more efficiently on the stack.
	ssmtEl = SumSoftmaxTheta[ fIdx ];

	// Scale sumsoftmaxtheta
	ssmtEl /= dSamples;

	ThetaElType & thetaBarEl = thetaBar[fIdx];
	if ( thetaBarEl.size() != K )
	    thetaBarEl.resize( K );


	//---------------------------------------------------------------------
	// DDF: for some reason, we want to ("want to" = copied logic from
	// AdmixedIndividual) put the softmax transformation of SumSoftmaxTheta
	// into thetabar, but only for those elements for which the
	// corresponding element of the current value of theta is greater than 0.
	// The simplest way to do this is to set the values of ssmt at those
	// indexes to 0.0 (since they won't be used anyway), then use
	// softmax_gt0.
	//---------------------------------------------------------------------

	const ThetaElType & thetaEl = getCurTheta()[ fIdx ];

	gp_assert_eq( thetaEl.size() , thetaBarEl.size() );
	gp_assert_eq( thetaEl.size() , ssmtEl.size() );

	ThetaElType::iterator ssmtIt = ssmtEl.begin();
	for ( ThetaElType::const_iterator tIt = thetaEl.begin() ; tIt != thetaEl.end() ; ++tIt, ++ssmtIt )
	    if ( *tIt <= 0.0 )
		*ssmtIt = 0.0;


	//---------------------------------------------------------------------
	// DDF: It might make more sense if pvector::softmax_gt0() just resized
	// the output vector to the number of elements that qualify and only
	// pushed those elements onto it.  For now, we will roughly replicate
	// the behaviour of AdmixedIndividual::getPosteriorMeans(), which leaves
	// the values of non-qualifying elements in thetaBar undefined [raw
	// new'd contents, indeterminate per ISO C++ 1997-5.3.4(14A)]; instead,
	// we will keep slots for them but initialise to 0.
	//---------------------------------------------------------------------
	ssmtEl.softmax_gt0( thetaBarEl, 0.0 );

	}

    #if SUM_PROBS_IS_IMPLEMENTED
	if ( getOptions().getLocusAncestryProbsIndicator() )
	    SumProbs /= double(samples);
    #endif

    }



//-----------------------------------------------------------------------------
// WritePosteriorMeans()
// DDF: code is adapted from code copied from AdmixedIndividual::WritePosteriorMeans().
//-----------------------------------------------------------------------------

void Pedigree::WritePosteriorMeans( ostream& os, unsigned int samples, bool globalrho ) const
    {

    ThetaType thetaBar	;
    RhoType   rhobar	;

    getPosteriorMeans( thetaBar, rhobar, samples );

    for ( ThetaType::iterator it = thetaBar.begin() ; it != thetaBar.end() ; ++it )
	copy( it->begin(), it->end(), ostream_iterator<double>(os, "\t") );

    if ( ! globalrho )
	copy( rhobar.begin(), rhobar.end(), ostream_iterator<double>(os, "\t") );

    }



} // ---- end namespace genepi
