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


//=============================================================================
/// Issues yet to resolve:
///	* Rho indexing.  On what is it indexed for pedigrees?  Which element to
///	  pass into TPC constructor?  Update: for the moment, we are restrictng
///	  pedigrees to globalrho=1, so it's not relevant.
//=============================================================================


#include "Pedigree.h"


#include <cmath> // exp(), log()

#include <bclib/misc.h>
#include <bclib/rand.h>

#include "config.h"	// AGGRESSIVE_RANGE_CHECK
#include "TransProbCache.h"
#include "HiddenMarkovModel.new.h"


#if 0
    /// The way it was for individuals.
    #define IS_ADMIXED()  options.isAdmixed( this_gamete )
#else
    /// For pedigrees, we assume that all founders are admixed.
    #define IS_ADMIXED()  true
#endif


#define SEPARATE_HMM_FOR_EACH_CHROM 0


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



#define NOT_NEEDED_FOR_FIRST_VERSION			    0
#define THIS_IS_HOW_IT_WAS				    0
#define USE_SINGLE_POPULATION_SPECIAL_CASE		    0
#define NON_GLOBAL_RHO_WORKS				    0
#define NOT_NEEDED_FOR_FIXEDALLELEFREQ_EQ_1		    0
#define NOT_NEEDED_UNLESS_CONJUGATE_UPDATE		    0
#define SUPPORT_ASSOCIATION_TESTS			    0
#define IMPLEMENTED_AFFECTED_ONLY_SCORE_TEST_FOR_PEDIGREES  0
#define TRACK_PEDIGREE_MISSING				    0
#define PEDIGREES_HAVE_PLOIDINESS			    0
#define POSTERIOR_MEANS_IS_IMPLEMENTED			    1
#define SUM_PROBS_IS_IMPLEMENTED			    0


#if ! NOT_NEEDED_FOR_FIRST_VERSION
    #define IAmUnderTest    0
#endif



namespace genepi { // ----



/// This is bad; we need the concept of context!
const AdmixOptions * Pedigree::optionsRef = 0;



#if 0
    static bool greater_than_zero( double x )
	{
	return x > 0.0;
	}
#endif



//-----------------------------------------------------------------------------
//
// getGenome() [static]
//
/// Context access.  This is bad.
//
//-----------------------------------------------------------------------------

inline static Genome & getGenome()
    {
    return Individual::getGenome();
    }



//-----------------------------------------------------------------------------
// getK() [perhaps could be static but currently isn't]
// Nasty hack for context access.
//-----------------------------------------------------------------------------

PopIdx Pedigree::K = 0;

void Pedigree::setK( PopIdx _K )
    {
    if ( K != 0 )
	throw std::runtime_error( genepi::estr("Set K to ") + _K
				    + " after setting to " + K );

    K = _K;
    }



#if 0
    //-----------------------------------------------------------------------------
    // getNumChromosomes() [perhaps could be static but currently isn't]
    // Context access.
    //-----------------------------------------------------------------------------
    ChromIdxType Pedigree::getNumChromosomes() const
	{
	#error Not yet implemented.
	}
#endif



//-----------------------------------------------------------------------------
/// This stuff should be in the constructor, when we have a derived class.
//-----------------------------------------------------------------------------

void Pedigree::InitialiseAdmixedStuff( const AdmixOptions & options )
    {
    // This is bad.  Where is the concept of context?
    if ( optionsRef == 0 )
	optionsRef = &options;
    else
	gp_assert( optionsRef == &options );


    logLikelihood.value = 0.0;
    logLikelihood.ready = false;
    logLikelihood.HMMisOK = false;
    const double step0 = step;
    ThetaTuner.SetParameters( step0, 0.0001, 10.0, 0.44 );


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

    #if 0 && ! IMPLEMENTED_AFFECTED_ONLY_SCORE_TEST_FOR_PEDIGREES
	if ( options.getScoreTestIndicator() )
	    throw std::runtime_error( "Sorry, pedigrees are not (yet) compatible with affected-only score test."
					" Re-run with something else." );
    #endif

    #if ! SUPPORT_ASSOCIATION_TESTS
	if ( options.hasAnyAssociationTests() )
	    throw std::runtime_error( "Sorry, pedigrees are not (yet) compatible with association tests."
					" Re-run with some options turned off." );
    #endif



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


    _rho  .assign( NumGametes, 0.0 ); // set to 0 for unadmixed gametes
    rhohat.assign( NumGametes, 0.0 ); // set to 0 for unadmixed gametes


    for( unsigned g = 0; g < NumGametes; ++g)
	{
	if ( options.isAdmixed(g) )
	    {
	    _rho  [g] = init;
	    rhohat[g] = init;
	    }
	}

    sumlogrho.assign( _rho.size(), 0.0 );
    }



//-----------------------------------------------------------------------------
// getNInheritedByAffected()
//-----------------------------------------------------------------------------

int Pedigree::getNInheritedByAffected( FounderIdx fIdx, PopIdx k, const AncestryVector & av, const InheritanceVector & iv ) const
    {

    // Is this optimization worth having?  Does it come up in real datasets?
    #if 0
	if ( getNAffected() == 0 )
	    return 0; // *** RETURN HERE ***
    #endif

    int rv = 0;

    struct Ancestry
	{
	bool pAncestry;
	bool mAncestry;

	Ancestry( bool p, bool m ) : pAncestry(p), mAncestry(m) {}
	Ancestry() : pAncestry(false), mAncestry(false) {}

	bool hasAncOfType( const InheritanceVector::SegInd & si ) const
	    {
	    return (si == InheritanceVector::SI_PATERNAL) ? pAncestry : mAncestry;
	    }
	};

    Ancestry ancStack[ getNMembers() ]; // Could just be (getNMembers()-fIdx), then subtract indexes

    PopIdx mAncestry;
    const PopIdx pAncestry = av.getBoth( fIdx, mAncestry );
    ancStack[fIdx].pAncestry = (pAncestry == k);
    ancStack[fIdx].mAncestry = (mAncestry == k);

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
		    ++rv;
		}
	    Then recurse here...
	    }

    #else

	// Or do something like this...

	// FIXME: if the founder itself is affected, should that contribute to m?
	// If so, use this code:
	#if 0
	    if ( memberAt(fIdx).getOutcome() != 0 ) // isAffected()
		{
		if ( pAncestry == k )
		    ++rv;
		if ( mAncestry == k )
		    ++rv;
		}
	#endif

	for ( MemberIdx cIdx = fIdx + 1 ; cIdx < getNMembers() ; ++cIdx )
	    {

	    const Organism & child	= memberAt( cIdx );
	    Ancestry &	     cAncestry	= ancStack[ cIdx ];

	    const Organism & father	= *child.getFather();
	    const Ancestry & pAncestry	= ancStack[ father.getPIdx() ];

	    if ( pAncestry.hasAncOfType( iv.paternal(cIdx) ) )
		{
		cAncestry.pAncestry = true;
		if ( child.getOutcome() != 0 ) // isAffected()
		    ++rv;
		}

	    // Even if we already have one inherited k-state, we look for another in
	    // case of incest-cycle in a multi-generational pedigree
	    const Organism & mother = *child.getMother();
	    const Ancestry & mAncestry = ancStack[ mother.getPIdx() ];
	    if ( mAncestry.hasAncOfType( iv.maternal(cIdx) ) )
		{
		cAncestry.mAncestry = true;
		if ( child.getOutcome() != 0 ) // isAffected()
		    ++rv;
		}

	    }

    #endif

    return rv;

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
    else AGGRESSIVE_ONLY( if ( nFromK == 2 ) )
	return nAffOver2 * (1.0 - mu);
    AGGRESSIVE_ONLY( else throw std::logic_error( "invalid nFromK" ); )
    }


// Static helper:
static inline double aa_info( int nFromK, double nNPlus3over8, double nNPlus3over16, double n3over16, double mu )
    {
    if ( (nFromK == 0) || (nFromK == 2) )
	return nNPlus3over8 * (1 - mu) * mu;
    else AGGRESSIVE_ONLY( if ( nFromK == 1 ) )
	return (nNPlus3over16 * (1 - mu) * mu) - n3over16;
    AGGRESSIVE_ONLY( else throw std::logic_error( "invalid nFromK" ); )
    }


// Method:
void Pedigree::accumAOScore( AffectedsOnlyTest & aoTest ) const
    {

    const MemberIdx nAff = getNAffected();
    // Is this optimization worth having?  Does it come up in real datasets?
    if ( nAff == 0 )
	return; // *** RETURN HERE ***

    // Precompute invariants outside the loop:
    const double nAffOver2	= double(nAff) / 2;
    const double negNAffOver2	= - nAffOver2;
    const double nAffOver4	= double(nAff) / 4;
    const double nAffOver16	= double(nAff) / 16;
    const double nNPlus3	= nAff * (nAff + 3);
    const double n3over16	= (3.0 * nAff) / 16;
    const double nNPlus3over8	= nNPlus3 / 8;
    const double nNPlus3over16	= nNPlus3 / 16;

    const SimpleLocusArray & loci = getSLoci();

    // If K==2, only evaluate for 1 population, otherwise for all of them.
    // Existing code uses k==1 (not 0), yet treats it within the AffectedsOnlyTest object as "k==0", so:
    const PopIdx k_begin = (getK() == 2) ? 1 : 0;

//    for ( SLocIdxType t = loci.size() ; t-- != 0 ; )
    for ( SLocIdxType t = 0 ; t < loci.size() ; ++t )
	{

	// If K==2, only evaluate for 1 population, otherwise for all of them.
	// Existing code uses k==1 (not 0), so:
	//for ( PopIdx k = (getK() == 2) ? 1 : getK() ; k-- != 0 ; )
	for ( PopIdx k = k_begin ; k < getK() ; ++k )
	    {
	    double scoreAvg   = 0.0;
	    double scoreSqAvg = 0.0;
	    double infoAvg    = 0.0;

	    const HiddenStateSpace & hss	    = getStateProbs( t );
	    const cvector<double> &  condStateProbs = getHMM().getCondStateProbsAtLocus( t );

	    for ( HiddenStateSpace::Iterator stIt( hss ) ; stIt ; ++stIt )
		{

		double stScore = 0.0;
		double stInfo  = 0.0;

		for ( FounderIdx fIdx = getNFounders() ; fIdx-- != 0 ; )
		    {

		    const int nFromK = stIt.getAV().nCopiesFromKAtFounder( fIdx, k );

		    if ( stIt.getAV().isHetrozygousForPop( fIdx, k ) )
			{
			const int m = getNInheritedByAffected( fIdx, k, stIt.getAV(), stIt.getIV() );
			stScore += 0.25 * (m - nAffOver2);
			stInfo += nAffOver16;
			}

		    const double mu = Theta[ fIdx ][ k ];

		    stScore += aa_score( nFromK, nAffOver2, negNAffOver2, nAffOver4, mu );
		    stInfo  += aa_info( nFromK, nNPlus3over8, nNPlus3over16, n3over16, mu );

		    }

		const double stWeight = condStateProbs[ stIt.getNon0Index() ];
		const double wtScore = stScore * stWeight;
		scoreAvg += wtScore;
		scoreSqAvg += wtScore * stScore;
		infoAvg += stInfo * stWeight;
		}

	    const double condVarianceOfU = scoreSqAvg - (scoreAvg * scoreAvg);

	    aoTest.getAffectedsScore	( t, k-k_begin ) += scoreAvg	    ;
	    aoTest.getAffectedsVarScore ( t, k-k_begin ) += condVarianceOfU ;
	    aoTest.getAffectedsInfo	( t, k-k_begin ) += infoAvg	    ;

	    }

	}

    }



//-----------------------------------------------------------------------------
// setRho()
// Virtual, overridden from PedBase.  Only makes sense if each pedigree keeps its own HMM and TPC.
//-----------------------------------------------------------------------------

void Pedigree::setRho( double nv )
    {
    getTPC().setRho( nv );
    getHMM().thetaChanged();
    }



//-----------------------------------------------------------------------------
// HMMIsBad()
// What does this mean?  Exposition of internals to the outside world.
//-----------------------------------------------------------------------------

void Pedigree::HMMIsBad( bool /*loglikisbad*/ )
    {
    #if 0
	// Copied directly from Individual, could be placed in some base class,
	// except seems to have no bearing on pedigrees since they track for
	// themselves whether the HMM is OK.
	logLikelihood.HMMisOK = false;
	if ( loglikisbad )
	    logLikelihood.ready = false;
    #endif
    }



//-----------------------------------------------------------------------------
// getHMM()
//-----------------------------------------------------------------------------

TransProbCache & Pedigree::getTPC() const
	{
	#if NON_GLOBAL_RHO_WORKS
	    #error Which rho to use here?
	#else
	    // Assume that _rho[0] is the global rho to use:
	    if ( tpCache == 0 )
		tpCache = new TransProbCache( *this, _rho[0], Theta );
	#endif
	return *tpCache;
	}

    /// Allocate the HiddenMarkovModel for the Pedigree object (only if have one HMM
    /// for each pedigree rather than one per chromosome).

    HiddenMarkovModel & Pedigree::getHMM() const
	{
	if ( hmm == 0 )
	    hmm = new HiddenMarkovModel( *this, getTPC(), &Theta );
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
		 double *			    SumLogTheta		,
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
	    logpratio += ProposeThetaWithRandomWalk( options, alpha );
	    }
	  #if NEEDED_ONLY_FOR_CONJUGATE_UPDATE
	    else
		ProposeTheta(options, alpha, SumLocusAncestry, SumLocusAncestry_X);
	  #endif
	}
    catch ( string & s )
	{
	throw std::runtime_error( "Error encountered while generating proposal individual admixture proportions:\n  " + s );
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
    Accept_Reject_Theta(logpratio, K, options.isRandomMatingModel(), RW );

    // Update the value of admixture proportions used in the regression model, but
    // only if pedigree is a single individual:
    if ( (options.getNumberOfOutcomes() > 0) && (getNMembers() == 1) )
	UpdateAdmixtureForRegression(K, NumCovariates, poptheta, options.isRandomMatingModel(), Covariates);

    if ( (! anneal) && (iteration > options.getBurnIn()) )
	{
	// Accumulate sums in softmax basis for calculation of posterior means

	double a[ K ]; // Or could use: bclib::pvector<double> a;

	for ( ThetaIdx tIdx = 0 ; tIdx < getNTheta() ; ++tIdx )
	    {
	    Theta[tIdx].inv_softmax_gt0( a ); // Equivalent: Theta[tIdx].inv_softmax( a, gt_0 );
	    SumSoftmaxTheta[tIdx] += a;
	    }
	}


    if ( ! IAmUnderTest )
	{
	for ( ThetaIdx tIdx = 0 ; tIdx < getNTheta() ; ++tIdx )
	    {
	    //SumLogTheta.addAfterFunction( Theta[tIdx], log );
	    const pvector<double> & th = Theta[ tIdx ];
	    for ( size_t k = 0 ; k < th.size() ; ++k )
		SumLogTheta[k] += log( th[k] );
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
	    avgtheta[k] = (ThetaProposal[k] + ThetaProposal[k + NumHiddenStates ])/ 2.0 - poptheta[k];
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
// UpdateAdmixtureForRegression()
// Update individual admixture values (mean of both gametes) used in the regression model
//-----------------------------------------------------------------------------

void Pedigree::UpdateAdmixtureForRegression( int NumHiddenStates, int NumCovariates,
                                               const PopThetaType & poptheta, bool /*RandomMatingModel*/,
					       bclib::DataMatrix * Covariates )
    {
    cvector<double> avgtheta( NumHiddenStates );


    for ( PopIdx k = getK() ; k-- != 0 ; )
	{
	double sum = 0.0;
	for ( ThetaIdx tIdx = getNTheta() ; tIdx-- != 0 ; )
	    sum += Theta[ tIdx ][ k ];
	avgtheta[ k ] = sum / getNTheta();
	}


    const PopIdx K = getK();
    const size_t NumberOfInputCovariates = NumCovariates - K;
    // Seems strange, why does this loop start at 1?
    for( PopIdx k = 1 , covIdx = NumberOfInputCovariates ; k < K ; k++, covIdx++ )
	Covariates->set( getIndex(), covIdx, avgtheta[ k ] - poptheta[ k ] );
    }



//-----------------------------------------------------------------------------
// ProposeThetaWithRandomWalk()
//-----------------------------------------------------------------------------

double Pedigree::ProposeThetaWithRandomWalk( const AdmixOptions & options, const AlphaType & alpha )
    {
    double LogLikelihoodRatio = 0.0;
    double LogPriorRatio      = 0.0;


//fprintf( stderr, "RNG check 2: %20.18lf\n", bclib::Rand::gennor( 20, 10 ) );

    //-------------------------------------------------------------------------
    // Generate proposals
    //-------------------------------------------------------------------------

    const PopIdx K = getK();

    // NOTE *M1*: Maybe this should be be allocated once for the class?
    pvector<double> a( K );

    for ( ThetaIdx tIdx = 0 ; tIdx < getNTheta() ; ++tIdx )
	{

	if ( IS_ADMIXED() )
	    {
	    const pvector<double> & th	    = Theta	   [ tIdx ];
	    pvector<double> &	    th_prop = ThetaProposal[ tIdx ];

	    th.inv_softmax_gt0( a );

	    // Random walk step - on all elements of array a
	    for ( PopIdx k = 0 ; k < K ; ++k )
		if ( a[k] != 0.0 )
		    a[k] = bclib::Rand::gennor( a[k], step );

	    // Reverse transformation from numbers on real line to proportions
	    a.softmax_gt0( th_prop );


	    // Compute contribution of this founder to log prior ratio
	    // Prior densities must be evaluated in softmax basis
	    for ( PopIdx k = 0; k < K; ++k )
		LogPriorRatio += alpha[tIdx][k] * (log(th_prop[k]) - log(th[k]));
	    }

	else // IS_ADMIXED()
	    Theta[ tIdx ] = ThetaProposal[ tIdx ];

	} // End loop over Theta's elements (i.e. founders)

    // Get log likelihood at current parameter values - do not force update, store result of update
    LogLikelihoodRatio -= getLogLikelihood( options, false, true );

    // Get log likelihood at proposal theta and current rho - force update
    // store result in loglikelihood.tempvalue, and accumulate loglikelihood ratio
    logLikelihood.tempvalue = getLogLikelihood( options, ThetaProposal, _rho, true );
    LogLikelihoodRatio += logLikelihood.tempvalue;

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

void Pedigree::Accept_Reject_Theta( double logpratio, int /*NumHiddenStates*/, bool /*RandomMatingModel*/, bool RW )
    {
    // Loop over populations: if any element of proposed Dirichlet parameter
    // vector is too small, reject update without test step
    bool should_do_test = true;
    for ( ThetaIdx tIdx = getNTheta() ; should_do_test && (tIdx-- != 0) ; )
	{
	const pvector<double> & th	= Theta	       [ tIdx ];
	const pvector<double> & th_prop = ThetaProposal[ tIdx ];

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
	    if ( bclib::Rand::myrand() < AccProb )
		accept = true;
	    }
	else
	    accept = true;
	}


    if ( accept )
	{
	// Set proposed values as new values
	setAdmixtureProps( ThetaProposal );
	if ( RW )
	    {
	    // If random-walk update, store the temp log-likelihood and set loglikelihood.HMMisOK to true
	    storeLogLikelihood( true );
	    }
	else
	    {
	    // conjugate update of parameters invalidates both HMM forward probs and stored loglikelihood
	    logLikelihood.HMMisOK = false;
	    logLikelihood.ready = false;
	    }
	}
    else
	{
	// If RW proposal is rejected, loglikelihood.HMMisOK is already set to
	// false, and stored log-likelihood is still valid
	}


    if ( RW )
	{
	// Update step size in tuner object every w updates
	if ( (NumberOfUpdates % w) == 0 )
	    step = ThetaTuner.UpdateStepSize( AccProb );
	}

    }



//-----------------------------------------------------------------------------
// storeLogLikelihood()
//
/// To call if a Metropolis proposal is accepted
//-----------------------------------------------------------------------------

void Pedigree::storeLogLikelihood( bool setHMMAsOK )
    {
    logLikelihood.value = logLikelihood.tempvalue;
    logLikelihood.ready = true;
    if ( setHMMAsOK )
	logLikelihood.HMMisOK = true;
    }



//-----------------------------------------------------------------------------
/// Set admixture proportions
//-----------------------------------------------------------------------------

void Pedigree::setAdmixtureProps( const ThetaType & rhs )
    {
    gp_assert( Theta.size() == getNTheta() );
    gp_assert( rhs  .size() == getNTheta() );

    for ( ThetaIdx tIdx = getNTheta() ; tIdx-- != 0 ; )	
{
	pvector<double> &	sub_theta = Theta [ tIdx ];
	const pvector<double> & sub_rhs	  = rhs	  [ tIdx ];

	gp_assert( sub_theta.size() == getK() );
	gp_assert( sub_rhs  .size() == getK() );

	for ( PopIdx k = getK() ; k-- != 0 ; )
	    sub_theta[ k ] = sub_rhs[ k ];
	}
    }



//*****************************************************************************
//***************************  Log-Likelihoods  *******************************
//*****************************************************************************

//-----------------------------------------------------------------------------
// public function:
// Calls private function to get log-likelihood at current parameter values, and
// stores it either as loglikelihood.value or as loglikelihood.tempvalue store
// should be false when calculating energy for an annealed run, or when
// evaluating proposal for global sum-intensities
// **** COPIED FROM AdmixedIndividual, NOT Individual ****
//-----------------------------------------------------------------------------

double Pedigree::getLogLikelihood( const Options & options, bool forceUpdate, bool store )
    {
    if ( (! logLikelihood.ready) || forceUpdate )
	{
	logLikelihood.tempvalue = getLogLikelihood( options, Theta, _rho, true );
	if ( store )
	    {
	    logLikelihood.value = logLikelihood.tempvalue;
	    logLikelihood.ready = false;    //true;
	    logLikelihood.HMMisOK = false;  //true; //because forward probs now correspond to current parameter values
	    }					    //and call to UpdateHMM has set this to false
	return logLikelihood.tempvalue;
	}
    else
	return logLikelihood.value; // nothing was changed
    }



//-----------------------------------------------------------------------------
//
// Private function: gets log-likelihood at parameter values specified as
// arguments, but does not update loglikelihoodstruct
//
// This is a combination of the AdmixedIndividual:: and Individual:: versions,
// because the AdmixedIndividual version called the superclass's version
// internally.
//
//-----------------------------------------------------------------------------

double Pedigree::getLogLikelihood(  const Options &	options	  ,
				    const ThetaType &	theta	  ,
				    const RhoType &	rho	  ,
				    bool		updateHMM ) const
    {
    double rv;

    #if USE_SINGLE_POPULATION_SPECIAL_CASE
      if ( getK() == 1 )
	rv = getLogLikelihoodOnePop();
      else
    #endif
	{
	// ==== This was the body of the Individual:: version, a superclass-call
	// ==== from AdmixedIndividual::getLogLikelihood()
	// ==== rv = Individual::getLogLikelihood(options, theta, rho, updateHMM);


	#if SEPARATE_HMM_FOR_EACH_CHROM
	    rv = 0.0;
	    for ( ChromIdxType chromIdx = 0 ; chromIdx < getNumChromosomes() ; chromIdx++ )
		{
		if ( updateHMM )
		    {// force update of forward probs
		    UpdateHMMInputs( chromIdx, options, theta, rho );
		    }

		// For individuals: is_diploid = !isHaploid && (!Loci->isXChromosome(j) || SexIsFemale)
		// Should be: getChromosome(j).isDiploid()
		// But in fact, since we have only one HMM per chromosome, it should
		// remember its state-space, so not need to have is-diploid passed
		// into each method-call.
		rv += getHMM( chromIdx ).getLogLikelihood();
		}
	#else
	    if ( updateHMM )
		updateHMMInputs( options, theta, rho );
	    rv = getHMM().getLogLikelihood();
	#endif

	// FIXME: under what circumstances does this happen?  Throw the proper
	// exception from within HMM instead.
	if ( isnan(rv) )
	    throw std::runtime_error( "HMM returns log-likelihood as nan (not a number)" );
	}

    return rv;
    }



//-----------------------------------------------------------------------------
// getLogLikelihoodAtPosteriorMeans()
//-----------------------------------------------------------------------------

double Pedigree::getLogLikelihoodAtPosteriorMeans( const Options & options )
    {

#if POSTERIOR_MEANS_IS_IMPLEMENTED

    double rv;

    // Should set allele freqs also to posterior means, and recalculate prob
    // genotypes at these freqs before calling getloglikelihood

    #if USE_SINGLE_POPULATION_SPECIAL_CASE
      if ( getK() == 1 )
	rv = getLogLikelihoodOnePop();
      else
    #endif
	{
	ThetaType thetabar( getNTheta() );
	RhoType	  rhobar;

	rhobar.reserve( sumlogrho.size() );
	for ( unsigned i = 0 ; i < sumlogrho.size() ; ++i )
	    rhobar.push_back( exp(sumlogrho[i]/ double(options.getTotalSamples() - options.getBurnIn())) );

	const double scale_factor = options.getTotalSamples() - options.getBurnIn();

	SumSoftmaxTheta /= scale_factor;

	// apply softmax transformation to obtain thetabar
	for ( ThetaIdx tIdx = 0 ; tIdx < getNTheta() ; ++tIdx )
	    SumSoftmaxTheta[tIdx].softmax_gt0( thetabar[tIdx] ); // Equivalent: SumSoftmaxTheta[tIdx].softmax( thetabar[tIdx], gt_0 );

	// Rescale sumsoftmaxtheta back
	SumSoftmaxTheta *= scale_factor;

	updateHMMInputs( options, thetabar, rhobar );
	rv = getHMM().getLogLikelihood();
	}

#else // POSTERIOR_MEANS_IS_IMPLEMENTED

    if ( &options == 0 ) {;} // Suppress compiler warning
    const double rv = 0.0;

#endif // ! POSTERIOR_MEANS_IS_IMPLEMENTED

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
// UpdateHMMInputs()
//-----------------------------------------------------------------------------

#if SEPARATE_HMM_FOR_EACH_CHROM

    void Pedigree::UpdateHMMInputs(
	    ChromIdxType	chromIdx,
	    const Options &	options ,
	    const ThetaType &	theta	,
	    const RhoType &	rho	) const
	{
	// Updates inputs to HMM for chromosome j also sets Diploid flag in
	// Chromosome (last arg of SetStateArrivalProbs)
	const bool diploid = !isHaploid && (j!=X_posn || SexIsFemale);
	const bool isRandomMating = options.isRandomMatingModel();

	const Chromosome & chrom = getChromosome( chromIdx );

	C->HMM->SetGenotypeProbs( GenotypeProbs[j], GenotypesMissing[j] );
	//   if(!SexIsFemale) cout << "chr " << j << " " << X_posn << " ";
	//   if(!diploid) cout << "haploid GenotypeProbs set for chromosome " << j << "\n";

	if ( ! options.isGlobalRho() )
	    {
	    // Set locus correlation, f, if individual- or gamete-specific rho
	    chrom.SetLocusCorrelation( rho, isRandomMating );
	    }

	// Set the state-arrival-probabilities in the HMM:
	// For individuals we do this:
	//	    "in the haploid case in random mating model, pass pointer to maternal admixture props."
	// For pedigrees, we will pass the entire Theta array into the HMM, regardless.
	#if THIS_IS_HOW_IT_WAS
	    if ( diploid || !isRandomMating )
		C->HMM->SetStateArrivalProbs(theta, (int)isRandomMating, diploid);
	    else
		C->HMM->SetStateArrivalProbs(theta + NumHiddenStates, true, false);
	#else
	    C->HMM->SetStateArrivalProbs(theta + NumHiddenStates, true, false);
	#endif

	logLikelihood.HMMisOK = false;//because forward probs in HMM have been changed
	}

#else

    /// Updates inputs to the HMM (which might seem odd considering the name).
    void Pedigree::updateHMMInputs( const Options & options, const ThetaType & theta, const RhoType & /*rho*/ ) const
	{

	// For some reason, rho is not updated here if global-rho is turned on
	// (it is updated in the Chromosome/TPC from PopAdmix).
	if ( ! options.isGlobalRho() )
	    {
	    #if NON_GLOBAL_RHO_WORKS
		getTPC().setRho( something else );
		getHMM().thetaChanged();
	    #else
		throw std::runtime_error( "non-global-rho not yet supported." );
	    #endif
	    }

	getHMM().setTheta( &theta );

	logLikelihood.HMMisOK = false; // Because forward probs in HMM have been changed
	}

#endif



//-----------------------------------------------------------------------------
// SetUniformAdmixtureProps() [protected]
// Set the initial values of theta.
//-----------------------------------------------------------------------------

void Pedigree::SetUniformAdmixtureProps()
    {
    bclib::pvector<double> uniform;
    uniform.resize( getK(), 1.0 / K );

    for ( size_t fIdx = getNTheta() ; fIdx-- != 0 ; )
	Theta[ fIdx ] = uniform;
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
		Outcome->set( getIndex(), k, bclib::Rand::gennor( R[k]->getExpectedOutcome(getIndex()), 1 / sqrt( R[k]->getlambda() ) ) );
	    else
		Outcome->set( getIndex(), k, (bclib::Rand::myrand() * R[k]->getExpectedOutcome(getIndex()) < 1) ? 1 : 0 );
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


bool Pedigree::isHaploidatLocus( unsigned int ) const
    {
    #if PEDIGREES_HAVE_PLOIDINESS
	throw std::runtime_error( "Pedigrees don't yet know what haploid is." );
    #else
	return false;
    #endif
    }


bool Pedigree::isHaploidIndividual() const
    {
    #if PEDIGREES_HAVE_PLOIDINESS
	throw std::runtime_error( "Pedigrees don't yet know what haploid is." );
    #else
	return false;
    #endif
    }



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
	    fprintf( stderr, "%s theta (%d,%zd): ", prefix, getMyNumber(), fIdx );
	    for ( size_t x = 0 ; x < Theta[fIdx].size() ; ++x )
		fprintf( stderr, " %18.15lf", Theta[fIdx][x] );
	    fprintf( stderr, "\n" );
	    }
	}
#endif


// draw initial values for admixture proportions theta from Dirichlet prior
void Pedigree::drawInitialAdmixtureProps( const AlphaType & alpha )
    {
      const PopIdx K = getK();

      for ( FounderIdx fIdx = 0 ; fIdx < getNFounders() ; ++fIdx )
	{

	double sum = 0.0;

	for ( PopIdx k = 0 ; k < K ; ++k )
	    sum += alpha[fIdx][k];

	for ( PopIdx k = 0 ; k < K ; ++k )
	    {
	    thetahat[fIdx][k] = alpha[fIdx][k] / sum; // set thetahat to prior mean
	    dirparams[k] = alpha[fIdx][k];
	    }

	// draw theta from Dirichlet with parameters dirparams
	bclib::Rand::gendirichlet( K, dirparams.getVector_unsafe(), Theta[fIdx] );

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
//
// UpdateScores()
//
/// DDF: Code is direct copy from AdmixedIndividual, but cleaned up a little.
///
/// Copied comment: "Merge with updatescoretests"
//
//-----------------------------------------------------------------------------

#if IMPLEMENTED_AFFECTED_ONLY_SCORE_TEST_FOR_PEDIGREES
    void Pedigree::UpdateScores( const AdmixOptions & options, bclib::DataMatrix * Outcome, bclib::DataMatrix * Covariates,
				const vector<bclib::Regression*> & R, AffectedsOnlyTest & affectedsOnlyTest,
				CopyNumberAssocTest& ancestryAssocTest )
	{
	// DDF: admixtureCovars is only needed if options.getTestForLinkageWithAncestry() is
	//	    turned on, but since it apparently is not too large, just allocate
	//	    it on the stack here.
	double admixtureCovars[ NumHiddenStates - 1 ];

	for ( ChromIdxType j = 0; j < numChromosomes; j++ )
	    {
	    Chromosome & C = getChromosome(j);

	    // Update of forward probs here is unnecessary if SampleTheta was called and proposal was accepted
	    // Update Forward/Backward probs in HMM
	    if ( ! logLikelihood.HMMisOK )
		UpdateHMMInputs(j, options, Theta, _rho);

	    // Update of score tests for linkage with ancestry requires update of backward probs
	    if ( options.getTestForLinkageWithAncestry() )
		for ( int t = NumHiddenStates ; t-- != 0 ; )
		    admixtureCovars[t] = Covariates->get( getIndex(), Covariates->nCols()-NumHiddenStates+1+t );

	    UpdateScoreTests( options, options.getTestForLinkageWithAncestry() ? admixtureCovars : 0,
				Outcome, C, R, affectedsOnlyTest, ancestryAssocTest);
	    } //end chromosome loop
	}
#else

    void Pedigree::UpdateScores( const AdmixOptions & options, bclib::DataMatrix *, bclib::DataMatrix *,
				const vector<bclib::Regression*> &, AffectedsOnlyTest & aoTest,
				CopyNumberAssocTest & )
	    {
	    if ( options.getTestForAffectedsOnly() )
		accumAOScore( aoTest );
	    }

#endif



#if POSTERIOR_MEANS_IS_IMPLEMENTED

//-----------------------------------------------------------------------------
// getPosteriorMeans() [private helper for WritePosteriorMeans()
// DDF: code is adapted from code copied from AdmixedIndividual::getPosteriorMeans().
//-----------------------------------------------------------------------------

void Pedigree::getPosteriorMeans( ThetaType & thetaBar, RhoType & rhoMean, unsigned int samples ) const
    {
    const double dSamples = samples;

    #if 0
	for ( size_t i = SumSoftmaxTheta.size() ; i-- != 0 ; )
	    SumSoftmaxTheta[i] /= dSamples;
    #else
	for ( ThetaType::iterator it = SumSoftmaxTheta.begin() ; it != SumSoftmaxTheta.end() ; ++it )
	    *it /= dSamples;
    #endif

    rhoMean.clear();
    for ( size_t i = sumlogrho.size() ; i-- != 0 ; )
	rhoMean.push_back( exp( sumlogrho[i] / dSamples ) );

    thetaBar.reserve( getNTheta() );

    // Apply softmax transformation to obtain thetabar
    for ( FounderIdx fIdx = 0 ; fIdx < getNTheta(); fIdx++ )
	{
	const bclib::pvector<double> & thetaEl = Theta[fIdx];
	if ( thetaBar.size() <= fIdx )
	    thetaBar.resize( fIdx + 1 );
	bclib::pvector<double> & thetaBarEl = thetaBar[fIdx];

	#if 0
	    // Keep track of qualifying elements in boolean array on the stack.
	    // I believe that this should not be necessary.
	    bool qualifies[ thetaEl.size() ];
	    for ( PopIdxType k = getK() ; k-- != 0 ; )
		qualifies[k] = (thetaEl[k] > 0.0);
	    thetaEl.softmax_gt0( thetaBarEl );
	#else
	    // DDF: It might make more sense if pvector::softmax_gt0() just
	    // resized the output vector to the number of elements that qualify
	    // and only pushed those elements onto it.  For now, we will roughly
	    // replicate the behaviour of
	    // AdmixedIndividual::getPosteriorMeans(), which leaves the values
	    // of non-qualifying elements in thetaBar undefined [raw new'd
	    // contents, indeterminate per ISO C++ 1997-5.3.4(14A)]; instead, we
	    // will keep slots for them but initialise to 0.
	    thetaEl.softmax_gt0( thetaBarEl );
	#endif

	// Rescale sumsoftmaxtheta back (DDF: why weren't the original values preserved?)
	SumSoftmaxTheta[ fIdx ] *= dSamples;
	}

    #if SUM_PROBS_IS_IMPLEMENTED
	if ( getOptions().getLocusAncestryProbsIndicator() )
	    SumProbs /= double(samples);
    #endif

    }

#endif // ! POSTERIOR_MEANS_IS_IMPLEMENTED



//-----------------------------------------------------------------------------
// WritePosteriorMeans()
// DDF: code is adapted from code copied from AdmixedIndividual::WritePosteriorMeans().
//-----------------------------------------------------------------------------

void Pedigree::WritePosteriorMeans( ostream& os, unsigned int samples, bool globalrho ) const
    {

#if POSTERIOR_MEANS_IS_IMPLEMENTED

    ThetaType thetaBar	;
    RhoType   rhobar	;

    getPosteriorMeans( thetaBar, rhobar, samples );

    for ( ThetaType::iterator it = thetaBar.begin() ; it != thetaBar.end() ; ++it )
	{
	copy( it->begin(), it->end(), ostream_iterator<double>(os, "\t") );
	os << '\n';
	}

    if ( ! globalrho )
	copy( rhobar.begin(), rhobar.end(), ostream_iterator<double>(os, "\t") );

#else // POSTERIOR_MEANS_IS_IMPLEMENTED

    //fprintf( stderr, "Warning: posterior-means is not implemented.\n" );
    if ( os	   ) {;} // Suppress compiler warning
    if ( samples   ) {;} // Suppress compiler warning
    if ( globalrho ) {;} // Suppress compiler warning

#endif // ! POSTERIOR_MEANS_IS_IMPLEMENTED

    }



} // ---- end namespace genepi
