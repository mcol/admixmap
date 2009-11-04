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


#include "config.h" // AGGRESSIVE_RANGE_CHECK


#include <cstddef>  // size_t

#include <bclib/cvector.h>
#include "Genotype.h"
#include "InheritanceVector.h"
#include "OrganismArray.h"
#include "PedBase.h"
#include "ProposalRingBuffer.h"
#include "SimpleLocusArray.h"

#include "AlleleArray.h" // AlleleProbTable, AlleleProbVect (for emission probabilities)

#include <vector>


// Need to sort this out:
class Chromosome;
class CompositeLocus;
namespace genepi
    {
    typedef size_t CLocIdxType;
    typedef size_t ChromIdxType;
    };



// For admix-model:
//------------------------------------------------------------------------
// See NOTE *4*:
// Used for the additional methods to run the admixmap model.  These should be
// moved into a derived class in admixmap, not here in libbase.
class AdmixOptions;
class Options;
class CopyNumberAssocTest;
class AdmixOptions;
class AffectedsOnlyTest;
#include "TwoDimArray.h" // aoCache
namespace bclib { class DataMatrix; }
#include "common.h" // for DataType
#include <bclib/StepSizeTuner.h>
#if 0 // Circular dependencies
    #include "HiddenMarkovModel.new.h"
    #include "TransProbCache.h"
#else
    namespace genepi { class HiddenMarkovModel; }
    namespace genepi { class TransProbCache; }
#endif

/// See getRNG().
#define PED_HAS_OWN_PRNG defined(_OPENMP)
#if PED_HAS_OWN_PRNG
    #include <bclib/NRand.h>
    #define RNG_NORMAL(M,S) getRNG().normal(M,S)
    #define RNG_UNIFORM()   getRNG().standard_uniform()
#else
    #define RNG_NORMAL(M,S) Rand::gennor(M,S)
    #define RNG_UNIFORM()   Rand::myrand()
#endif
//------------------------------------------------------------------------


namespace genepi { // ----

/** \addtogroup base
 * @{ */



// Needed for emission probability computation:
class AncestryVector;
class HiddenStateSpace;



//-----------------------------------------------------------------------------
//
//  PEDIGREE CLASS
//
/// A class to represent a pedigree as read in from the pedigree file; the list
/// of the @ref Organism "Organisms" in the pedigree, as well as aggregate and
/// probability information about the pedigree.
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
///	<TD STYLE="padding-top: 6pt; padding-bottom: 0pt;">
///	We implement an ordering of the Organism within the Pedigree for
///	iteration such that no member of a family will be visited prior to both
///	of its parents being visited.  This ordering is achieved by assigning
///	each non-founder (or even founders too) a "depth" which will be the
///	<I>largest</I> number of steps required to traverse the parent-tree
///	starting at that node.  It is similar to "generation", except that all
///	founders of a given pedigree have the same depth of 0 (despite the fact
///	that some may be of a different generation than others), and a given
///	individual's depth is then the <I>maximum</I> path-length to a founder.
///	Another way to say this is that every sib's depth is the maximum of its
///	two parents' depths, plus 1.  Thus, all founders have depth 0; a child
///	both of whose parents are founders has depth 1; and a child one of whose
///	parents is a founder and the other of whom is the child of two founders
///	has depth 2.
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
///	The Pedigree does not "own" the @link Organism Organisms @endlink within
///	it (i.e. will not delete them when it is deleted).  They remain
///	controlled by the OrganismArray, a reference to which the Pedigree
///	keeps.
///	</TD>
///  </TR>
///
/// </TABLE>
///
/// <A name="note-4"></A>
/// <TABLE STYLE="border: groove 3pt aqua;">
///
///  <TR>
///	<TD><B>NOTE *4*</B></TD>
///	<TD>
///	Additional data members and methods to run the admixmap model. These
///	should not be here in the base class, but rather should be moved into a
///	derived class in admixmap, not here in libbase (the implementation of
///	those methods is already kept in a separate source-code file in
///	admixmap, AdmixPedigree.cc).  This also applies to deriving from
///	PedBase, which should take place as multiple inheritance in
///	AdmixPedigree, like this:
///
///	<CODE>
///	class AdmixPedigree : public Pedigree , public PedBase
///	    { ... };
///	</CODE>
///	</TD>
///  </TR>
///
/// </TABLE>
//
//-----------------------------------------------------------------------------

class Pedigree : public PedBase // See NOTE *4*
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
					SLocIdxType		  sLocIdx	  ,
					const AncestryVector &	  av		  ,
					const InheritanceVector & iv		  ,
					const Haplotype *	  founderHapState ,
					double			  emProb	  );


    private:
	const OrganismArray & memberPool;

	IdType	   id		 ;
	MemberIdx  nMembers	 ;
	FounderIdx nFounders	 ;
	Member * * sortedMembers ;

	InheritanceSpace * ivSpace;

	// Mendelian error detection and tracking:
	mutable int    nMendelErrs;
	mutable bool * mendelErrsByLocus; ///< Array, indexed on SLocIdxType
	bool & mendelErrAt( SLocIdxType t ) const;


	/// Array of hidden-state-spaces, one for each locus.  This belongs in a
	/// subclass.  We use a pointer-to-array (indexed on locus-index) rather
	/// than std::vector to avoid including the full class definition here.
	mutable HiddenStateSpace * stateProbs;


	void recurseSib( SLocIdxType		sLocIdx		,
			 Haplotype *		memberHapState	,
			 const AncestryVector & ancestry	,
			 InheritanceVector &	iv		,
			 MemberIdx		memDepth	,
			 double			emProbTerm	,
			 StateReceiver		receiver	) const;

	void recurseFounder(	SLocIdxType		sLocIdx		,
				PopIdx			K		,
				const AlleleProbTable & alProbTab	,
				Haplotype *		memberHapState	,
				AncestryVector &	ancestry	,
				MemberIdx		memDepth	,
				double			probProdSoFar	,
				StateReceiver		receiver	) const;


    protected:

	/// The constructor is protected; use generatePedigrees() to create the
	/// Pedigree objects.
	/// @param max_n If non-0, pedigrees with sizes exceeding max_n with be
	///	reduced to size max_n by eliminating non-founders; if there are
	///	more than max_n founders, an exception will be thrown.
	Pedigree( const OrganismArray &	      pool	,
		  OrganismArray::ConstPedIter firstM	,
		  OrganismArray::ConstPedIter endM	,
		  MemberIdx		      max_n = 0 );


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
	///
	/// Create a vector of Pedigree objects from the raw genotype data (as
	/// read in from a pedfile by GenotypeParser).
	///
	/// This is effectively the public "constructor" for Pedigree objects.
	/// The Pedigree objects are created in @a rv (output parameter).
	///
	/// @param max_n See Pedigree::Pedigree()
	///
	//---------------------------------------------------------------

	static void generatePedigrees( const OrganismArray & organisms	,
				       cvector<Pedigree> &   rv		,
				       MemberIdx	     max_n	);


	//---------------------------------------------------------------
	// Access to our "context","domain" (locus file, populations, etc.).
	// Perhaps these could be static.
	// Many are convenience methods (simply derived from others).
	//---------------------------------------------------------------

	private: static PopIdx K; public:
	PopIdx	     getK	      () const	///< The number of populations
						{ gp_assert_ne( K, 0 ); return K; }
	ChromIdxType getNumChromosomes() const; ///< The number of chromosomes
	static void setK( PopIdx _K );

	const Chromosome & getChromosome( ChromIdxType chromIdx ) const;

	/// Return the pool (pedfile) from which members are drawn
	const OrganismArray &	 getMemberPool() const { return memberPool; }
	const SimpleLocusArray & getSLoci     () const { return memberPool.getSLoci(); } ///< Convenience
	SLocIdxType		 getNSLoci    () const { return getSLoci().size()    ; } ///< Convenience

	#if 0 // Unimplemented; composite-loci do not as yet work with pedigrees
	    const CompositeLocus & getCLocus( CLocIdxType ) const;
	    CLocIdxType		   getNCLoci() const;
	#endif


	//---------------------------------------------------------------
	// Access to data members:
	//---------------------------------------------------------------

	const IdType & getId() const { return id; }
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
	size_t getNNonFndrs() const { return (nMembers - nFounders); }	///< Number of non-founders

	const InheritanceSpace & getIVSpace() const { return *ivSpace; }
	size_t getNFounderGametes() const { return (nFounders << 1   ); } ///< Number of founder gametes
	size_t getNMeiosis	 () const { return (getNNonFndrs()<<1); } ///< Number of meiosis
							///< NB: see haploid note in recurseSib().



	/// Array-style access: @a mIdx must be between 0 and getNMembers()-1.
	const Member & memberAt( MemberIdx mIdx ) const
	    {
	    #if AGGRESSIVE_RANGE_CHECK
		if ( mIdx >= getNMembers() )
		    throwMRange( mIdx );
	    #endif

	    return *(sortedMembers[mIdx]);
	    }


	/// Array-style access: @a fIdx must be between 0 and getNFounders()-1.
	const Member & founderAt( FounderIdx fIdx ) const
	    {
	    #if AGGRESSIVE_RANGE_CHECK
		if ( fIdx >= getNFounders() )
		    throwFRange( fIdx );
	    #endif

	    return *(sortedMembers[fIdx]);
	    }



	//---------------------------------------------------------------
	// Mendelian error access:
	//---------------------------------------------------------------

	int  getNMendelErrs() const { return nMendelErrs; }
	bool haveMendelErrAt( SLocIdxType t ) const;



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
	/// using genPossibleStatesInternal() prior to calling.  This could be
	/// in a separate class.
	const HiddenStateSpace & getStateProbs( SLocIdxType sLocIdx ) const;

	/// Release any memory consumed by generating the states.  After calling
	/// this, genPossibleStatesInternal() must then be re-called before
	/// getStateProbs() can again be called.
	void releaseStateProbs() const;


	// These are not guaranteed to produce any output if ostream support is
	// not compiled into Pedigree, InheritanceVector, AncestryVector, etc.:
	static void dbgRecursion( bool nv ); ///< Request that recursion debugging information output to cout
	static void dbgEmission ( bool nv ); ///< Request that emission-probability debugging output to cout



    //=========================================================================
    //
    // See NOTE *4*:
    // Additional data members and methods to run the admixmap model.  These
    // should not be here in the base class, but rather should be moved into a
    // derived class in admixmap, not here in libbase.  This also applies to
    // deriving from PedBase, which should take place as multiple
    // inheritance in AdmixPedigree, like this:
    //
    // class AdmixPedigree : public Pedigree , public PedBase
    //  { ... };
    //
    // See also NOTE X1 above.
    //
    //=========================================================================

    protected:
	void SetUniformAdmixtureProps(); // Set the initial values of theta


    public:

	/// This should be the constructor of the derived AdmixedPedigree class.
	void InitialiseAdmixedStuff( const AdmixOptions & options );

	void SampleTheta(
		 int				    iteration	      , // in
		 double *			    SumLogTheta       , // out (accumulate with existing value)
		 const bclib::DataMatrix *	    Outcome	      , // in
		 const DataType *		    OutcomeType       , // in
		 const std::vector<double> &	    lambda	      , // in
		 int				    NumCovariates     , // in
		 bclib::DataMatrix *		    Covariates	      , //
		 const std::vector<const double*> & beta	      , //
		 const PopThetaType &		    poptheta	      , // in
		 const AdmixOptions &		    options	      , // in
		 const AlphaType &		    alpha	      , // in?
		 double				    DInvLink	      , // in
		 double				    dispersion	      , // in
		 CopyNumberAssocTest &		    ancestryAssocTest , // Not used by pedigrees
		 bool				    RW		      , // in (is-random-walk)
		 bool				    anneal	      );// in

	void HMMIsBad( bool loglikisbad );

	void SampleHiddenStates ( const Options & options );
	void SampleJumpIndicators ( bool sampleArrivals );


    private:

	size_t nAffected;

	/// Returns size of Theta, used to allocate storage for the various
	/// thetas and SumSoftmaxTheta.  Was: NumGametes --> K*NumGametes
	typedef FounderIdx ThetaIdx;
	ThetaIdx getNTheta() const { return getNFounders(); }



	//---------------------------------------------------------------------
	/// Admixture proportions, one vector of size K for each founder (in the
	/// case of pedigrees, we assume that the admixture proportion is the
	/// same for both gametes in each founder).  Rather than keeping a
	/// separate "Theta" and "ThetaProposal", and shifting new values from
	/// the proposed to the current vectors (as in the individual version),
	/// we'll use pointers and create a sort of ring buffer with ability to
	/// propose a new theta, run the HMM for it, and if we decide to accept
	/// it just update the pointer and free up the old "current" vector to
	/// be re-used for the next proposal, but if we reject it, just free up
	/// the proposed vector for re-use as the next proposal.  In addition to
	/// eliminating the shifting around of theta values in memory, it
	/// _vastly_ simplifies the code, and makes it much more robust.
	//---------------------------------------------------------------------

	ProposalRingBuffer<ThetaType> thetas;
	const ThetaType & getCurTheta() const { return thetas.getCurrent(); } ///< Convenience.
	ThetaType &	  getCurTheta()	      { return thetas.getCurrent(); } ///< Convenience.
	double getCurThetaVal( FounderIdx f, PopIdx k ) const { return getCurTheta()[f][k]; } ///< Convenience.


	//---------------------------------------------------------------------
	/// Structure to hold the cached values of log-likelihood from the HMM.
	/// Makes explicit the database concepts of "propose", "reject",
	/// "accept".  It would be nice if we could integrate this with the
	/// ProposalRingBuffer used for theta and rho, although it currently
	/// contains more logic for refreshing the values from the HMM when
	/// necessary.
	//---------------------------------------------------------------------

	class LLBuf
	    {
	    private:

		// Since the HMM can in theory change (and is also not yet allocated at the time
		// that the llCache is constructed), we will rather keep a reference to the
		// pedigree and get the HMM from it.
		#if 0
		    genepi::HiddenMarkovModel hmm;
		    genepi::HiddenMarkovModel & getHMM() { return hmm; }
		#else
		    Pedigree & ped;
		    genepi::HiddenMarkovModel & getHMM() { return ped.getHMM(); }
		#endif

		double accVal	    ; ///< Last accepted value
		bool   accValIsValid; ///< Does curVal correspond to the last accepted
				      ///< {rho,theta} parameters (it should always
				      ///< be so after the first calculation
				      ///< if the LogLik is calculated to decide if
				      ///< the proposed values will be accepted, and new values
				      ///< are always proposed before being accepted).

		bool   propActive   ; ///< Is a proposal currently active?
		double propVal	    ; ///< The proposed value, if (propActive && propIsValid)

		bool   propIsValid  ; ///< Does propVal correspond to the current values of
				      ///< {rho,theta}? (in fact, only if propActive also turned
				      ///< on).  This is a little hard to track because it
				      ///< requires that we know when rho and/or theta becomes
				      ///< dirty, yet not too bad because we can mark it as dirty
				      ///< when a proposal begins, and once it is computed for a
				      ///< proposal, we rarely change the parameters before
				      ///< accepting or rejecting it.  In fact, we rarely ask for
				      ///< it more than once, so it really doesn't even need to be
				      ///< cached here.

	    public:

		LLBuf( Pedigree & _ped );
		LLBuf( Pedigree & _ped, const LLBuf & rhs );
		LLBuf & operator=( const LLBuf & rhs );

		void startProposal();
		void acceptProposal();
		void rejectProposal();
		double getAcceptedVal() /*const*/;
		double getProposedVal() /*const*/;


		/// Get the log-likelihood of the "current" parameters (rho,theta): if a proposal
		/// is active, then the proposed value, otherwise the last accepted value.
		double getCurVal() /*const*/
		    {
		    return propActive ? getProposedVal() : getAcceptedVal();
		    }


		/// Is a proposal currently underway?
		bool isProposalActive() const { return propActive; }

		void parmsChanged();
	    };




	mutable ThetaType SumSoftmaxTheta; ///< Mutable so that WritePosteriorMeans() can modify it.


	ProposalRingBuffer<RhoType> rhos;
	const RhoType & getCurRho() const { return rhos.getCurrent(); } ///< Convenience.
	RhoType &	getCurRho()	  { return rhos.getCurrent(); } ///< Convenience.


	RhoType sumlogrho ; ///< foo.

	/// Accessor for global-rho, which is stored in the first element of the rho array.
	/// NOTE: This may be directly accessed in PopAdmix, we also need to convert it to use this accessor method.
	/// NOTE: Because it may be important for performance, we implement this
	///	as an inline method here; but since the "Options" class (which
	///	is needed to assert that global-rho is turned on) is a forward
	///	definition, this method is _not_ inline in the case of
	///	aggressive assertions.
	#if AGGRESSIVE_RANGE_CHECK
	    double & getGlobalRhoVal();
	#else
	    double & getGlobalRhoVal() { return getCurRho()[0]; }
	#endif

	/// Const version of getGlobalRhoVal()
	// This is a little scary in that could easily result in infinite
	///	recursion; perhaps we should give explicitly different names to
	///	the const and non-const versions rather than relying on
	///	overloading.
	double getGlobalRhoVal() const { return const_cast<Pedigree*>(this)->getGlobalRhoVal(); }


	double		     step	     ;
	unsigned int	     NumberOfUpdates ;
	unsigned int	     w		     ;
	bclib::StepSizeTuner ThetaTuner	     ;
	unsigned int	     NumGametes	     ;
	unsigned int	     myNumber	     ;

	#if PED_HAS_OWN_PRNG
	    genepi::NRand rng; ///< See getRNG().
	#endif

	// For affected-only test: maintain precomputed factors
	typedef float AAProbType;
	cvector<AAProbType> aaInfo  ; // Indexed on HSS-non0-index
	cvector<AAProbType> aaScore ;
	void accumAAScore();


	mutable genepi::TransProbCache * tpCache;
	genepi::TransProbCache & getTPC() const;

	mutable genepi::HiddenMarkovModel * hmm;
	genepi::HiddenMarkovModel & getHMM() const;
	void freeHMM() const;

	LLBuf llCache;

    public:

	MemberIdx    getNAffected() const { return nAffected; }
	unsigned int getMyNumber () const { return myNumber  ; } ///< "number" of this pedigree, counting from 1
	unsigned int getIndex    () const { return myNumber-1; } ///< "number" of this pedigree, counting from 0
	void setMyNumber( unsigned int nv );

	#if PED_HAS_OWN_PRNG
	    /// This is one way to solve the non-thread-safety of bclib::Rand:
	    /// keep a separate PRNG for each pedigree.  We really only need one
	    /// per thread rather than one per pedigree, so we could also try
	    /// making Rand's global instance threadprivate; yet a PRNG is
	    /// typically not a large object, so there is little downside to
	    /// this approach.
	    NRand & getRNG() { return rng; }
	#endif



	//--------------------------------------------------------------------
	// Combined proposal methods
	//--------------------------------------------------------------------

	ThetaType & startThetaProposal();
	void	    rejectThetaProposal();
	void	    acceptThetaProposal();


    private:

    double ProposeThetaWithRandomWalk( const AlphaType & alpha );

    double LogAcceptanceRatioForRegressionModel( RegressionType RegType, bool RandomMatingModel,
					       PopIdx K, int NumCovariates,
					       const bclib::DataMatrix * Covariates, const double * beta,
					       double Outcome, const PopThetaType & poptheta, double lambda) const;
    void Accept_Reject_Theta( double logpratio, bool isRandomMatingModel, bool isRandomWalk );
    void updateAdmixtureForRegression( int NumCovariates, const PopThetaType & poptheta,
				     bool ModelIndicator, bclib::DataMatrix * Covariates);

    public:
	/// Overridden from PedBase, public method, manages cached values and
	/// calls calcLogLikelihood() if necessary.
	double getLogLikelihood( const Options & options, bool forceUpdate, bool store );

	/// Called if a Metropolis proposal is accepted.  The parameter is
	/// ignored in the case of pedigrees.
	virtual void storeLogLikelihood( const bool setHMMAsOK );

	//--------------------------------------------------------------------------
	// Public rho-proposal methods, overridden from PedBase, called from PopAdmix, ignored for individuals.
	virtual void setRho( double nv );
	virtual void startRhoProposal();
	virtual void acceptRhoProposal();
	virtual void rejectRhoProposal();
	//--------------------------------------------------------------------------

    private:

    double getLogLikelihoodAtPosteriorMeans( const Options & options );
    double getLogLikelihoodOnePop() const;
    void setAdmixtureProps( const ThetaType & rhs );



    // Methods (overridden from PedBase) from AdmixedIndividual:
    virtual void SetGenotypeProbs(int j, int jj, unsigned locus, bool chibindicator);
    virtual void drawInitialAdmixtureProps(const AlphaType &alpha);

    // ResetSufficientStats() is only needed if we are doing conjugate updates
    // (which we don't do for pedigrees), as well as when we do not do conjugate
    // updates.
    virtual void ResetSufficientStats();

    void SampleHapPair(unsigned chr, unsigned jj, unsigned locus, AlleleFreqs *A,
			  bool skipMissingGenotypes, bool annealthermo, bool UpdateCounts);
    void SampleMissingOutcomes( bclib::DataMatrix * Outcome, const vector<bclib::Regression*> & R );

    void UpdateScores( const AdmixOptions & options, bclib::DataMatrix * Outcome,
			bclib::DataMatrix * Covariates,
			const vector<bclib::Regression*> & R, AffectedsOnlyTest & affectedsOnlyTest,
			CopyNumberAssocTest & ancestryAssocTest );

    bool GenotypeIsMissing	( unsigned int clIdx ) const; ///< locus is a comp locus
    bool simpleGenotypeIsMissing( SLocIdxType  slIdx ) const; ///< locus is a simple locus
    bool isHaploidatLocus(unsigned j)const;
    bool isHaploidIndividual()const;


    /// Copied from AdmixedIndividual because it is used by WritePosteriorMeans().
    void getPosteriorMeans( ThetaType & ThetaMean, RhoType & rhoMean /* output parameter? */, unsigned int samples ) const;


    /// Overridden from PedBase (body is adapted from code copied from
    /// AdmixedIndividual).  Perhaps should be moved up to base class.
    void WritePosteriorMeans( ostream& os, unsigned samples, bool globalrho ) const;


    /// Cache for getNInheritedByAffected()
    mutable TwoDimArray<short,PopIdx,FounderIdx> * aoCache;
    /// Helper method for getNInheritedByAffected (affected-only test computations).
    short calcNInheritedByAffected( PopIdx k, FounderIdx fIdx, const AncestryVector & av, const InheritanceVector & iv ) const;
    /// Helper method for affected-only test computations.
    int getNInheritedByAffected( PopIdx k, FounderIdx fIdx, const AncestryVector & av, const InheritanceVector & iv ) const;

    /// Supports affected-only test computations.
    void accumAOScore( AffectedsOnlyTest & aoTest ) const;


    /// This is very bad; we need the concept of context.
    static const AdmixOptions * optionsRef;
    public:
	/// This is very bad; we need the concept of context.
	static const AdmixOptions & getOptions() { gp_assert(optionsRef != 0); return *optionsRef; }
	static void setOptions( const AdmixOptions & opts );


    private:

    // ====== DEBUGGING METHODS (overridden from PedBase) ======
    #if PEDBASE_DEBUG_METHODS
	virtual void dumpTheta( const char * prefix ) const;
    #endif
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
