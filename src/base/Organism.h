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
/// \file Organism.h
/// Definition of the genepi::Organism class.
//=============================================================================

#ifndef __base_Organism_h
#define __base_Organism_h



#include "SimpleLocusArray.h"
#include "bclib/estr.h"

#include <cstddef>	// size_t
#include <climits>	// SHRT_MAX
#include <list>



namespace genepi { // ----

/** \addtogroup base
 * @{ */



class Genotype;
class GenotypeParser;


typedef unsigned int OrgIdType ;
typedef std::string  FamIdType ;
typedef size_t	     PopIdx    ;



//-----------------------------------------------------------------------------
//
/// Class to hold an organism (that is, one row in a genotype or pedigree input file).
///
/// Organisms which were read in from a pedigree file contain a family ID (getFamId()) and
/// an organism ID; the pair uniquely identifies an individual in the file,
/// i.e. the organism ID is only unique within the pedigree.
///
/// <A name="note-1"></A>
/// <TABLE STYLE="border: groove 3pt aqua;">
///  <TR>
///	<TD><B>NOTE *1*</B></TD>
///	<TD>
///	We implement an ordering for iteration such that no member of a family
///	will be visited prior to both of its parents being visited (see
///	<A HREF="classgenepi_1_1Pedigree.html#note-1">NOTE *1* in Pedigree</A>).
///	</TD>
///  </TR>
/// </TABLE>
///
/// <A name="note-2"></A>
/// <TABLE STYLE="border: groove 3pt aqua;">
///  <TR>
///	<TD><B>NOTE *2*</B></TD>
///	<TD>
///	Once the list of pointers to records in this sorted-order is created (in
///	the Pedigree object -- see
///	<A HREF="classgenepi_1_1Pedigree.html#note-2">NOTE *2* in Pedigree</A>),
///	the depth field is no longer needed, so the data member here is re-used
///	as the index within the sorted-list in the Pedigree.  This makes it much
///	easier to manipulate the digraph of the members within the pedigree in
///	that it avoids the need to replicate the father and mother pointers in
///	the pedigree's external data structure (ordered list), but we still can
///	easily follow a parent-pointer and then find the resulting index in the
///	sorted-list, as well as any "parallel arrays" that share the same
///	indexing.
///	</TD>
///  </TR>
/// </TABLE>
///
/// <A name="note-3"></A>
/// <TABLE STYLE="border: groove 3pt aqua;">
///  <TR>
///	<TD><B>NOTE *3*</B></TD>
///	<TD>
///	Because the genotype data is dynamically allocated without
///	reference-counting, copying these objects is an expensive operation.
///	This is not a problem because they are currently managed by the
///	OrganismArray class; they are read in from the input file by the
///	GenotypeParser class which places them into the OrganismArray, which
///	will destroy them when it is destroyed; in the meantime, other object
///	which wish to use them should keep references or pointers, rather than
///	making copies.
///	</TD>
///  </TR>
/// </TABLE>
///
//
//-----------------------------------------------------------------------------

class Organism
    {

    friend class GenotypeParser;


    public:
	enum SexType
	    {
	    SEX_UNKNOWN = 0,
	    SEX_MALE	= 1,
	    SEX_FEMALE	= 2
	    };


	enum OutcomeType
	    {
	    OUTCOME_UNKNOWN	,
	    OUTCOME_UNAFFECTED	,
	    OUTCOME_AFFECTED	,
	    };


    private:

	/// Input file row number, for formatting error messages
	int	    lineNum  ;

	FamIdType   famId    ;
	OrgIdType   orgId    ;

	SexType	    sex	     ;
	OutcomeType outcome  ;
	Genotype *  gtypes   ;
	bool	    gtypedFlg;

	/// Model this organism as a single haploid gamete.
	bool singleGameteModel;


	// If we allow alphanumeric organism-IDs, and store as
	// std::string's (POD type), we cannot put them in a union with
	// the parent-pointer (punning).  In that case, we might leave
	// them always allocated here, and clear the dynamic storage
	// once no longer needed.
	union {
	    Organism * father	 ; ///< 0 for blank (missing, unknown)
	    OrgIdType  fatherId; };
	union {
	    Organism * mother	 ; ///< 0 for blank (missing, unknown)
	    OrgIdType  motherId; };

	short depth; ///< See NOTE *1*

	/// Keep a linked-list of the children; while this is not
	/// strictly necessary, it simplifies some algorithms
	std::list<Organism*> children;


	/// It is useful, given an organism, to obtain a reference to
	/// the input-file from which it was parsed, e.g. to format
	/// error messages.  There is typically only one
	/// genotype/pedigree input file, so we could use a quasi-global
	/// reference to the input file object, but instead keep a
	/// reference here for better encapsulation.  we can easily pass
	/// around references organisms yet access the GenotypeParser.
	const GenotypeParser & inFile;


    protected:

	// Because we use two passes, we initially store the
	// father/mother IDs in the place of the pointers (AKA punning);
	// a second pass dereferences them.  This is potentially
	// breakable on certain architectures, but the compiler should
	// notice it.
	void setFatherId( const OrgIdType & _fatherId ) { fatherId = _fatherId; }
	void setMotherId( const OrgIdType & _motherId ) { motherId = _motherId; }
	OrgIdType getFatherId() const { return fatherId; }
	OrgIdType getMotherId() const { return motherId; }

	// throw() is removed to avoid #include <ParseError>
	/// Throw a syntax-error exception
	void throwError( const std::string & msg ) const /*throw(ParseError)*/ GP_NO_RETURN;

    public:

	int getLineNum() const { return lineNum; } ///< Input file row number, for formatting error messages

	const FamIdType & getFamId() const { return famId ; } ///< Get the family-ID
	const OrgIdType & getOrgId() const { return orgId ; } ///< The organism-ID (unique within family)

	SexType getSex	() const { return sex; }
	bool	isFemale() const { return (sex == SEX_FEMALE ); } ///< true if female, false if male/unknown
	bool	isMale	() const { return (sex == SEX_MALE   ); } ///< true if male, false if female/unknown
	bool	sexKnown() const { return (sex != SEX_UNKNOWN); }

	OutcomeType getOutcome() const { return outcome; }
	void setOutcome( OutcomeType nv ) { outcome = nv; }

	/// Returns reference to the genotype at locus locIdx.
	const Genotype & getGType( SLocIdxType sLocIdx ) const;

	/// Is the organism modeled by a single gamete?
	bool isHaploid( IsXChromType is_xchrom ) const
	    {
	    return ((is_xchrom == CHR_IS_X) && isMale()) || singleGameteModel;
	    }
	void setSingleGameteModel( bool nv ) { singleGameteModel = nv; }

	/// Returns pointer to the father's record, or null if none/unknown.
	const Organism * getFather() const { return father; }

	/// Returns pointer to the mother's record, or null if none/unknown.
	const Organism * getMother() const { return mother; }

	bool hasFather() const { return (father != 0); }
	bool hasMother() const { return (mother != 0); }
	bool isFounder() const { return ((father == 0) && (mother == 0)); }

	static const short UNKNOWN_DEPTH = SHRT_MAX;
	void  setDepth( short nv ) { depth = nv; }  ///< See <A HREF="#note-1">NOTE *1*</A>
	short getDepth() const { return depth; }    ///< See <A HREF="#note-1">NOTE *1*</A>
	void  setPIdx( short nv ) { depth = nv; }   ///< See <A HREF="#note-1">NOTE *1*</A>
	size_t getPIdx() const { return depth; }    ///< Index within pedigree (see <A HREF="#note-1">NOTE *2*</A>)
	//^^^^ really Pedigree::MemberIdx

	const std::list<Organism*> & getChildren() const { return children; }
	std::list<Organism*> &	     getChildren()	 { return children; }
	typedef std::list<Organism*>::const_iterator ChConstIter;
	ChConstIter childrenBegin() const { return children.begin(); }
	ChConstIter childrenEnd	 () const { return children.end	 (); }
	bool hasChildren() const { return ! children.empty(); }


	/// Are all of the genotypes missing or do we have data for this
	/// individual: this can be inferred from the genotype data, or cached
	/// for the organism.
	bool isGenotyped() const { return gtypedFlg; }
	bool isMissing() const { return ! isGenotyped(); }
	///< Convenience: true if no genotyped data at all, i.e. &equiv; <CODE>!isGenotyped()</CODE>.

	const GenotypeParser & getInFile() const { return inFile; }


	/// Get the container of simple loci
	const SimpleLocusArray & getSLoci() const; // { return inFile.getSLoci(); }

	/// Convenience: returns the number of loci in the locus file
	SLocIdxType getNSimpleLoci() const;// { return inFile.getNSimpleLoci(); }

	// Constructors:
	//Organism( const GenotypeParser & gp ) : gtypes( 0 ), inFile( gp ) {}
	Organism( const GenotypeParser & gp );
	Organism( const Organism & rhs ); ///< See <A HREF="#note-3">NOTE *3*</A>
	~Organism();

	/// This violates encapsulation in that it doesn't copy the
	/// inFile reference; effectively requires that the lhs and rhs
	/// came from the same input file.
	Organism & operator=( const Organism & rhs );

	/// A human-readable identifier of the family/organism ID
	estr idDesc() const;

	/// A human-readable indicator of the file-name and row # from
	/// whence this organism came:
	estr inLineDesc() const;

    };



/** @} */

} // ---- end namespace genepi



#endif // ! __base_Organism_h
