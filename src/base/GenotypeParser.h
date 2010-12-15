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
/// \file GenotypeParser.h
/// Definition of the genepi::GenotypeParser class.
//=============================================================================


#ifndef __base_GenotypeParser_h
#define __base_GenotypeParser_h



#include <cstddef>	// size_t
#include <map>		// multimap

#include "GFileLexer.h"
#include "Organism.h"
#include "PairSecIter.h"
#include "SimpleLocusArray.h"



namespace genepi { // ----



/** \addtogroup base
 * @{ */



//-----------------------------------------------------------------------------
//
/// Class to read and parse a genotype or pedigree file, cooking the data a
/// little.  It also serves as a container for the data while the haplotype
/// probabilities are extracted.  This of course should be broken into two
/// classes, a parser and a container.
///
/// See also
/// genepi::convert(const Organism&org,const Genome&loci,std::vector<genotype> &genotypes,bool**missing)
//
//-----------------------------------------------------------------------------

class GenotypeParser : public GFileLexer
    {

    public:

	typedef std::vector<Organism>	 RowsType     ;
	typedef RowsType::iterator	 RowIter      ;
	typedef RowsType::const_iterator ConstRowIter ;


	static const OrgIdType MISSING_PARENT_ID  = 0; ///< The input-file indicator for no-parent
	static const OrgIdType NON_PEDFILE_ORG_ID = 0; ///< The organism-ID used internally for non-pedfiles


    private:
	const SimpleLocusArray & simpleLoci;

	int	   gtypeOffset	 ; ///< Column index of first genotype
	bool	   hasSexFlag	 ;
	bool	   isPedFileFlag ;
	RowsType   rows		 ;

	/// A multi-valued map from family-ID to the members of that family.  In
	/// addition to implementing rapid search for individual
	/// (family-ID,individual-ID) tuples, this allows us to locate all of
	/// the members of a given pedigree.
	/// Note: each insertion places the element according to the default
	/// comparison function over FamIdType (that is, the order is defined
	/// by string comparison).
	std::multimap<FamIdType,Organism*> pmap;


	// Private methods used to cook the data after loading:
	void buildAndCheckPedGraphs();


    protected:
	Organism::SexType lexSex();


	void insert_pmap( Organism & org )
	    { pmap.insert( std::multimap<FamIdType,Organism*>::value_type( org.getFamId(), &org ) ); }


	void validateOrg( Organism & org, bool maleXHomozygWarn ) const;


    public:

	// File-format queries:
	bool hasSexColumn () const { return hasSexFlag    ; }
	bool isPedFile    () const { return isPedFileFlag ; }


	//---------------------------------------------------------------
	// Container-style methods: access to Organism records
	// (iteration/search): All "array-style" access is 0-based.
	//---------------------------------------------------------------

	/// Search by family/organism ID; throws an exception if not found.
	const Organism & findById( const FamIdType & famId, const OrgIdType & orgId ) const;

	/// Search by family/organism ID; returns null pointer if not found.
	const Organism * findByIdIfExists( const FamIdType & famId, const OrgIdType & orgId ) const;

	/// Test existence of a given family/organism ID:
	bool exists( const FamIdType & famId, const OrgIdType & orgId ) const
		{ return (findByIdIfExists( famId, orgId ) != 0); }


	/// Range-checked array-index access to organisms; for extra-fast
	/// access, see getOrganismUnsafe().  For easier syntax, see
	/// operator[]().
	const Organism & getOrganism( size_t orgIdx ) const { return rows.at(orgIdx); }

	// Currently, this is only needed for InputData to "merge" the
	// outcome-file back into the organisms; not needed for pedigree files.
	Organism & getOrganism( size_t orgIdx ) { return rows.at(orgIdx); }

	/// Range-checked array-index access to organisms; for extra-fast
	/// access, see getOrganismUnsafe().
	const Organism & operator[]( size_t orgIdx ) const { return rows.at(orgIdx); }

	/// Use this for extra-fast (inline and non-range-checked) array access
	/// to the organism records &mdash; e.g. inside inner loops, but in
	/// general, the range-checked operator[] is a better choice for
	/// non-performance-critical applications.
	const Organism & getOrganismUnsafe( size_t orgIdx ) const { return rows[orgIdx]; }

	size_t getNumOrganisms	     () const { return rows.size(); }
	size_t getNumberOfIndividuals() const { return getNumOrganisms(); } ///< compatibility


	// Access to pedigrees as a collection of organisms:

	// Pedigree subset manipulation:

	//typedef MMVIter<FamIdType,Organism*> PedIter;
	typedef MMVIter (FamIdType,Organism*) PedIter;
	typedef MMVCIter(FamIdType,Organism*) ConstPedIter;

	// Iterate over entire file in Pedigree order:
	ConstPedIter beginByPed () const { return ConstPedIter(pmap.begin()); }
	ConstPedIter endByPed	() const { return ConstPedIter(pmap.end  ()); }

	std::pair<ConstPedIter,ConstPedIter> findOrgsInPed( const FamIdType & famId ) const
	    {
	    std::pair<std::multimap<FamIdType,Organism*>::const_iterator,
		      std::multimap<FamIdType,Organism*>::const_iterator> p = pmap.equal_range(famId);
	    return std::pair<ConstPedIter,ConstPedIter>( ConstPedIter(p.first), ConstPedIter(p.second) );
	    }

	void findOrgsInPed( const FamIdType & famId, ConstPedIter & begin, ConstPedIter & end ) const
	    {
	    std::pair<std::multimap<FamIdType,Organism*>::const_iterator,
		      std::multimap<FamIdType,Organism*>::const_iterator> p = pmap.equal_range(famId);
	    begin = p.first;
	    end = p.second;
	    }

	/// Returns a ConstPedIter pointing to the first member of the pedigree,
	/// and puts an iterator pointing to the last member of the pedigree in
	/// <CODE>end</CODE>.
	ConstPedIter findOrgsInPed( const FamIdType & famId, ConstPedIter & end ) const
	    {
	    std::pair<std::multimap<FamIdType,Organism*>::const_iterator,
		      std::multimap<FamIdType,Organism*>::const_iterator> p = pmap.equal_range(famId);
	    end = p.second;
	    return p.first;
	    }


	/// Returns the number of organisms in pedigree <B>famId</B>
	size_t numOrgsInPed( const FamIdType & famId ) const
	    {
	    return pmap.count( famId );
	    }


	/// Get the container of simple loci
	const SimpleLocusArray & getSLoci() const { return simpleLoci; }

	/// Convenience: returns the number of loci in the locus file
	size_t getNSimpleLoci() const { return getSLoci().size(); }


	/// Run through the input data after the initial parsing has finished,
	/// checking the dataset as a whole for some basic validation.
	void validate() const;


	//---------------------------------------------------------------
	// Constructors/destructor:
	//---------------------------------------------------------------

	/// Constructs the object, opens the file, and parses the data:
	GenotypeParser( const char * fileName, const SimpleLocusArray & sLoci,
			bool outcomeIsBinary, bool maleXHomozygWarn );
	~GenotypeParser();

    };



/** @} */



} // ---- end namespace genepi



#endif // ! __base_GenotypeParser_h
