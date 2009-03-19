//=============================================================================
//
// Copyright (C) 2009  David D. Favro  gpl-copyright@meta-dynamic.com
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
/// \file AlleleFreqParser.h
/// Definition of the AlleleFreqParser class.
//=============================================================================

#ifndef __base_AlleleFreqParser_h
#define __base_AlleleFreqParser_h



#include "AlleleArray.h"
#include "GFileLexer.h"
#include "SimpleLocusArray.h"

#include <vector>



namespace genepi { // ----



/** \addtogroup base
 * @{ */



class AlleleFreqParser : public GFileLexer
    {
    public:
	typedef std::vector<std::string> PopArray;

    private:
	const SimpleLocusArray & loci	     ; ///< Keep a reference to the loci
	AlleleProbVect &	 outProb     ; ///< The output: normalized probabilities
	PopArray &		 populations ; ///< The header line (output)


    protected:
	void checkLocusComplete( Genotype::AlleleType lastAllele,
				 const SimpleLocus *  lastLocus ) const;

    public:

	const SimpleLocusArray & getLoci   () const { return loci   ; }
	const AlleleProbVect &	 getOutProb() const { return outProb; }


	//---------------------------------------------------------------
	// Constructors/destructor:
	//---------------------------------------------------------------

	/// Constructs the object, but doesn't yet open the file; use parse() to
	/// process the file.  See also the static method @sa{parse(const char
	/// *,PopArray&,AlleleProbVect&)} for a more convenient method.
	AlleleFreqParser( const char * fileName, const SimpleLocusArray & _loci,
			    PopArray & _outPop, AlleleProbVect & _outProb );

	/// Opens the file, parses the data, normalizes the frequencies into the
	/// output AlleleProbVect, and closes the file.
	void parse();

	/// Static convenience method to construct the parser object and parse the file.
	static void parse( const char * fileName, const SimpleLocusArray & _loci,
			    PopArray & _outPop, AlleleProbVect & _outProb )
	    { AlleleFreqParser(fileName,_loci,_outPop,_outProb).parse(); }
    };



} // ---- end namespace genepi



#endif // ! __base_AlleleFreqParser_h
