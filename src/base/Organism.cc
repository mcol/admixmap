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
/// \file Organism.cc
/// Implementation of the Organism class.
//=============================================================================

#include "Organism.h"

#include "GenotypeParser.h"

#include <cstring>  // memcpy()



namespace genepi { // ----



//=============================================================================
// These are not inlined in the header to avoid circular-dependecy with the parser:
//=============================================================================

SLocIdxType Organism::getNSimpleLoci() const
    {
    return inFile.getNSimpleLoci();
    }

const SimpleLocusArray & Organism::getSLoci() const
    {
    return inFile.getSLoci();
    }



//=============================================================================
// Organism
//=============================================================================

Organism::Organism( const Organism & rhs ) :
	lineNum ( rhs.lineNum	) ,
	famId	( rhs.famId	) ,
	orgId	( rhs.orgId	) ,
	sex	( rhs.sex	) ,
	outcome ( rhs.outcome	) ,
	father	( rhs.father	) ,
	mother	( rhs.mother	) ,
	depth	( rhs.depth	) ,
	inFile	( rhs.inFile	)
    {
    const SLocIdxType nLoci = inFile.getNSimpleLoci();
    gp_assert( nLoci != 0 );
    gtypes = new Genotype[ nLoci ];

    // NB: This will not be safe if Genotype ever has a constructor (see NOTE *3*
    // in Genotype.h)
    memcpy( gtypes, rhs.gtypes, nLoci * sizeof(*gtypes) );
    }


Organism::Organism( const GenotypeParser & gp ) :
	inFile ( gp )
    {
    gp_assert( inFile.getNSimpleLoci() != 0 );
    gtypes = new Genotype[ inFile.getNSimpleLoci() ];
    }


Organism::~Organism()
    {
    delete[] gtypes;
    }


Organism & Organism::operator=( const Organism & rhs )
    {
    famId   = rhs.famId	  ;
    orgId   = rhs.orgId	  ;
    father  = rhs.father  ;
    mother  = rhs.mother  ;
    depth   = rhs.depth   ;
    sex     = rhs.sex     ;
    outcome = rhs.outcome ;

    // We can't copy a reference, so we assert that they already refer to the
    // same object.  This violates encapsulation and effectively requires that
    // the lhs and rhs came from the same input file.
    gp_assert( &inFile == &rhs.inFile );

    const SLocIdxType nLoci = inFile.getNSimpleLoci();
    gp_assert( nLoci != 0 );

    gtypes = new Genotype[ nLoci ];

    // NB: This will not be safe if Genotype ever has a constructor (see NOTE
    // *3* in GFileLexer.h)
    memcpy( gtypes, rhs.gtypes, nLoci * sizeof(*gtypes) );

    return *this;
    }



//-----------------------------------------------------------------------------
// throwError() [protected]
//-----------------------------------------------------------------------------

void Organism::throwError( const std::string & msg ) const
	//throw( ParseError )
    {
    throw GFileLexer::ParseError( msg, getInFile().getFileName(), getLineNum() );
    }



//-----------------------------------------------------------------------------
// getGType()
//-----------------------------------------------------------------------------

const Genotype & Organism::getGType( SLocIdxType sLocIdx ) const
    {
    if ( sLocIdx >= getNSimpleLoci() )
	throw std::invalid_argument( estr("Genotype index ") + sLocIdx +
			" exceeds number of loci (" + getNSimpleLoci() + ')' );

    return gtypes[ sLocIdx ];
    }



//-----------------------------------------------------------------------------
// idDesc()
//-----------------------------------------------------------------------------

estr Organism::idDesc() const
    {
    estr rv;
    if ( inFile.isPedFile() )
	rv << "{fam:" << famId << ";org:" << orgId << '}';
    else
	rv << famId;
    return rv;
    }


//-----------------------------------------------------------------------------
// inLineDesc()
//-----------------------------------------------------------------------------

estr Organism::inLineDesc() const
    {
    estr rv( getInFile().getFileName() );
    return (rv << ':' << getLineNum());
    }



} // ---- end namespace genepi
