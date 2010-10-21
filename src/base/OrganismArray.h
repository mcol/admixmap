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
/// \file OrganismArray.h
/// Definition of the genepi::OrganismArray typedef.
//=============================================================================

#ifndef __base_OrganismArray_h
#define __base_OrganismArray_h



#include "GenotypeParser.h"



namespace genepi { // ----




/** \addtogroup base
 * @{ */



/// We temporarily used the GenotypeParser as the container to store the
/// Organism; they need their own container class, and the parser should only
/// be a reader. This typedef may facilitate that conversion by allowing
/// separate names to be used for the parser and container.

typedef GenotypeParser OrganismArray;




} // ---- end namespace genepi



/** @} */



#endif // ! __base_OrganismArray_h
