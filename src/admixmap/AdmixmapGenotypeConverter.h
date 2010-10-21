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
/// \file AdmixmapGenotypeConverter.h
/// Definition of the
/// genepi::convert(const Organism&org,const Genome&loci,std::vector<genotype> &genotypes,bool**missing)
/// function.
//=============================================================================

#ifndef __base_AdmixmapGenotypeConverter_h
#define __base_AdmixmapGenotypeConverter_h



#include "GenotypeParser.h"
#include "Genome.h"

#include <vector>



/** \addtogroup base
 * @{ */



namespace genepi { // ----



//-----------------------------------------------------------------------------
//
/// Convert an Organism's genotype data to the vector of genotypes and
/// missing-array that is used by Model / AdmixMapModel.
///
/// This is a temporary function to bridge between the "old" (GenotypeLoader)
/// and the "new" (GenotypeParser) method of external and internal data
/// representation (input-data-file parsing).  It replaces and is equivalent to
/// GenotypeLoader::GetGenotype(), but is a separate function rather than a
/// method, in part because it is admixmap-specific, so should not be included
/// in hapmixmap by being placed in the base library
//
//-----------------------------------------------------------------------------

void convert( const Organism & org, const Genome & loci, std::vector<genotype> & genotypes, bool** missing );



} // ---- end namespace genepi



/** @} */



#endif // ! __base_AdmixmapGenotypeConverter_h
