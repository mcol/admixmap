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
/// \file PedigreeAddl.h
/// Some additional functions for Pedigrees, outside of the class.
//=============================================================================

#ifndef __base_PedigreeAddl_h
#define __base_PedigreeAddl_h


#include <iosfwd>
#include <vector>

#include "Pedigree.h"



namespace genepi { // ----

/** \addtogroup base
 * @{ */



/// Output a summary of stats about a pedigree -- the number of members, number of founders, etc.
std::ostream & print_aggregate_summary( std::ostream & os, const std::vector<Pedigree> & peds );



/** @} */

} // ---- end namespace genepi



#endif // ! __base_PedigreeAddl_h
