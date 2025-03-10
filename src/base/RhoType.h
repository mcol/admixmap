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
/// \file RhoType.h
/// Definition of the genepi::RhoType typedef.
//=============================================================================


#ifndef __base_RhoType_h
#define __base_RhoType_h



#include <bclib/cvector.h>



namespace genepi { // ----


/** \addtogroup base
 * @{ */



/// Type definition for "rho".  Vector of doubles, indexed on ???
typedef cvector< double > RhoType;




/** @} */


} // ---- end namespace genepi



#endif // ! __base_RhoType_h
