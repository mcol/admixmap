//=============================================================================
//
// Copyright (C) 2009  David D. Favro  gpl@meta-dynamic.com
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
/// \file Exceptions.cc
//=============================================================================

#pragma implementation
#include "exceptions.h"



#include <cerrno>
#include <cstring>  // strerror()
#include <stdexcept>
#include <string>



namespace genepi { // ----



void throwSysErr( const char * sc )
    {
    throw std::runtime_error( std::string(sc) + " failed: " + strerror(errno) );
    }



} // end namespace genepi ----
