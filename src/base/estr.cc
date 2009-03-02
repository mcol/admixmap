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
/// \file estr.cc
//=============================================================================

#include "estr.h"


#include <cstdio>
#include <cstring>  // strcasecmp()

#include "exceptions.h"



#if defined(__GLIBC__) && ! (_BSD_SOURCE || _XOPEN_SOURCE >= 500 || _ISOC99_SOURCE)
    #error snprintf() not available
#endif



namespace genepi { // ----



static const int BSIZE = 30;


estr::estr( short x )
    {
    char buf[BSIZE];
    const int rc = snprintf( buf, sizeof(buf), "%hd", x );
    gp_assert( (rc>=0) && (rc<int(sizeof(buf))) );
    std::string::operator=( buf );
    }


estr::estr( unsigned short x )
    {
    char buf[BSIZE];
    const int rc = snprintf( buf, sizeof(buf), "%hu", x );
    gp_assert( (rc>=0) && (rc<int(sizeof(buf))) );
    std::string::operator=( buf );
    }


estr::estr( int x )
    {
    char buf[BSIZE];
    const int rc = snprintf( buf, sizeof(buf), "%d", x );
    gp_assert( (rc>=0) && (rc<int(sizeof(buf))) );
    std::string::operator=( buf );
    }


estr::estr( unsigned int x )
    {
    char buf[BSIZE];
    const int rc = snprintf( buf, sizeof(buf), "%u", x );
    gp_assert( (rc>=0) && (rc<int(sizeof(buf))) );
    std::string::operator=( buf );
    }


estr::estr( long x )
    {
    char buf[BSIZE];
    const int rc = snprintf( buf, sizeof(buf), "%ld", x );
    gp_assert( (rc>=0) && (rc<int(sizeof(buf))) );
    std::string::operator=( buf );
    }


estr::estr( unsigned long x )
    {
    char buf[BSIZE];
    const int rc = snprintf( buf, sizeof(buf), "%lu", x );
    gp_assert( (rc>=0) && (rc<int(sizeof(buf))) );
    std::string::operator=( buf );
    }


estr::estr( double x )
    {
    char buf[BSIZE];
    const int rc = snprintf( buf, sizeof(buf), "%lf", x );
    gp_assert( (rc>=0) && (rc<int(sizeof(buf))) );
    std::string::operator=( buf );
    }



bool estr::equalsCaseInsens( const estr & rhs ) const
    {
    return (strcasecmp( c_str(), rhs.c_str() ) == 0);
    }



bool equalsCaseInsens( const std::string & lhs, const std::string & rhs )
    {
    return (strcasecmp( lhs.c_str(), rhs.c_str() ) == 0);
    }




} // ---- end namespace genepi
