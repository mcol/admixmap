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
/// \file estr.cc
//=============================================================================

#include "bclib/estr.h"


#include <cstdio>
#include <cstring>  // strcasecmp()

#include "bclib/exceptions.h"



#if defined(__GLIBC__) && ! (_BSD_SOURCE || _XOPEN_SOURCE >= 500 || _ISOC99_SOURCE)
    #error snprintf() not available
#endif



namespace genepi { // ----


static const size_t ES_BSIZE = 28;


#define ECON(T,FMT) \
    estr::estr( T x ) \
	{ \
	char buf[ES_BSIZE]; \
	const int rc = snprintf( buf, sizeof(buf), "%" #FMT, x ); \
	gp_assert( (rc>=0) && (rc<int(sizeof(buf))) ); \
	std::string::operator=( buf ); \
	}


ECON( short	    , hd )
ECON( unsigned short, hu )
ECON( int	    , d  )
ECON( unsigned int  , u  )
ECON( double	    , lf )

// See https://bugs.launchpad.net/ubuntu/+bug/221979
#if defined(__WIN32)
    estr::estr(		 long x ) { operator=( static_cast<	    int>( x ) ); }
    estr::estr( unsigned long x ) { operator=( static_cast<unsigned int>( x ) ); }
#else
    ECON( long int	    , ld )
    ECON( unsigned long int , lu )
#endif



bool estr::equalsCaseInsens( const estr & rhs ) const
    {
    return (strcasecmp( c_str(), rhs.c_str() ) == 0);
    }



bool equalsCaseInsens( const std::string & lhs, const std::string & rhs )
    {
    return (strcasecmp( lhs.c_str(), rhs.c_str() ) == 0);
    }



void stripWhitespace( std::string & str )
    {
    std::string::iterator iter;
    for ( iter = str.begin() ;
	    (iter != str.end()) && isblank(*iter) ;
	    ++iter )
	; // Do nothing
    if ( iter != str.begin() )
	str.erase( str.begin(), iter );

    std::string::reverse_iterator endIter;
    for ( endIter = str.rbegin() ;
	    (endIter != str.rend()) && isblank(*endIter) ;
	    ++endIter )
	; // Do nothing
    const std::string::size_type cutChars = (endIter - str.rbegin());
    if ( cutChars != 0 )
	str.resize( str.length() - cutChars );
    }




} // ---- end namespace genepi
