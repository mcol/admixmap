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
/// \file CodeTimer.cc
/// Implementation of the CodeTimer class.
//=============================================================================

#include "CodeTimer.h"

#include <stdexcept>
#include <cstring>	// strlen(), strerror()
#include <cerrno>
#include <sstream>
#include <iomanip>

using namespace std;



#if USE_TIME_T
	#define TM(x)	    (x)
	#define diff(x,y)   (x - y)
	#define DIFF_T	    time_t
#else
	#define TM(x)	    (x.tv_sec)
	#define diff(x,y)   (x.tv_sec - y.tv_sec + (double(x.tv_usec-y.tv_usec)/1000000.0))
	#define DIFF_T	    double
#endif



namespace genepi { // ----



//-----------------------------------------------------------------------------
// asctime_without_newline() [static, protected]
//-----------------------------------------------------------------------------

string CodeTimer::asctime_without_newline( const TimeType & tm )
   {
    #if (_POSIX_C_SOURCE || _XOPEN_SOURCE || _BSD_SOURCE || _SVID_SOURCE)
	char BUFFER[31];
	const char * const rv = asctime_r( localtime(&TM(tm)), BUFFER );
    #else
	const char * const rv = asctime( localtime(&TM(tm)) );
    #endif

    const size_t len = strlen( rv );

    #if 1 || USE_TIME_T
	return string( rv, len - 1 );
    #else
	ostringstream str;
	str.write( rv, len - 1 ) << '.' << setw(3) << setfill('0') << ((tm.tv_usec + 500) / 1000);
	return str.str();
    #endif
    }



//-----------------------------------------------------------------------------
// getNow() [static, protected]
//-----------------------------------------------------------------------------

CodeTimer::TimeType CodeTimer::getNow()
    {
    #if USE_TIME_T
	const TimeType rv = time(0);
	if ( rv == TimeType(-1) )
	    throw runtime_error( string("time() failed") + strerror(errno) );
    #else
       TimeType rv;
       if ( gettimeofday( &rv, 0 ) != 0 )
	    throw runtime_error( string("gettimeofday() failed") + strerror(errno) );
    #endif

    return rv;
    }



//-----------------------------------------------------------------------------
// restart()
//-----------------------------------------------------------------------------

void CodeTimer::restart()
    {
    start = getNow();
    }



//-----------------------------------------------------------------------------
// elapsedSeconds()
//-----------------------------------------------------------------------------

double CodeTimer::elapsedSeconds() const
    {
    return diff( getNow(), start );
    }



//-----------------------------------------------------------------------------
// local_now()
//-----------------------------------------------------------------------------

string CodeTimer::local_now()
    {
    return asctime_without_newline( getNow() );
    }



//-----------------------------------------------------------------------------
// local_elapsed()
//-----------------------------------------------------------------------------

string CodeTimer::local_elapsed() const
    {
    #if USE_TIME_T
	const unsigned long elapsed_ms = (getNow() - start) * 1000;
    #else
	const TimeType now = getNow();
	const unsigned long elapsed_ms = (now.tv_sec-start.tv_sec) * 1000 +
					 (now.tv_usec-start.tv_usec+500)/1000;
    #endif

    stringstream str;

    str << (elapsed_ms/1000)
	<< '.' << setw(3) << setfill('0') << (elapsed_ms%1000) << 's';

    if ( elapsed_ms >= 60000 )
	{
	const unsigned int hours     = elapsed_ms / 3600000;
	const unsigned int mod_hours = elapsed_ms % 3600000;
	const unsigned int minutes   = mod_hours / 60000;
	const unsigned int mod_mins  = mod_hours % 60000;
	const unsigned int seconds   = mod_mins / 1000;
	const unsigned int ms	     = mod_mins % 1000;

	str << " = " << hours
	    << ':' << setw(2) << setfill('0') << minutes
	    << ':' << setw(2) << setfill('0') << seconds
	    << '.' << setw(3) << setfill('0') << ms;
	}

    return str.str();
    }



//-----------------------------------------------------------------------------
// local_started()
//-----------------------------------------------------------------------------

string CodeTimer::local_started () const
    {
    return asctime_without_newline( start );
    }



} // ---- end namespace genepi
