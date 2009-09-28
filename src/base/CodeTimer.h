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
/// \file CodeTimer.h
/// Interface of the CodeTimer class.
//=============================================================================

#ifndef __base_CodeTimer_h
#define __base_CodeTimer_h



#define USE_TIME_T	0 // defined(__WIN32)

#if ! USE_TIME_T
    #include <sys/time.h>
#endif

#include <ctime>
#include <string>



namespace genepi { // ----




//-----------------------------------------------------------------------------
/// A very simple timer class to report how long has elapsed since it was constructed.
//-----------------------------------------------------------------------------
class CodeTimer
    {
    private:
	#if USE_TIME_T
	    typedef time_t	   TimeType;
	#else
	    typedef struct timeval TimeType;
	#endif

	TimeType start;

    protected:
	static TimeType getNow();
	static std::string asctime_without_newline( const TimeType & tm );


    public:
	CodeTimer() { restart(); }

	void restart();

	double elapsedSeconds() const;

	static std::string local_now	 ();
	std::string	   local_elapsed () const;
	std::string	   local_started () const;
    };




} // ---- end namespace genepi



#endif // ! __base_CodeTimer_h
