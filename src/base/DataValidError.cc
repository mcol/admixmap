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
/// \file DataValidError.cc
/// Implementation of the genepi::DataValidError class.
//=============================================================================

#include "DataValidError.h"

#include "bclib/estr.h"


using namespace std;



namespace genepi { // ----



//-----------------------------------------------------------------------------
// Class (static) variables:
//-----------------------------------------------------------------------------

bool DataValidError::emacsStyleFmt = true;



//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------

DataValidError::DataValidError( const string & msg, const string & fn, int ln ) :
	message	 ( msg ) ,
	fileName ( fn  ) ,
	lineNum	 ( ln  )
    {
    }



//-----------------------------------------------------------------------------
// Destructor
//-----------------------------------------------------------------------------

DataValidError::~DataValidError() throw()
    {
    }



//-----------------------------------------------------------------------------
/// Get a diagnostic message, ready to deliver to the user (equivalent to
/// what(), but returns a std::string rather than char pointer).
//-----------------------------------------------------------------------------

std::string DataValidError::getFormattedMsg() const
    {
    estr rv;

    if ( emacsStyleFmt )
	(rv = fileName) << ':' << lineNum << ": " << message;
    else
	(rv = message) << " at " << fileName << ':' << lineNum;

    return rv;
    }



//-----------------------------------------------------------------------------
// what() [virtual, overridden from std::exception]
//-----------------------------------------------------------------------------

const char * DataValidError::what() const throw()
    {
    // Yuck!  But what else can one do?
    const int NBUFS = 3;
    static estr bufs[ NBUFS ];
    static int bufN = 0;

    estr & rv = bufs[ bufN ];
    if ( ++bufN == NBUFS )
	bufN = 0;

    rv = getFormattedMsg();

    return rv.c_str();
    }



} // ---- end namespace genepi
