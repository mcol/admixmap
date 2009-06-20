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
/// \file DataValidError.h
//=============================================================================

#ifndef __base_DataValidError_h
#define __base_DataValidError_h



#include <exception>
#include <string>



namespace genepi { // ----



/** \addtogroup base
 * @{ */



/// Exception class for validation errors in user data files
class DataValidError : public std::exception
    {
    private:
	std::string message ;
	std::string fileName;
	int	    lineNum ;

	static bool emacsStyleFmt;

    public:

	static bool getEmacsStyleFmt()		{ return emacsStyleFmt; }
	static void setEmacsStyleFmt( bool nv ) { emacsStyleFmt = nv  ; }

	const std::string & getFileName() const { return fileName; }
	int		    getLineNum () const { return lineNum ; }
	const std::string & getMessage () const { return message ; }

	std::string & getMessageWritable() { return message ; }

	std::string getFormattedMsg() const;

	/// Overridden from std::exception
	virtual const char * what() const throw();

	DataValidError( const std::string & msg, const std::string & fn, int ln );
	~DataValidError() throw();
    };



} // ---- end namespace genepi



/** @} */



#endif // ! __base_DataValidError_h
