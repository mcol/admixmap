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
/// \file Exceptions.h
//=============================================================================

#ifndef __base_Exceptions_h
#define __base_Exceptions_h



#include <stdexcept>
#include "config.h"	// AGGRESSIVE_RANGE_CHECK

#include "estr.h"


#define GP_NDEBUG (! AGGRESSIVE_RANGE_CHECK)



#if __GNUC__
    #define GP_NO_RETURN __attribute__(( noreturn ))
#else
    #define GP_NO_RETURN
#endif



namespace genepi { // ----



/** \addtogroup base
 * @{ */



#if GP_NDEBUG
    #define gp_assert(x)
    #define gp_assert_eq(x,y)
#else

    #define gp_assert(x) \
	do  { \
	    if ( ! (x) ) \
		throw std::logic_error( genepi::estr("Assertion failed at " __FILE__ ":") + __LINE__ + \
			": " #x ); \
	    } while (0)

    #define gp_assert_eq(x,y) \
	do  { \
	    if ( (x) != (y) ) \
		throw std::logic_error( genepi::estr("Assertion failed at " __FILE__ ":") + __LINE__ + \
			": " #x " (" + x + ") == " #y " (" + y + ')' ); \
	    } while (0)

    #define gp_assert_ne(x,y) \
	do  { \
	    if ( (x) == (y) ) \
		throw std::logic_error( genepi::estr("Assertion failed at " __FILE__ ":") + __LINE__ + \
			": " #x " (" + x + ") != " #y " (" + y + ')' ); \
	    } while (0)

    #define gp_assert_gt(x,y) \
	do  { \
	    if ( (x) <= (y) ) \
		throw std::logic_error( genepi::estr("Assertion failed at " __FILE__ ":") + __LINE__ + \
			": " #x " (" + x + ") > " #y " (" + y + ')' ); \
	    } while (0)

    #define gp_assert_ge(x,y) \
	do  { \
	    if ( (x) < (y) ) \
		throw std::logic_error( genepi::estr("Assertion failed at " __FILE__ ":") + __LINE__ + \
			": " #x " (" + x + ") >= " #y " (" + y + ')' ); \
	    } while (0)

    #define gp_assert_lt(x,y) \
	do  { \
	    if ( (x) >= (y) ) \
		throw std::logic_error( genepi::estr("Assertion failed at " __FILE__ ":") + __LINE__ + \
			": " #x " (" + x + ") < " #y " (" + y + ')' ); \
	    } while (0)

    #define gp_assert_le(x,y) \
	do  { \
	    if ( (x) > (y) ) \
		throw std::logic_error( genepi::estr("Assertion failed at " __FILE__ ":") + __LINE__ + \
			": " #x " (" + x + ") <= " #y " (" + y + ')' ); \
	    } while (0)

#endif



void throwSysErr( const char * sc );
#define GP_SAFE_SYS_CALL(sc) if ( (sc) != 0 ) throwSysErr(#sc);



} // end namespace genepi ----



/** @} */



#endif // ! __base_Exceptions_h
