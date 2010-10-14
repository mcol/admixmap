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
/// \file exceptions.h
/// Several macros to throw exceptions, primarily assertions of equality or
/// inequality.  Note that this file was originally in the "base" library, but
/// was moved to bclib because these macros are needed by some of the utility
/// classes that are in bclib.  Because of this, its filename now differs from
/// Exceptions.h (in the same directory) only in case, a condition that is
/// somewhat undesirable and should be corrected.
//=============================================================================

#ifndef __bclib_exceptions_h
#define __bclib_exceptions_h



#include <stdexcept>

#ifdef HAVE_CONFIG_H
#include "config.h"	// AGGRESSIVE_RANGE_CHECK
#endif

#include "estr.h"


#define GP_NDEBUG (! AGGRESSIVE_RANGE_CHECK)



#if __GNUC__
    #define GP_NO_RETURN __attribute__(( noreturn ))
#else
    #define GP_NO_RETURN
#endif



namespace genepi { // ----



/** \addtogroup bclib
 * @{ */



#if GP_NDEBUG
    template <typename T	    > void t_gp_assert( const T & /*x*/ ) {}
    template <typename T, typename U> void t_gp_assert_eq( const T & /*x*/, const U & /*y*/ ) {}
    template <typename T, typename U> void t_gp_assert_ne( const T & /*x*/, const U & /*y*/ ) {}
    template <typename T, typename U> void t_gp_assert_gt( const T & /*x*/, const U & /*y*/ ) {}
    template <typename T, typename U> void t_gp_assert_ge( const T & /*x*/, const U & /*y*/ ) {}
    template <typename T, typename U> void t_gp_assert_lt( const T & /*x*/, const U & /*y*/ ) {}
    template <typename T, typename U> void t_gp_assert_le( const T & /*x*/, const U & /*y*/ ) {}

    #define gp_assert(x)      genepi::t_gp_assert(x)
    #define gp_assert_eq(x,y) genepi::t_gp_assert_eq(x,y)
    #define gp_assert_ne(x,y) genepi::t_gp_assert_ne(x,y)
    #define gp_assert_gt(x,y) genepi::t_gp_assert_gt(x,y)
    #define gp_assert_ge(x,y) genepi::t_gp_assert_ge(x,y)
    #define gp_assert_lt(x,y) genepi::t_gp_assert_lt(x,y)
    #define gp_assert_le(x,y) genepi::t_gp_assert_le(x,y)
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



#endif // ! __bclib_exceptions_h
