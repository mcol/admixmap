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
/// \file estr.h
//=============================================================================

#ifndef __base_estr_h
#define __base_estr_h


#include <string>



namespace genepi { // ----



/** \addtogroup base
 * @{ */



//-----------------------------------------------------------------------------
/// A subclass of <A
/// HREF="http://www.icce.rug.nl/documents/cplusplus/cplusplus04.html"><CODE>std::string</CODE></A>
/// that attempts to be slightly more convenient to use.
///
/// In particular,
///	- Integral types may be used as initializers to constructors, and as rhs
///	  to binary operators such as <B>+</B> or <B>+=</B>.
///	- operator<<() is overloaded to behave as one might expect: equivalent
///	  to <B>+=</B> but with left grouping and identical syntax to
///	  <A HREF="http://www.icce.rug.nl/documents/cplusplus/cplusplus05.html"><CODE>ostream</CODE>s</A>.
///
/// More generally, see <B>std::stringstream</B>.
//-----------------------------------------------------------------------------

class estr : public std::string
    {

    public:


	//------------------------------------------------------------------
	// Constructors (specializations and catch-all template):
	//------------------------------------------------------------------
	estr( short	     x );
	estr( unsigned short x );
	estr( int	     x );
	estr( unsigned int   x );
	estr( long	     x );
	estr( unsigned long  x );
	estr( double	     x );

	estr() {}

	template< typename T > estr( const T & x ) : std::string( x ) {}



	//------------------------------------------------------------------
	// Operators:
	//------------------------------------------------------------------


	/// Template to generate any arbitrary <B>operator=()</B>.
	template< typename T > estr & operator=( const T & x )
	    { return reinterpret_cast<estr&>(std::string::operator=(x)); }
	#define ESTR_EQ(T) estr & operator=( T x ) \
	    { return reinterpret_cast<estr&>(std::string::operator=(estr(x))); }
		ESTR_EQ(short);
		ESTR_EQ(unsigned short);
		ESTR_EQ(int);
		ESTR_EQ(unsigned int);
		ESTR_EQ(long);
		ESTR_EQ(unsigned long);
		ESTR_EQ(double);
	#undef ESTR_EQ


	/// Template to generate any arbitrary <B>operator+=()</B>.
	template< typename T > estr & operator+=( const T & x )
	    { return reinterpret_cast<estr&>(std::string::operator+=(x)); }
	#define ESTR_PLEQ(T) estr & operator+=( T x ) \
	    { return reinterpret_cast<estr&>(std::string::operator+=(estr(x))); }
		ESTR_PLEQ(short);
		ESTR_PLEQ(unsigned short);
		ESTR_PLEQ(int);
		ESTR_PLEQ(unsigned int);
		ESTR_PLEQ(long);
		ESTR_PLEQ(unsigned long);
		ESTR_PLEQ(double);
	#undef ESTR_PLEQ


	/// Equivalent to operator+=(), but groups to the left and is
	/// syntactically compatible with ostreams.  This is the template
	/// version, which is a catch-all for all types for which we have no
	/// specific specialization.
	template< typename T > estr & operator<<( const T & x )
	    { return operator+=( x ); }


	/// Template to generate any arbitrary <B>operator+()</B> for which we
	/// haven't a specialization.
	///
	/// Yikes:
	/// -# We hope that this neither creates unnecessary or excessive
	///    copying of strings, nor accesses references to destroyed
	///    temporaries.
	/// -# We hope both that the specializations below are used in
	///    preference to the template, and that the template is used
	///    whenever no specialization exists.
	template< typename T > estr operator+ ( const T & x )
	    { return ((*reinterpret_cast<std::string*>(this)) + x); }

	#define ESTR_PLUS(T) estr operator+( T x ) { return (*this) + estr(x); }
	    ESTR_PLUS(short);
	    ESTR_PLUS(unsigned short);
	    ESTR_PLUS(int);
	    ESTR_PLUS(unsigned int);
	    ESTR_PLUS(long);
	    ESTR_PLUS(unsigned long);
	    ESTR_PLUS(double);
	#undef ESTR_PLUS


	bool equalsCaseInsens( const estr & rhs ) const;
    };


bool equalsCaseInsens( const std::string & lhs, const std::string & rhs );



} // ---- end namespace genepi



/** @} */



#endif // ! __base_estr_h
