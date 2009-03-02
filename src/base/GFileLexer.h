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
/// \file GFileLexer.h
/// Definition of the GFileLexer class
//=============================================================================


#ifndef __base_GFileLexer_h
#define __base_GFileLexer_h


#include <fstream>
#include <list>
#include <stdexcept>
#include <string>
#include <vector>

#include "exceptions.h"		// GP_NO_RETURN
#include "estr.h"
#include "DataValidError.h"
#include "Genotype.h"



namespace genepi { // ----



/** \addtogroup base
 * @{ */



//-----------------------------------------------------------------------------
//
/// Class to lexically analyze various genetics input files.
///
/// Used to read the various input files for <CODE>admixmap</CODE> and
/// <CODE>hapmixmap</CODE>, stripping comments, and breaking the input lines
/// into tokens as follows:
///
/// - Tokens are separated by whitespace (i.e. space or tab)
/// - Newlines are significant (they are lexed as a token, T_EOL).
/// - Other whitespace is folded into a token separator and otherwise ignored.
/// - Each token is of one of the following types: string, integer,
///   floating-point, or genotype.
/// - Integers consist of only digits and an optional preceding minus sign
/// - Floating-point consist of digits, a single decimal point which may appear
///   at any point in the number, including as its first character,
/// - Genotypes are a pair of integers, separated by a delimiter character, typically
///   ','.  Either or both of the integers may be missing, in which case a
///   special "missing" marker is stored in its place in the Genotype object.
///   As a special case, a genotype may be surrounded by double-quotes.  Another
///   special case: a value of 0 for either of the alleles is equivalent to
///   being missing.
/// - Strings are essentially any other token; they may be quoted or unquoted.
///   Strings may be quoted with either single or double quotes, but the closing
///   quote character must match the opening one.  If the same quote-character
///   that delimits the string also appears within it, the character must be
///   escaped by preceding with a backslash (\).  Backslashes must also be
///   escaped, and no other characters may be escaped.  Strings that contain
///   whitespace must be quoted; unquoted strings end at the first whitespace,
///   newline, or end-of-file that is encountered.
/// - Comments which begin with '#' and end at the end-of-line are ignored (the
///   newline is not considered part of the comment and is returned as a T_EOL
///   token; thus comments that start at the beginning of a line appear as blank
///   lines.  It is up to the upper-level parser how to interpret blank lines.
///   Comments may begin anywhere on a line, but a comment character inside a
///   string (quoted or not) does not begin a comment.
///
/// The main operation of the class is via lexToken(), which is called
/// repeatedly by the parser, and returns a Token object of arbitrary type.  If
/// the expected token-type is known in advance, lexInteger(), lexFloat() and
/// lexGenotype() are provided for convenience, and an entire line can be parsed
/// via lexLine() or lexLineStrings().  A token pushback-stack is provided to
/// the upper-level parser via pushback(const Token &).
///
/// This class does not have any concept of the semantic meaning of the tokens;
/// it can be used either by deriving a class which is the upper-level parser
/// and calls the lexing methods on itself, or by being instantiated and used
/// (called) from a separate class.
///
/// Currently there are no hooks or virtual methods to easily allow derived
/// classes to easily define new token-types, but it can be done by
/// re-interpreting string tokens, or by adding virtual methods here.
//
//-----------------------------------------------------------------------------

class GFileLexer
    {
    public:

	typedef std::string string; ///< convenience

	typedef Genotype GType;


	//---------------------------------------------------------------------
	/// Exception class for syntax errors in input files.
	///
	/// The lexical-analyzer will throw one of these when it encounters a
	/// syntax-error in an input file; upper-level parsers are welcome to do
	/// the same, most easily done by calling throwError().  These should
	/// be caught at some higher level (or in the global uncaught-exception
	/// handler), an appropriate error-message printed (presumably just by
	/// calling what()), and appropriate action taken (typically abort
	/// execution, or a more robust handler).
	//---------------------------------------------------------------------

	class ParseError : public DataValidError
	    {
	    public:
		ParseError( const string & msg, const string & fn, int ln ) :
		    DataValidError( msg, fn, ln ) {}
	    };



	//---------------------------------------------------------------------
	/// Enumeration of the possible types of tokens.  Derived classes that
	/// wish to define additional token-types should begin at N_TOKEN_TYPES.
	//---------------------------------------------------------------------

	enum TokenType
	    {
	    T_STRING	,
	    T_INTEGER	,
	    T_FLOAT	,
	    T_GTYPE	,  /// Pair of integers (unsigned short)
	    T_EOL	,  /// End-of-line
	    T_EOF	,  /// End-of-file
	    N_TOKEN_TYPES  /// Always keep last
	    };



	//---------------------------------------------------------------------
	/// Class to hold one token: essentially a type tag (TokenType) and a
	/// union of the possible values (which can only be externally accessed
	/// via type-safe methods, of course).
	//---------------------------------------------------------------------
	class Token
	    {
	    friend class GFileLexer;
	    private:
		TokenType  type	    ;
		string	   strVal   ;
		union
		    {
		    long   intVal   ;
		    double floatVal ;
		    GType  gtypeVal ; // See NOTE *3* in Genotype.h
		    };

		/// Throw an exception indicating an unexpected token-type
		void throwExpectType( TokenType t ) const GP_NO_RETURN;

	    public:
		TokenType    getType	() const { return type; }
		const char * getTypeName() const;

		bool isType( TokenType t ) const { return (type == t); }

		void assert_type( TokenType t ) const
		    { if ( ! isType(t) ) throwExpectType( t ); }

		const string &	getStrVal   () const { assert_type( T_STRING  ); return strVal  ; }
		long		getIntVal   () const { assert_type( T_INTEGER ); return intVal  ; }
		double		getFloatVal () const { assert_type( T_FLOAT   ); return floatVal; }
		GType		getGTypeVal () const { assert_type( T_GTYPE   ); return gtypeVal; }

		/// Is this a semantically significant token (i.e. <I>not</I>
		/// the "separator" end-of-line or end-of-file tokens).
		bool isToken() const { return (type != T_EOL) && (type != T_EOF); }

		/// Create a string representation of the token.  This is not
		/// guaranteed to be the token exactly as it appeared in the
		/// input file, but often will be.  In the case of string
		/// tokens, it is the string itself (no conversion required); in
		/// other cases, the value is converted to string, a relatively
		/// inefficient operation, so use with caution.
		estr asString() const;

		/// Value as a floating-point.  All of the following will produce 1.0d:
		/// 1 1.0 "1" "1.0"
		/// Except that the strings are not yet implemented.
		double asFloat() const;
	    };



    private:

	// Private data:
	const char *	 fileName	;
	std::ifstream	 istr		;
	int		 lineNum	;
	std::list<char>  pb_cstack	;
	std::list<Token> pb_tstack	;
	char		 gtypeDelim	;

	// Token-parsing private helper methods:
	void doParseGType( long val1	, Token & tok );
	void doParseFloat( long intPart , Token & tok );
	void doParseInt	 ( char ch	, Token & tok );
	void doParseStr	 ( char ch	, Token & tok );


    protected:

	//---------------------------------------------------------------------
	// Methods used by the token-parsing methods:
	//---------------------------------------------------------------------

	void setVals( Genotype & gt, long nv1, long nv2 ) const;
	void setHaploid( Genotype & gt, long nv1 ) const;

	void checkForComments( char & ch );

	char parseOneInt( char ch, long & val, bool & endOfToken, int & nDigs );
	long parseFinalInt( char ch );

	/// Returns the last pushback()ed character; if none, gets a character
	/// from the input stream; if end-of-file reached, returns EOF_MARKER.
	char getchar();

	/// Returned by getchar() when EOF is encountered
	static const char EOF_MARKER = '\0';

	/// Push a character back onto the input stream.
	void pushback( char ch );

	/// Throw a syntax-error exception.  This should probably be public if
	/// the upper-level parsing classes instantiate and use this class
	/// rather than deriving from it.  The error message will contain the
	/// current filename and line number.
	void throwError( const string & msg ) const throw(ParseError) GP_NO_RETURN;
	void assert_type( TokenType t, const Token & tok, const char * fieldName = 0 ) const throw(ParseError);


    public:


	//---------------------------------------------------------------------
	// Constructor/destructor:
	//---------------------------------------------------------------------

	/// If the upper-level parsers are derived classes (rather than
	/// instantiating and using this class), this should probably be
	/// protected.
	GFileLexer( const char * fn );
	~GFileLexer();



	int	     curLineNum () const { return lineNum ; }
	const char * getFileName() const { return fileName; }



	//---------------------------------------------------------------------
	// Get/set of attributes:
	//---------------------------------------------------------------------

	char getGTypeDelim(	    ) const { return gtypeDelim; }
	void setGTypeDelim( char nv )	    { gtypeDelim = nv;	 }


	//---------------------------------------------------------------------
	// Public lexing methods:
	//---------------------------------------------------------------------

	/// Lex the next token of arbitrary type and return it: this is the main
	/// lexing method which all others call.  Any pushed-back tokens will be
	/// returned prior to lexing new ones.
	Token  lexToken();
	void   pushback( const Token & tok );

	/// Skips past blank lines and comments to the next "real" token.
	/// Returns true if token is found, false if EOF before token.
	bool   skipToToken();

	string lexString  (); ///< Accepts an integer, which will be interpreted as a string

	/// Throws an exception if the token is not an integer
	long   lexInteger ( const char * fieldName = 0 );

	double	 lexFloat (); ///< Throws an exception if the next token is not a floating-point
	Genotype lexGType (); ///< Throws an exception if the next token is not a genotype

	void   lexLine	     ( std::vector<Token > & tokens  );
	void   lexLineStrings( std::vector<string> & strings );
    };



} // ---- end namespace genepi



/** @} */



#endif // ! __base_GFileLexer_h
