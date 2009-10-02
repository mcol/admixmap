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
/// \file GFileLexer.cc
/// Implementation of the GFileLexer class.
//=============================================================================


#include "GFileLexer.h"


#include <cctype>   // isblank(), isdigit()
#include <iostream> // cerr (for warn())

#include "bclib/estr.h"


#define SUPPORT_QUOTED_GENOTYPE		1
#define SUPPORT_QUOTED_INTEGER		1
#define SUPPORT_SINGLE_HAPLOID		1
#define ACCEPT_INT_AS_STRING		1
#define ALLOW_EXPLICIT_MISSING_VALS	1


#if SUPPORT_QUOTED_INTEGER
    #include <cstdlib>	// strtol()
    #include <cerrno>
#endif


static const char ESCAPE_CHAR	= '\\'	;
static const char COMMENT_CHAR	= '#'	;
static const char ALT_GT_DELIM	= '/'	;



/// Yet more insanity: if turned off, we don't make any assumptions about the
/// lexical type of tokens at this level, but let the upper-level parser decide
/// how to interpret things based on the semantics; i.e. more-or-less everything
/// is a string at this level.  The biggest effect of turning it off is that
/// "1X" (without quotes) will parse as a string rather than throwing a
/// BadIntErr exception.
#define CONTEXT_FREE_SYNTAX	0



using namespace std;



namespace genepi { // ----



//-----------------------------------------------------------------------------
// static helper functions:
//-----------------------------------------------------------------------------

#if defined(__GNUC__) && ((_XOPEN_SOURCE >= 600) || _ISOC99_SOURCE)
    #define isbl(x) isblank(x)
#else
    static bool isbl( char ch ) { return (ch == ' ') || (ch == '\t'); }
#endif


static inline int dval( char ch ) { return (ch - '0'); }


static long long_pow10( int exponent )
    {
    static const long p10tab [] =
	{
	1		,
	10		,
	100		,
	1000		,
	10000		,
	100000		,
	1000000		,
	10000000	,
	100000000	,
	1000000000	,
      #if __LONG_MAX__ > 2147483647L // 64-bit only:
	10000000000	,
	100000000000	,
	1000000000000	,
	10000000000000	,
	100000000000000	,
      #endif
	};
    const int N_ENT = sizeof(p10tab) / sizeof(*p10tab);

    if ( (exponent < 0) || (exponent >= N_ENT) )
	throw std::overflow_error( estr("decimal part has too many digits: ") + exponent );

    return p10tab[ exponent ];
    }



static const char * typeName( GFileLexer::TokenType t )
    {
    static const char * const NAMES [] =
	{
	"string"	,
	"integer"	,
	"float"		,
	"genotype"	,
	"eol"		,
	"eof"		,
	};

    gp_assert( (sizeof(NAMES)/sizeof(*NAMES)) == GFileLexer::N_TOKEN_TYPES );
    gp_assert( t >= 0 );
    gp_assert( t < GFileLexer::N_TOKEN_TYPES );

    return NAMES[ t ];
    }



//-----------------------------------------------------------------------------
// BadIntErr()
//-----------------------------------------------------------------------------

GFileLexer::BadIntErr::BadIntErr( const string & msg, const string & fn, int ln ) :
    ParseError( "badly formed integer: " + msg, fn, ln ) {}
GFileLexer::BadIntErr::BadIntErr( char invalidChar, const string & fn, int ln ) :
    ParseError( std::string("invalid character `") + invalidChar + "' in integer.", fn, ln ) {}



//-----------------------------------------------------------------------------
//  Token::asString()
/// Not terribly efficient, avoid when possible
//-----------------------------------------------------------------------------

estr GFileLexer::Token::asString() const
    {
    if ( type == T_STRING )
	return strVal;
    else if ( type == T_INTEGER )
	return estr( intVal );
    else if ( type == T_FLOAT )
	return estr( floatVal );
    else if ( type == T_GTYPE )
	return gtypeVal.desc();
    else if ( type == T_EOL )
	return "<end-of-line>";
    else if ( type == T_EOF )
	return "<end-of-file>";
    else
	throw logic_error( estr("Unknown token type ") + int(type) );
    }



//-----------------------------------------------------------------------------
//
//  Token::asFloat()
//
/// Supposed to support strings-as-floats ("1.2"), but currently does not.
///
//
//-----------------------------------------------------------------------------

double GFileLexer::Token::asFloat( const char * fieldName ) const
    {
    if ( isType( T_FLOAT ) )
	return floatVal;
    else if ( isType( T_INTEGER ) )
	return intVal;
    else
	{
	estr msg( "Can't convert " );
	if ( fieldName != 0 )
	    msg << fieldName << " from type ";
	msg << getTypeName() << " to floating-point";
	throw logic_error( msg );
	}
    }



//-----------------------------------------------------------------------------
// Token::throwExpectType() [protected]
//-----------------------------------------------------------------------------

void GFileLexer::Token::throwExpectType( TokenType t ) const
    {
    throw runtime_error( string("expected ") + typeName(t) + " but found " +
			    getTypeName() + " (" + asString() + ')' );
    }



//-----------------------------------------------------------------------------
// Token::getTypeName()
//-----------------------------------------------------------------------------

const char * GFileLexer::Token::getTypeName() const
    {
    return typeName( type );
    }



//-----------------------------------------------------------------------------
// getchar() [protected]
//-----------------------------------------------------------------------------

char GFileLexer::getchar()
    {
    char rv;

    if ( pb_cstack.empty() )
	{
	if ( ! istr.get( rv ) )
	    rv = EOF_MARKER;
	else
	    {
	    // Get a new character from the input stream, translating CR-LF into
	    // a single LF:
	    if ( rv == '\r' )
		{
		if ( ! istr.get( rv ) )
		    rv = EOF_MARKER;
		if ( rv != '\n' )
		    {
		    pushback( rv );
		    rv = '\r';
		    }
		}
	    }
	}
    else
	{
	rv = pb_cstack.back();
	pb_cstack.pop_back();
	}

    if ( rv == '\n' )
	++lineNum;

    return rv;
    }



//-----------------------------------------------------------------------------
// pushback(char) [protected]
//
/// In most (all?) cases, this can be implemented via the token pushback stack.
//-----------------------------------------------------------------------------

void GFileLexer::pushback( char ch )
    {
    if ( ch == '\n' )
	--lineNum;

    pb_cstack.push_back( ch );
    }



//-----------------------------------------------------------------------------
// throwError() [protected]
//-----------------------------------------------------------------------------

void GFileLexer::throwError( const string & msg, bool atLastTokenLine ) const
	throw( ParseError )
    {
    throw ParseError( msg, fileName, atLastTokenLine ? getLastTokenLine() : getCurLineNum() );
    }



//-----------------------------------------------------------------------------
// assert_type() [protected]
//-----------------------------------------------------------------------------

void GFileLexer::assert_type( TokenType t, const Token & tok, const char * fieldName ) const
	throw( ParseError )
    {
    if ( tok.getType() != t )
	{
	estr msg;
	if ( fieldName != 0 )
	    {
	    msg = "in field \"";
	    msg << fieldName << "\": ";
	    }
	msg << "expected " << typeName(t) << " but found " << tok.getTypeName()
	    << " (" << tok.asString() << ')';
	throwError( msg );
	}
    }



//-----------------------------------------------------------------------------
// warn() [protected]
//-----------------------------------------------------------------------------

void GFileLexer::warn( const string & msg ) const
    {
    if ( nWarnings == 0 )
	std::cerr << '\n';
    std::cerr << getFileName() << ':' << getLastTokenLine() << ": warning: " << msg << '\n';
    ++nWarnings;
    }



//-----------------------------------------------------------------------------
// checkForComments() [protected]
/// ch, passed in by reference, is a freshly getchar()d character; we check here
/// for a comment character, and if found, skip until the end-of-line, in which
/// case the terminal '\n' is left in ch.
//-----------------------------------------------------------------------------

void GFileLexer::checkForComments( char & ch )
    {
    if ( ch == COMMENT_CHAR )
	do
	    {
	    ch = getchar();
	    } while ( (ch != '\n') && (ch != EOF_MARKER) );
    }



//-----------------------------------------------------------------------------
// parseOneInt() [protected]
//
/// Parse an integer from the input stream, stopping at the first non-integer
/// character that is encountered, which is returned as the return value (and
/// pushed back onto the input stream if it is a newline so that it can later be
/// relexed as a T_EOL).  The result is returned in val.  The first character
/// is assumed to have been already read from the input stream, and is passed in
/// as ch.  endOfToken will be set if the terminating character represents the
/// end of the token (i.e. whitespace or end-of-file).  The number of digits in
/// the integer is returned in nDigs.
//-----------------------------------------------------------------------------

char GFileLexer::parseOneInt( char ch, long & val, bool & endOfToken, int & nDigs )
    {
    long accum;
    int  aNDigs;

    const bool isNegative = (ch == '-');
    if ( isNegative )
	{
	accum = 0;
	aNDigs = 0;
	}
    else
	{
	if ( ! isdigit(ch) )
	    throw BadIntErr( ch, fileName, getLastTokenLine() );
	accum = dval( ch );
	aNDigs = 1;
	}


    while ( isdigit( ch = getchar() ) )
	{
	accum = (accum * 10) + dval( ch ); // TODO: check for integer overflow
	++aNDigs;
	}

    // Skip comments:
    checkForComments( ch );

    endOfToken = (ch == '\n') || (ch == EOF_MARKER) || isbl( ch );

    if ( ch == '\n' )
	pushback( ch );

    val = isNegative ? -accum : accum;
    nDigs = aNDigs;

    return ch;
    }



//-----------------------------------------------------------------------------
//  parseFinalInt() [protected]
//
/// Lex an integer which cannot be part of a larger token (must be terminated by
/// whitespace, EOL, EOF, etc.)
//-----------------------------------------------------------------------------

long GFileLexer::parseFinalInt( char ch )
    {
    bool endOfToken;
    long rv;
    int  nDigs = 0; // Initialize to suppress compiler warning

    ch = parseOneInt( ch, rv, endOfToken, nDigs );

    if ( ! endOfToken )
	throw BadIntErr( ch, fileName, getLastTokenLine() );

    return rv;
    }



//-----------------------------------------------------------------------------
// setVals() [protected]
//-----------------------------------------------------------------------------

inline void GFileLexer::setVals( Genotype & gt, long nv1, long nv2 ) const
    {
    if ( nv1 == 0 ) // 0 is an alias for missing
	gt.val1 = GType::MISSING_VAL;
    else
	{
	gt.val1 = static_cast<unsigned short>( nv1 );
	if ( nv1 != static_cast<long>( gt.val1 ) )
	    throwError( "overflow: first allele of genotype exceeds unsigned short" );
	}

    if ( nv2 == 0 ) // 0 is an alias for missing
	gt.val2 = GType::MISSING_VAL;
    else
	{
	gt.val2 = static_cast<unsigned short>( nv2 );
	if ( nv2 != static_cast<long>( gt.val2 ) )
	    throwError( "overflow: second allele of genotype exceeds unsigned short" );
	}

    // Since the pair of integers are unordered, we will assure that val2 is
    // missing before val1 is.
    gt.condense();
    }



//-----------------------------------------------------------------------------
// GType::setHaploid() [protected]
//-----------------------------------------------------------------------------

inline void GFileLexer::setHaploid( Genotype & gt, long nv1 ) const
    {
    setVals( gt, nv1, Genotype::HAPLOID_VAL );
    }



//-----------------------------------------------------------------------------
// doParseGType() [private]
//
/// Is called just after the delimiter is parsed; val1 contains the first integer.
//-----------------------------------------------------------------------------

void GFileLexer::doParseGType( long val1, Token & tok )
    {
    char ch = getchar();

    if ( isdigit(ch) )
	{
	const long val2 = parseFinalInt( ch );
	#if ! ALLOW_EXPLICIT_MISSING_VALS
	    if ( val2 == GType::MISSING_VAL )
		throwError( "second allele uses the missing-flag value" );
	#endif
	setVals( tok.gtypeVal, val1, val2 );
	}
    else
	{
	checkForComments( ch );
	if ( ch == '\n' )
	    pushback( ch );
	else if ( ! isbl( ch ) )
	    throwError( "non-numeric as second component of gtype" );
	setVals( tok.gtypeVal, val1, GType::MISSING_VAL );
	}

    tok.type = T_GTYPE;
    }



//-----------------------------------------------------------------------------
// doParseFloat() [private]
//
/// Is called just after the decimal point is parsed; intPart contains the
/// integer part parsed prior to the decimal point.
//-----------------------------------------------------------------------------

void GFileLexer::doParseFloat( long intPart, Token & tok )
    {
    long decPart;
    int	 nDecDigs = 0; // Initialize to suppress compiler warning
    bool endOfToken;

    const char ch = parseOneInt( getchar(), decPart, endOfToken, nDecDigs );

    if ( ! endOfToken )
	throwError( string("invalid character '") + ch + "' in floating-point" );

    // This is the result of something like: 21.-03
    if ( decPart < 0 )
	throwError( "invalid character '-' after decimal point" );

    tok.type = T_FLOAT;
    tok.floatVal = double(intPart) + (double(decPart) / long_pow10( nDecDigs ));
    }



//-----------------------------------------------------------------------------
// doParseInt() [private]
//
/// NOTE! If CONTEXT_FREE_SYNTAX is not turned on, it is possible to return a
///	  token of type string here!!!
//-----------------------------------------------------------------------------

void GFileLexer::doParseInt( char ch, Token & tok )
    {
    bool endOfToken;
    long val;
    int  nDigs = 0; // Initialize to suppress compiler warning

    ch = parseOneInt( ch, val, endOfToken, nDigs );

    if ( (ch == gtypeDelim) || (ch == ALT_GT_DELIM) )
	{
	#if ! ALLOW_EXPLICIT_MISSING_VALS
	    if ( val == GType::MISSING_VAL )
		throwError( "first allele uses the missing-flag value" );
	#endif
	doParseGType( val, tok );
	}
    else if ( ch == '.' )
	doParseFloat( val, tok );
    else if ( endOfToken )
	{
	tok.type = T_INTEGER;
	tok.intVal = val;
	}
    else
	{
	// Part-way through an integer, we encounter a non-digit character:
	#if CONTEXT_FREE_SYNTAX
	    throw BadIntErr( ch, fileName, getLastTokenLine() );
	#else
	    doParseStr( ch, tok );
	    tok.strVal.insert( 0, estr(val) ); // .prepend()
	#endif
	}
    }



//-----------------------------------------------------------------------------
// doParseStr() [private]
//-----------------------------------------------------------------------------

void GFileLexer::doParseStr( char ch, Token & tok )
    {
    const char NO_QUOTE = '\0';
    char quoteChar;

    if ( (ch == '\'') || (ch == '"') )
	{
	quoteChar = ch;
	ch = getchar();
	}
    else
	quoteChar = NO_QUOTE;

    tok.strVal.clear();

    while ( (ch != '\n') && (ch != EOF_MARKER) &&
	    ((quoteChar == NO_QUOTE) ? (! isbl(ch)) : (ch != quoteChar)) )
	{
	if ( ch == ESCAPE_CHAR )
	    {
	    ch = getchar();
	    if ( (ch != ESCAPE_CHAR) && (ch != quoteChar) )
		throwError( "only quote-characters and escape-characters "
				"may be escaped in quoted strings" );
	    }
	tok.strVal += ch;
	ch = getchar();
	}

    // Disallow EOL and EOF inside quoted strings:
    if ( quoteChar != NO_QUOTE )
	{
	if ( ch == '\n' )
	    throwError( "end-of-line inside quoted string" );
	else if ( ch == EOF_MARKER )
	    throwError( "end-of-file inside quoted string" );
	}

    if ( ch == '\n' )
	pushback( ch );

    tok.type = T_STRING;
    }



//-----------------------------------------------------------------------------
// lexToken()
//-----------------------------------------------------------------------------

GFileLexer::Token GFileLexer::lexToken()
    {
    Token rv;

    if ( ! pb_tstack.empty() )
	{
	rv = pb_tstack.back();
	pb_tstack.pop_back();
	}

    else
	{

	char ch = getchar();

	// Skip leading whitespace:
	while ( isblank(ch) )
	    ch = getchar();

	// Skip comments:
	checkForComments( ch );

	rv.lineNum = getCurLineNum();

	if	( ch == '\n'	    ) rv.type = T_EOL;
	else if ( ch == EOF_MARKER  ) rv.type = T_EOF;
	else if ( ch == gtypeDelim  ) doParseGType( GType::MISSING_VAL, rv );
	else if ( ch == ALT_GT_DELIM) doParseGType( GType::MISSING_VAL, rv );
	else if ( ch == '.'	    ) doParseFloat( 0, rv );
	else if ( isdigit( ch )	    ) doParseInt( ch, rv );
	else if ( ch == '-'	    ) doParseInt( ch, rv );
	else			      doParseStr( ch, rv );

	}

    lastTokenLine = rv.getLineNum();

    return rv;
    }



//-----------------------------------------------------------------------------
// pushback(Token)
//-----------------------------------------------------------------------------

void GFileLexer::pushback( const Token & tok )
    {
    // It's not clear that this is desirable; perhaps better to just leave
    // lastTokenLine alone:
    if ( ! pb_tstack.empty() )
	lastTokenLine = pb_tstack.back().getLineNum();

    pb_tstack.push_back( tok );
    }



//-----------------------------------------------------------------------------
// skipToToken()
//-----------------------------------------------------------------------------

bool GFileLexer::skipToToken()
    {
    Token tok;

    for ( tok = lexToken(); ! tok.isToken(); tok = lexToken() )
	if ( tok.isType( T_EOF ) )
	    return false; // **** RETURN HERE ****

    pushback( tok );

    return true;
    }



//-----------------------------------------------------------------------------
// lexString()
//-----------------------------------------------------------------------------

string GFileLexer::lexString()
    {

    #if ACCEPT_INT_AS_STRING

	Token tok = lexToken();
	if ( tok.isType( T_INTEGER ) )
	    {
	    tok.strVal = tok.asString();
	    tok.type = T_STRING;
	    }
	else
	    assert_type( T_STRING, tok );

    #else

	const Token & tok = lexToken();
	assert_type( T_STRING, tok );

    #endif

    return tok.strVal;
    }



//-----------------------------------------------------------------------------
// lexInteger()
//-----------------------------------------------------------------------------

long GFileLexer::lexInteger( const char * fieldName )
    {
    const Token & tok = lexToken();

    // Inexplicably, we allow this:
    #if SUPPORT_QUOTED_INTEGER
	if ( tok.isType( T_STRING ) )
	    {
	    char * endPtr;
	    errno = 0;
	    const char * tPtr = tok.strVal.c_str();

	    // Skip leading whitespace:
	    while ( isspace(*tPtr) )
		++tPtr;

	    const long val = strtol( tok.strVal.c_str(), &endPtr, 10 );
	    if ( errno == ERANGE )
		{
		if ( val == LONG_MAX )
		    throwError( tok.strVal + ": integer overflow (too large)" );
		if ( val == LONG_MIN )
		    throwError( tok.strVal + ": integer underflow (too large negative)" );
		gp_assert( false ); // WTF, undocumented here
		}
	    gp_assert_eq( errno, 0 );

	    // Skip trailing whitespace:
	    while ( (endPtr != 0) && isspace(*endPtr) )
		++endPtr;

	    if ( (endPtr != 0) &&    // no digits at all
		 (*endPtr == '\0') ) // had trailing non-digits
		return val;
	    }
    #endif

    assert_type( T_INTEGER, tok, fieldName );
    return tok.intVal;
    }



//-----------------------------------------------------------------------------
// lexFloat()
//-----------------------------------------------------------------------------

double GFileLexer::lexFloat( const char * fieldName )
    {
    const Token & tok = lexToken();

    // Accept an integer as a float:
    if ( tok.isType( T_INTEGER ) )
	return tok.intVal;

    assert_type( T_FLOAT, tok, fieldName );
    return tok.floatVal;
    }



//-----------------------------------------------------------------------------
// lexGType()
//-----------------------------------------------------------------------------

Genotype GFileLexer::lexGType( const char * fieldName )
    {
    Token tok = lexToken();


    // Special-case: single integer is a haploid:
    #if SUPPORT_SINGLE_HAPLOID
	if ( tok.isType( T_INTEGER ) )
	    {
	    tok.type = T_GTYPE;
	    setHaploid( tok.gtypeVal, tok.intVal );
	    }
    #endif


    // Special-cast: allow genotype data (pairs of integers) to be quoted:
    #if SUPPORT_QUOTED_GENOTYPE
	if ( tok.isType( T_STRING ) )
	    {
	    if ( (tok.strVal.size() == 3) &&
		    isdigit(tok.strVal[0]) &&
		    ((tok.strVal[1] == gtypeDelim) || (tok.strVal[1] == ALT_GT_DELIM)) &&
		    isdigit(tok.strVal[2]) )
		{
		tok.type = T_GTYPE;
		setVals( tok.gtypeVal, tok.strVal[0] - '0', tok.strVal[2] - '0' );
		}
	    }
    #endif


    assert_type( T_GTYPE, tok, fieldName );
    return tok.gtypeVal;
    }



//-----------------------------------------------------------------------------
// lexLine()
//-----------------------------------------------------------------------------

void GFileLexer::lexLine( std::vector<Token> & tokens )
    {
    for ( Token tok = lexToken() ; tok.isToken() ; tok = lexToken() )
	tokens.push_back( tok );
    }



//-----------------------------------------------------------------------------
// lexLineStrings()
//-----------------------------------------------------------------------------

void GFileLexer::lexLineStrings( std::vector<string> & strings )
    {
    for ( Token tok = lexToken() ; tok.isToken() ; tok = lexToken() )
	strings.push_back( tok.asString() );
    }



//-----------------------------------------------------------------------------
// Constructor
/// \param fn The filename to open for reading
//-----------------------------------------------------------------------------

GFileLexer::GFileLexer( const char * fn ) :
	fileName	( fn	) ,
	istr		( fn	) ,
	lineNum		( 1	) ,
	lastTokenLine	( 0	) ,
	nWarnings	( 0	) ,
	gtypeDelim	( ','	)   // Default delimiter
    {
    if ( istr.bad() || (! istr.is_open()) )
	{
	string msg( "opening file(\"" );
	msg += fn;
	msg += "\")";
	throwSysErr( msg.c_str() );
	}

    // Set istr to throw an exception on errors (this does not include attempt
    // to read past end-of-file):
    istr.exceptions( std::ios_base::badbit );
    }



//-----------------------------------------------------------------------------
// Destructor
//-----------------------------------------------------------------------------

GFileLexer::~GFileLexer()
    {
    }




} // ---- end namespace genepi
