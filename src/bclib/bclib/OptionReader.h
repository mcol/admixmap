//=============================================================================
//
// Copyright (C) 2007  David O'Donnell
// Portions Copyright (C) 2009  David D. Favro
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
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
/// \file OptionReader.h
/// Definition of the bclib::OptionReader class.
//=============================================================================
 
/*
 *   Acts as a utility class for reading
 *   commandline options as well as options from
 *   an optionfile with delimited type value pairs.
 *   Can be extended (used as a base class) to provide more customised option handling, if required.
 *
 *   Handles the POSIX style single character options ( -w )
 *   as well as the newer GNU long options ( --width )
 *   the '-' signs are not mandatory, however. In particular, this is not intended 
 *   to be used in programs that take an argument, unless the arguments are recast as (mandatory) options.
 *
 *  The option file assumes the traditional format of
 *  first character based comment lines and type value
 *  pairs with a delimiter , and flags which are not pairs (not yet implemented)
 *  The delimiter character can be changed.
 *
 *  	# this is a coment
 *  	# next line is an option value pair
 *  	width = 100
 *     	# next line is a flag
 *      noimages
 *
 * Does not support printing out Help and Usage since it is just as easy to write your own as supplying many long strings.
 *
 * Why not just use getopt() ?
 *
 *   getopt is C, this is pure C++ [DDF: actually, non-type-safe code like this
 *   is not very "pure" C++] and a lot easier to use Also, getopt is a POSIX
 *   standard not part of ANSI-C.  So it may not be available on platforms like
 *   Windows.
 */

#ifndef __bclib_OptionReader_h
#define __bclib_OptionReader_h

#include "bclib/bclib.h"
#include "bclib/cvector.h"
#include <map>
#include <string>
#include <vector>

using namespace::std;

BEGIN_BCLIB_NAMESPACE


/** \addtogroup bclib
 * @{ */


/// The type of the parameter of an option
enum OptionType{ nullOption       ,
                 boolOption       ,
                 intOption        ,
                 longOption       ,
                 floatOption      ,
                 doubleOption     ,
                 charOption       ,
                 stringOption     ,
                 fvectorOption    ,
                 dvectorOption    ,
                 uivectorOption   ,
                 rangeOption      ,
                 outputfileOption ,
                 oldOption        };

/// Pair to identify types of data members
typedef std::pair<void*, OptionType> OptionPair;

/// Map to match a user option to a data member
typedef std::map<std::string, OptionPair> OptionMap;

/// Map to match a program flag to its setting
typedef std::map<std::string, bool> FlagMap;

/// Class to read program options and flags
class OptionReader{
public:
  OptionReader();
  virtual ~OptionReader();
  ///toggle verbosity. true complains to cerr when there are problems
  void setVerbose(bool);
  ///read options, autodetecting file or command args, and set values
  virtual bool ReadUserOptions( int argc, char** argv, const char* fileargIndicator = 0);
  ///read command-line options
  bool ReadCommandLineArgs( int argc, char** argv, const char* fileargIndicator= 0);
  ///read options from file
  bool ReadArgsFromFile(const char* filename);
  ///check if all required options have been specified
  bool CheckRequiredOptions();
  ///set option variables to option values. returns false if any unrecognized options.
  virtual bool SetOptions();
  ///print all available options to stream
  void PrintAllOptions(std::ostream& os)const;

  ///add a long option
  void addOption(const string&, OptionType, void*, bool required = false);


    //-------------------------------------------------------------------------
    // D. Favro: adding a few type-safe versions here -- above non-type-safe
    //	version should not be public.  Actually, the whole thing needs to be
    //	redesigned from the ground up, but for the moment as a hacked-up
    //	workaround, these will help.

    /// Add a long option of type unsigned (type-safe):
    void addOption( const string & name, unsigned & value, unsigned defaultValue = 0 );

    void addOption( const string & name, long & value, long defaultValue = 0 ); ///< Add a long option of type long (type-safe)
    void addOption( const string & name, long & value, bool required ); 	///< Add a long option of type long (type-safe)
    void addOption( const string & name, int & value, int defaultValue = 0 );	///< Add a long option of type int (type-safe)
    void addOption( const string & name, int & value, bool required );		///< Add a long option of type int (type-safe)
    void addOption( const string & name, bool & value, bool defaultValue = false ); ///< Add a long option of type bool (type-safe)

    /// Add an option of type vector (type-safe).
    void addOption( const string & name, std::vector<float> & value, bool required = false )
	{ addOption( name, fvectorOption, &value, required ); }

    /// Add an option of type cvector (type-safe).
    void addOption( const string & name, genepi::cvector<float> & value, bool required = false )
	{ addOption( name, value.getVector_unsafe(), required ); }

    // ...and so on...
    //-------------------------------------------------------------------------


  ///add a short option
  void addOption(char, OptionType, void*, bool required = false);
  ///add a dual(long and short) option
  void addOption(char, const string&, OptionType, void*, bool required = false);
  ///add a flag with short name
  void addFlag(char);
  ///add a flag with long name
  void addFlag(const string&);
  ///add a flag with both long and short names
  void addFlag(char, const string&);

  ///set a use options manually (good for defaults)
  void setUserOption(const string&, const string&);
  ///print user options to file
  virtual void PrintUserOptions(const char* filename);
  //print user options to a stream
  virtual void PrintUserOptions(std::ostream& os);
  ///print user-specified flags to a stream
  void PrintUserFlags(std::ostream& os);
  ///set separator for options file. eg '=' for a = b
  void setFileSeparators(const string& sep);
  ///free some memory by clearing maps
  void clear();
  ///return specified option value
  const string& getUserOption(const string&);
  const string& getUserOption(char);
  ///check if a flag has been specified
  bool getFlag(const string&);
  bool getFlag(char);
  ///check if any options have been specified
  bool hasOptions()const;

protected:
  /// map to hold user options
  std::map<std::string, std::string> useroptions;
  OptionMap ProgOptions;

  /// Setting of each program flag
  FlagMap Flags;
  std::vector<std::string>    userflags;
  std::map<char, std::string> Short2LongMap;
  std::vector<std::string>    RequiredOptions;

private:
  string separators;
  bool hasoptions;
  bool Verbose;

  ///returns the long version of a short option/flag name. If the option/flag doesn't exist, returns a string version of short name.
  string short2Long(char shortname);
  ///set a flag. returns false if any unrecognized flags.
  bool SetFlag(const string&);
  ///assigns a options value to a variable
  int assign(OptionPair& opt, const string& value);
  ///reduces an option to standard form (lower case, no non-graph chars)
  void ParseOptionName(string& name);
  void ReportBadUserOption(std::ostream& os, string& line, unsigned linenum, const char* filename)const;

  // UNIMPLEMENTED: to avoid use
  OptionReader(const OptionReader&);
  OptionReader& operator=(const OptionReader&);
};


/** @} */

END_BCLIB_NAMESPACE


#endif // ! __bclib_OptionReader_h
