/** 
 *   Options.cc 
 *   This file is part of bcppcl
 *   Class to read program options and flags. 
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "bclib/OptionReader.h"
#include "bclib/StringConvertor.h"
#include "bclib/StringSplitter.h"
#include <string.h>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>

using namespace std;

BEGIN_BCLIB_NAMESPACE

OptionReader::OptionReader(){
  separators = "=";
  Verbose = true;
  hasoptions = false;
}

bool OptionReader::ReadUserOptions(int argc,  char** argv, const char* fileargIndicator){
  bool allClear = true;
  if(argc == 1)//no args
    allClear = true;
  //one arg, unless with a - prefix
  else if(argc == 2 && strncmp(argv[1], "-", 1)){
    allClear = ReadArgsFromFile(argv[1]);
  }
  else if(argc > 1){
    allClear = ReadCommandLineArgs(argc, argv, fileargIndicator);
  }
  if(!SetOptions())
    allClear = false;
  return allClear;
}

void OptionReader::setVerbose(bool v){
  Verbose = v;
}

//helper function to convert a char array to a string
string charArray2String(const char* a){
  string s(a);
  return s;
}
//helper function to create a sring with a single char
string char2String(char c){
  stringstream ss;
  ss << c;
  return ss.str();
}



//-----------------------------------------------------------------------------
// D. Favro: slight improvement.  Needs a lot more work.
//-----------------------------------------------------------------------------

void OptionReader::addOption( const string & optName, long & value, long defaultValue )
    {
    value=defaultValue;
    addOption( optName, longOption, &value, false );
    }

void OptionReader::addOption( const string & optName, long & value, bool required )
    {
    addOption( optName, longOption, &value, required );
    }

void OptionReader::addOption( const string & optName, int & value, int defaultValue )
    {
    value=defaultValue;
    addOption( optName, intOption, &value, false );
    }

void OptionReader::addOption( const string & optName, int & value, bool required )
    {
    addOption( optName, intOption, &value, required );
    }

void OptionReader::addOption( const string & optName, bool & value, bool defaultValue )
    {
    value=defaultValue;
    addOption( optName, boolOption, &value, false );
    }

//-----------------------------------------------------------------------------



void OptionReader::addOption(const string& name, OptionType otype, void* address, bool required){
  ProgOptions[name] = OptionPair(address, otype);
  if(required)RequiredOptions.push_back(name);
}
void OptionReader::addOption(const char* name, OptionType otype, void* address, bool required){
  addOption(charArray2String(name), otype, address, required);
}
void OptionReader::addOption(char shortname, const string& longname, OptionType otype, void* address, bool required){
  ProgOptions[longname] = OptionPair(address, otype);
  Short2LongMap[shortname] = longname;
  if(required)RequiredOptions.push_back(longname);
}
void OptionReader::addOption(char shortname, const char* longname, OptionType otype, void* address, bool required){
  addOption(shortname, charArray2String(longname), otype, address, required);
}
void OptionReader::addOption(char shortname, OptionType otype, void* address, bool required){
  stringstream opt;
  opt << shortname;
  const string longname = opt.str();
  ProgOptions[longname] = OptionPair(address, otype);
  Short2LongMap[shortname] = longname;
  if(required)RequiredOptions.push_back(longname);
}
void OptionReader::addFlag(char shortname){
  stringstream opt;
  opt << shortname;
  Flags[opt.str()] = false;
  Short2LongMap[shortname] = opt.str();
}
void OptionReader::addFlag(const string& longname){
  Flags[longname] = false;
}
void OptionReader::addFlag(const char* longname){
  addFlag(charArray2String(longname));
}
void OptionReader::addFlag(char shortname, const string& longname){
  Flags[longname] = false;
  Short2LongMap[shortname] = longname;
}
void OptionReader::addFlag(char shortname, const char* longname){
  addFlag(shortname, charArray2String(longname));
}

void OptionReader::setFileSeparators(const string& sep){
  separators = sep;
}

OptionReader::~OptionReader(){
}

void OptionReader::clear(){
  useroptions.clear();
  ProgOptions.clear();
  Short2LongMap.clear();
  Flags.clear();
  userflags.clear();
}

///convert option name to lower case and remove underscores
void OptionReader::ParseOptionName(string& name){
  StringConvertor::toLower(name);
  StringConvertor::RemoveCharFromString(name, '_');
}

void OptionReader::ReportBadUserOption(ostream& os, string& line, unsigned linenum, const char* filename)const{
  os << "** Invalid option: \"" << line << "\" on line " << linenum << " of " << filename << endl;
}


/// Strip leading and trailing whitespace from a string (static helper
/// function).  Returns true if entire string was whitespace.
static bool stripWS( string & str )
    {
    static const char wsChars [] = " \t\n\r";
    const string::size_type firstCh = str.find_first_not_of( wsChars, 0, sizeof(wsChars)-1 );
    if ( firstCh == string::npos )
	return true; // **** RETURN HERE ****
    const string::size_type lastCh = str.find_last_not_of( wsChars );
    if ( (lastCh == string::npos) || (lastCh < firstCh) )
	throw string( "insanity here" );
    str = str.substr( firstCh, lastCh + 1 - firstCh );
    return false;
    }



static istream& mgetline( istream& in, string & str, char delim )
    {
    char ch;
    //str.erase();
    while ( in.get( ch ) )
	{
	if ( ch == delim )
	    return in;
	else
	    str += ch;
	}
    return in;
    }


///read user options from file
bool OptionReader::ReadArgsFromFile(const char* filename){
  if ( (filename == 0) || (strlen(filename) ==0) )
    return false;

  ifstream fin(filename);

  if (!fin.is_open()) {
    if(Verbose){
      cerr << "Cannot open file \"" << filename << "\"." << endl;
      //exit(1);
    }
    return false;
  }

  std::string str;
  bool badoptions = false;
  unsigned linenum = 0;
  //read in line from file
  while (mgetline(fin,str,'\n')){// ## apparent memory leak 

    ++linenum;

    //ignore #comments
    // ****WARNING**** THIS WILL BREAK IF A '#' APPEARS WITHIN A QUOTED STRING!
    string::size_type comChPos = str.find_first_of("#");
    if( comChPos != string::npos )
      str.erase( comChPos );

    //skip blank lines
    if(str.find_first_not_of(" \t\n\r") < str.length() )
      {   
	// 	 //check there is an '=' on line
	// 	 if(str.find("=") == string::npos || ){
	// 	   ReportBadUserOption(cerr, str, linenum, filename);
	// 	   badoptions = true;
	// 	   continue;
	// 	 }
	//split on '='
	vector<string> tokens;
	StringSplitter::Tokenize(str, tokens, separators);
	//check we have 2 tokens, option name and value
	if(tokens.size() != 2){
	  if(Verbose)
	    ReportBadUserOption(cerr, str, linenum, filename);
	  badoptions = true;
	  continue;
	}

	// Check that tokens have graphical characters (not just whitespace),
	// and strip leading and trailing whitespace:
	if ( stripWS(tokens[0]) || stripWS(tokens[1]) ) {
	  if ( Verbose )
	    ReportBadUserOption(cerr, str, linenum, filename);
	  badoptions = true;
	  continue;
	}
	hasoptions = true;//we have at least one option

	//convert option name to required format, lowercase and no underscores
	ParseOptionName(tokens[0]);
	
	//now all whitespace has been removed, add option to map of user options
	useroptions[tokens[0]] = tokens[1];
      }
    //clear str ready for next time
    str.clear();
  }
  //close file
  fin.close();
  
  return !badoptions;
}

bool OptionReader::ReadCommandLineArgs(const int argc, char** argv, const char* fileargIndicator){
  const string delims = "-" + separators;
  bool badOptions = false;
  for(int i = 1; i < argc; ++i){
    
    if(fileargIndicator){
      //check for options file flag and read options from file
      if(!strncmp(argv[i], fileargIndicator, 2)){
	const char* optionsFilename= 0;
	if(strlen(argv[i])==strlen(fileargIndicator)){//space between indicator and the filename
	  ++i;//move to the next arg and interpret as a filename
	if(argv[i][0] == '-'){//oops, this is actually a new arg
	  cerr << "ERROR: " << fileargIndicator << " requires an options filename as argument\n";
	  exit(1);
	}
	optionsFilename = argv[i];
	}
	else{//no space
	  optionsFilename = argv[i]+2;
	}
	//cout << "Reading options from " << optionsFilename << endl;
	if(!ReadArgsFromFile(optionsFilename))
	  badOptions = true;
	continue;
      }
    }

    vector<string> args;
    string argstring = argv[i];
    argstring.erase(0, argstring.find_first_not_of("-"));//remove leading '-'s
    if(argstring.size()==0)
      continue;

    StringSplitter::Tokenize(argstring, args, separators);
    if(args.size()==1){//parse as flag
      hasoptions = true;//we have at least one option
      if(!SetFlag(args[0]))
	badOptions = true;
    }
    else if(args.size()==2){//parse as option
      hasoptions = true;//we have at least one option
      if(args[0].size()==1)//short name
	//store long version
	useroptions[short2Long(args[0][0])] = args[1];
      else//long name
	useroptions[args[0]] = args[1];
    }
    else{
      if(Verbose)
	cerr << "ERROR: mismatched arguments: " << argv[i] << endl;
      //exit(1);
      badOptions = true;
    }
  }
  return !badOptions;
}

string OptionReader::short2Long(char shortname){
  if(Short2LongMap.find(shortname) == Short2LongMap.end())
    //not recognized - turn into a string. Unrecognized options/flags will be deleted
    return char2String(shortname);
  else
    return Short2LongMap[shortname];

}

//return codes for option assignment
#define OPTION_ASSIGNED 0
#define UNRECOGNIZED_OPTION_TYPE 1
#define ERASE_OPTION 2

/// sets the value of a data member corresponding to a user option.
/// converts the value string to an appropriate type, with method indicated by opt.second
/// and assigns converted value to the data member with address opt.first
int OptionReader::assign(OptionPair& opt, const string& value){
  if(opt.second == intOption)
    *((int*)opt.first) = atoi(value.c_str());
  else if(opt.second ==longOption)
    *((long*)opt.first) = atoi(value.c_str());
  else if(opt.second == doubleOption)
    *((double*)opt.first) = atof(value.c_str());
  else if(opt.second ==floatOption)
    *((float*)opt.first) = atof(value.c_str());
  else if(opt.second == boolOption){
    const int b = atoi(value.c_str());
    if ( b == 1 ) *((bool*)opt.first) = true;
    else if(b == 0) *((bool*)opt.first) = false;
    else return ERASE_OPTION;
  }
  else if(opt.second == charOption)
    *((char*)opt.first) = value[0];

  else if(opt.second == stringOption)
    *((string*)opt.first) = value;

  else if(opt.second == dvectorOption){
    StringConvertor::StringToVec(value, *((vector<double>*)opt.first));
  }
  else if(opt.second == fvectorOption){
    StringConvertor::StringToVec(value, *((vector<float>*)opt.first));
  }
  else if(opt.second == uivectorOption){
    StringConvertor::StringToVec(value, *((vector<unsigned>*)opt.first));
  }
  else if(opt.second == rangeOption){//may be range or list of numbers
    vector<unsigned>* range = (vector<unsigned>*)opt.first;
    string::size_type colon = value.find(":" , 0);
    if(colon != string::npos){//read as range of numbers
      for(unsigned i = (unsigned)atoi((value.substr(0, colon)).c_str()); i <= (unsigned)atoi((value.substr(colon+1)).c_str()); ++i)
	range->push_back(i);
    }
    else StringConvertor::StringToVec(value, *range);//read as list
  }
  else if(opt.second == oldOption){//deprecated option - return signal to erase
    return ERASE_OPTION;
  }
  else if(opt.second != nullOption && opt.second != outputfileOption){
    //skipping output file names as they are set later, prefixing resultsdir
    return UNRECOGNIZED_OPTION_TYPE;//unrecognised option
  }
  return OPTION_ASSIGNED;//success
}

bool OptionReader::SetOptions(){
  //parse user options
  bool badOptions = false;
  vector<string> OptionsToErase;
  for(map<string, string>::iterator i = useroptions.begin(); i != useroptions.end(); ++i){
    if(ProgOptions.find(i->first)!=ProgOptions.end()){//if option has been defined
      int status = assign(ProgOptions[i->first], i->second);
      if (status==ERASE_OPTION) OptionsToErase.push_back(i->first);
    }
    else{
      if(Verbose){
	cerr << "Unknown option: " << i->first
	     << " with arg: " << i->second
	     << endl;
      }
      badOptions = true;
    }
  }
  for(vector<string>::const_iterator e = OptionsToErase.begin(); e != OptionsToErase.end(); ++e)
    useroptions.erase(*e);
  return !badOptions;
}

bool OptionReader::SetFlag(const string& f){
  string longname;
  bool badOptions = false;
  if(f.size()==1){//short name
    //store long version
    longname = short2Long(f[0]);
  }
  else//long name
    longname = f;

  if(Flags.find(longname) != Flags.end()){
    Flags[longname] = true;
    userflags.push_back(f);
  }
  else{//unrecognized flag
    if(Verbose)
      cerr << "Unknown flag: " << f << endl;
    badOptions = true;
  }
  return !badOptions;
}

void OptionReader::setUserOption(const string& name, const string& value){
  useroptions[name] = value;
}

bool OptionReader::CheckRequiredOptions(){
  bool allClear = true;
  for(vector<string>::const_iterator i = RequiredOptions.begin(); i!= RequiredOptions.end(); ++i){
    if(useroptions.find(*i) == useroptions.end()){//option not specified
      allClear = false;
      if(Verbose){
	cerr << "Required option '" << *i << "' not specified" << endl;
      }
    }
  }
  return allClear;
}

void OptionReader::PrintUserOptions(const char* filename){
//   string ss;
//   ss = ResultsDir + "/args.txt";
//   ofstream argstream(ss.c_str());

  ofstream argsfile(filename);
  PrintUserOptions(argsfile);
  argsfile.close();
}
///output Options table to stream
void OptionReader::PrintUserOptions(ostream& os){
  for( map<string, string>::iterator p= useroptions.begin(); p!=useroptions.end(); p++) {
    os << p->first << "=" << p->second <<endl;
  }
}
void OptionReader::PrintUserFlags(ostream& os){
  for(vector<string>::const_iterator f = userflags.begin(); f != userflags.end(); ++f){
    os << *f << endl;
  }
}

//print all available options to cout
void OptionReader::PrintAllOptions(ostream& os)const{
  const string typesAsStrings[] = {"null", "bool", "int", "long", "float", "double", "char", "string", "fvector", "dvector", "uivector", "range", "outputfile", "old"};

  os << "This is a list of all valid options:" << endl
       << "----------------------------------------" << endl;
  for( OptionMap::const_iterator p = ProgOptions.begin(); p != ProgOptions.end(); ++p){
    os << p->first << " ( " << typesAsStrings[p->second.second] << " )" << endl;
  }
  for(map<string, bool>::const_iterator f = Flags.begin(); f != Flags.end(); ++f){
    os << f->first << " ( flag ) " << endl;
  }

}

static string NullString = "";

const string& OptionReader::getUserOption(const string& name){
  if(ProgOptions.find(name) != ProgOptions.end() && useroptions.find(name) != useroptions.end())
    return useroptions[name];
  else//not a valid option or not specified
    return NullString;
}
const string& OptionReader::getUserOption(char name){
  if(Short2LongMap.find(name) != Short2LongMap.end())
    return getUserOption(Short2LongMap[name]);
  else//not a valid option 
    return NullString;
}
bool OptionReader::getFlag(const string& name){
  if(Flags.find(name) != Flags.end())
    return Flags[name];
  else
    return false;
}
bool OptionReader::getFlag(char name){
  if(Short2LongMap.find(name) != Short2LongMap.end())
    return getFlag(Short2LongMap[name]);
  else
    return false;
}
bool OptionReader::hasOptions()const{
  return hasoptions;
}

END_BCLIB_NAMESPACE
