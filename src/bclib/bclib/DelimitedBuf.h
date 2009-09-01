// *-*-C++-*-*
/** 
 *   \file DelimitedBuf.h
 *   Extensions of streambuf filebuf class to write delimited output
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef DELIMITEDBUF_H
#define DELIMITEDBUF_H

#include "bclib/bclib.h"
#include <streambuf>
#include <fstream>
#include <string>

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */



///extension of std::streambuf to write delimited output
class DelimitedStreamBuf : public std::streambuf
{
public:
  /**
     Constructor.
     Sets the delimiter character, which defaults to a tab
  */
  DelimitedStreamBuf(char _sep = '\t'):sep(_sep){}

  /**
     Write a character.
     This function in base class must be overridden.
     Just calls putchar.
  */
  virtual int_type overflow (int_type c) {
    return putchar(c);
  }

  /**
   *  @brief  Multiple character insertion.
   *  @param  s  A buffer area.
   *  @param  num  Maximum number of characters to write.
   *  @return  The number of characters written.
   *  @see std::streambuf::xsputn

   Before writing, adds delimiter and/or new line as necessary.
   */
  virtual
  std::streamsize xsputn (const char* s, std::streamsize num) ;

protected:
  ///delimiter character
  const char sep;

};

///extension of std::filebuf to write delimited files
class DelimitedFileBuf : public std::filebuf
{
public:
  /**
     Default constructor; optionally sets delimiter character.
  */
  DelimitedFileBuf(char _sep = '\t'):sep(_sep), needSep(false), needNewLine(false), delim(true){
  }

  /**
     Open-file constructor.
     Opens output file and optionally sets delimiter character.
  */
  DelimitedFileBuf(const char* filename, char _sep = '\t'): sep(_sep), needSep(false), needNewLine(false), delim(true){
    open(filename);
  }

  ///destructor
  virtual ~DelimitedFileBuf();

  /**
     Opens output file.
  */
  virtual void open(const char* filename);

  /**
     Toggles delimiting on/off
  */
  void delimit(bool b){
    delim = b;
  }

  ///friend function for use as a manipulator like std::endl
  friend void newline(DelimitedFileBuf& R );

protected:
  const char sep;
  bool needSep;
  bool needNewLine;
  bool delim;

  /**
     Write multiple characters.
     This is the function that makes the class work by writing newlines and delimiters where necessary.
  */
  virtual
  std::streamsize xsputn (const char* s,
			  std::streamsize num) ;

};

/**
   Manipulator to tell DelimitedFileBuf objects to start a new line on next output.
   Cannot just pass '\n' or std::endl as a delimiter may also be pending.
*/
inline void newline(DelimitedFileBuf& tb ){
  tb.needNewLine = true;
}


/** @} */

END_BCLIB_NAMESPACE


#endif
