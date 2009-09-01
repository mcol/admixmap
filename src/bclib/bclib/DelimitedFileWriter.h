// *-*-C++-*-*
/** 
 *   \file DelimitedFileWriter.h
 *   Class to write delimited files, an alternative to std::ofstream to enable polymorphic file writers
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef DELIMITEDFILEWRITER_H
#define DELIMITEDFILEWRITER_H

#include "bclib/bclib.h"
#include "bclib/Delimitedostream.h"

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */



///extension of Delimited ostream to write delimited files
class DelimitedFileWriter : public Delimitedostream{
public:
  /**
     Default constructor.
     Sets delimiter character, which defaults to a tab.
  */
  DelimitedFileWriter(char sep='\t'):Delimitedostream(&buf), buf(sep){};

  /**
     Constructor with filename.
     Opens output file and sets delimiter.
  */
  DelimitedFileWriter(const char* filename, char sep='\t' ): Delimitedostream(&buf), buf(filename, sep){}

  ///destructor
  virtual ~DelimitedFileWriter(){ }

  /**
     Open output file.
  */
  virtual void open(const char* filename){
    buf.open(filename);
  }

  /**
     Returns true if file is open.
  */
  bool is_open()const{
    return (buf.is_open());
  }

  /**
     Implements newline manipulator by calling the other version of newline on buf.
  */
  virtual void putNewLine(){
    newline(buf);
  }

  /**
     close output file.
  */
  virtual void close(){
    buf.close();
  }

  ///toggle delimiting
  void delimit(bool b){
    buf.delimit(b);
  }

protected:
  DelimitedFileBuf buf;

private:
  ///private to prevent undesired copying
  DelimitedFileWriter(const DelimitedFileWriter& );
  void operator=(const DelimitedFileWriter& );

};


/** @} */

END_BCLIB_NAMESPACE

#endif
