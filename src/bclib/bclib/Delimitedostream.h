// *-*-C++-*-*
/** 
 *   \file Delimitedostream.h
 *   Extension of std::ostream to write delimited streams
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef DELIMITEDSTREAM_H
#define DELIMITEDSTREAM_H

#include "bclib/bclib.h"
#include "bclib/DelimitedBuf.h"

BEGIN_BCLIB_NAMESPACE

///An extension of std::ostream to write delimited output
class Delimitedostream : public std::ostream{
public:
  /**
     Default constructor
  */
  Delimitedostream(){};
  /**
     Constructor with streambuf pointer as argument.
     Otherwise, the std::ostream default is used.
  */
  Delimitedostream(std::streambuf* buf):std::ostream(buf){};

  ///destructor
  virtual ~Delimitedostream(){};

  /**
     Template stream insertion operator.
     Returns reference to Delimitedostream instead of std::ostream
     so the newline manipulator will work.
  */
  template<class T>
  Delimitedostream& operator<<(T t){
    *((std::ostream*)this) << t;
    return *this; 
  }

  /**
     function called by newline.
     Just adds a newline character by default.
     Derived classes should override this with appropraite behaviour.
  */
  virtual void putNewLine(){
    *((std::ostream*)this) << '\n';
  }

  ///set decimal precision to p places
  void setDecimalPrecision(unsigned p);

};

///class to write delimited output to screen
class Delimitedstdout : public Delimitedostream{
public:
  /**
     Constructor.
     Sets delimiter character, which default to a tab.
  */
  Delimitedstdout(char _sep='\t'):Delimitedostream(&buf), buf(_sep){ };

  ///destructor
  virtual ~Delimitedstdout(){ }

private:
  DelimitedStreamBuf buf;
};

/**
   Manipulator to tell Delimitedostream objects to start a new line on next output.
   Cannot just pass '\n' or std::endl as a delimiter may also be pending.
*/
void newline(Delimitedostream& FW);

///stream insertion for newline
Delimitedostream& operator<<(Delimitedostream& FW, void (*manip)(Delimitedostream& ));

///stream insertion for std::endl, std::ends, std::flush
Delimitedostream& operator<<(Delimitedostream& FW, std::ostream& (*manip)(std::ostream& ));

END_BCLIB_NAMESPACE
#endif
