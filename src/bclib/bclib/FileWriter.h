// *-*-C++-*-*
/** 
 *   \file FileWriter.h
 *   This class acts as a wrapper for std::ofstream to enable polymorphic file writers
 *   log-concave distributions. 
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef FILEWRITER_H
#define FILEWRITER_H

#include "bclib/bclib.h"
#include <fstream>
#include <string>

BEGIN_BCLIB_NAMESPACE

///wrapper class for std::ofstream
class FileWriter{
public:
  FileWriter();
  FileWriter(const char*);
  virtual ~FileWriter();
  virtual void open(const char*);
  virtual void open(const std::string&);
  virtual void close();
  bool is_open();

  //set decimal precision
  void setDecimalPrecision(unsigned );

  virtual FileWriter& operator<<(const int);
  virtual FileWriter& operator<<(const unsigned);
  virtual FileWriter& operator<<(const long);
  virtual FileWriter& operator<<(const float);
  virtual FileWriter& operator<<(const double);
  virtual FileWriter& operator<<(const char);
  virtual FileWriter& operator<<(const bool);
  virtual FileWriter& operator<<(const char*);
  virtual FileWriter& operator<<(const std::string);

  friend void newline(FileWriter&  );
 
protected:
  std::ofstream file;
  bool needNewLine;///< used to determine if a new line is needed
private:
  template<class T>
  FileWriter& write(const T t);
};

//declarations of friend functions
void newline(FileWriter& FW);

//stream insertion for newline
FileWriter&  operator<<(FileWriter& FW, void (*man)(FileWriter&));

END_BCLIB_NAMESPACE
#endif
