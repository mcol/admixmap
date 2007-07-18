// *-*-C++-*-*
/** 
 *   \file TableWriter.h
 *   Class to write a table to file
 *   log-concave distributions. 
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef TABLE_WRITER_H
#define TABLE_WRITER_H

#include "FileWriter.h"

BEGIN_BCLIB_NAMESPACE

///class to write tables to file
class TableWriter : public FileWriter{
public:
  TableWriter(char='\t');
  TableWriter(const char*, char = '\t');
  ~TableWriter();
  void close();

  TableWriter& operator<<(const int);
  TableWriter& operator<<(const unsigned);
  TableWriter& operator<<(const long);
  TableWriter& operator<<(const float);
  TableWriter& operator<<(const double);
  TableWriter& operator<<(const char);
  TableWriter& operator<<(const bool);
  TableWriter& operator<<(const char*);
  TableWriter& operator<<(const std::string);

  friend void newline(TableWriter& R );

private:
  char sep;///<separator character
  bool needSep;///<used to determine if a separator is needed before a new element

  template<class T>
  TableWriter& write(const T t);

};
//declarations of friend functions
void newline(TableWriter& TW);

//stream insertion for newline
TableWriter&  operator<<(TableWriter& TW, void (*man)(TableWriter&));

END_BCLIB_NAMESPACE
#endif
