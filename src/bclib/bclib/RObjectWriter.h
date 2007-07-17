// *-*-C++-*-*
/** 
 *   \file RObjectWriter.h
 *   This class encapsulates the details of writing an R object
 *   log-concave distributions. 
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef ROBJECTWRITER_H
#define ROBJECTWRITER_H

#include "bclib/FileWriter.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

///struct used to write a comment in an R object
struct Rcomment{
  Rcomment(const char* c): str(c){};
  std::string str;
};

///class to write an R object
class RObjectWriter : public FileWriter{

public:

  ///no-arg c'tor
  RObjectWriter();
  ///c'tor with filename, equivalent to default c'tor + open
  RObjectWriter(const char*);
  ///c'tor with filename, equivalent to default c'tor + open
  RObjectWriter(const std::string&);
  ///open file and start writing R object
  void open(const char*);
  ///open file and start writing R object
  void open(const std::string&);
  ///finish writing R object by writing dimensions and dimnames
  ///supply a vector of dimensions and a vector of labels for at least the first dimension
  void close(const std::vector<int>& dims, const std::vector< std::vector<std::string> >& dimnames);
  //write a comment, prefixed by a comment character

  /**
     manipulator to start new line.
     differs from endl in that the object does not write the newline until there is more output and it remembers if it needs a comma first.
  */
  friend void newline(RObjectWriter& R );
  ///add a comment
  void comment(const char* );
  ///add a comment
  void comment(const std::string&);
  ///write to file
  RObjectWriter& operator<<(const int);
  RObjectWriter& operator<<(const unsigned);
  RObjectWriter& operator<<(const long);
  RObjectWriter& operator<<(const float);
  RObjectWriter& operator<<(const double);
  RObjectWriter& operator<<(const char);
  RObjectWriter& operator<<(const bool);
  RObjectWriter& operator<<(const char*);
  RObjectWriter& operator<<(const std::string);

private:
  bool needComma;///<used to determine if a comma is needed before a new element
  std::string comment_cache;

  ///write dimensions of R object
  void WriteDimensions(const std::vector<int>&, const std::vector< std::vector<std::string> >&);
  template <class T>
  RObjectWriter& write(const T);

};

//declarations of friend functions
void newline(RObjectWriter& R);

//stream insertion for newline
RObjectWriter&  operator<<(RObjectWriter& R, void (*man)(RObjectWriter&));
//stream insertion for comments
RObjectWriter&  operator<<(RObjectWriter& R, const Rcomment&);

#endif
