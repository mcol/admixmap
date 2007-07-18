/** 
 *   \file TableWriter.cc
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
#include "bclib/TableWriter.h"

BEGIN_BCLIB_NAMESPACE

//manipulators
void newline(TableWriter& TW){
  TW.needNewLine = true;
}

TableWriter& operator<<(TableWriter& TW, void (*manip)(TableWriter& )){
  manip(TW);
  return TW;
}

TableWriter::TableWriter(char c){
  sep = c;
  needNewLine = false;
  needSep = false;
}

TableWriter::~TableWriter(){

}

void TableWriter::close(){
  file << std::endl;
  file.close();
}

TableWriter::TableWriter(const char* filename, char c){
  this->open(filename);
  sep = c;
  needNewLine = false;
  needSep = false;
}

template<class T>
TableWriter& TableWriter::write(const T t){
  if(needSep) file << sep;
  if(needNewLine){
    file << std::endl;
    needNewLine = false;
  }
  file << t;
  needSep = true;
  return *this;
}

template TableWriter& TableWriter::write(const unsigned);
template TableWriter& TableWriter::write(const int);
template TableWriter& TableWriter::write(const long);
template TableWriter& TableWriter::write(const float);
template TableWriter& TableWriter::write(const double);
template TableWriter& TableWriter::write(const std::string);
template TableWriter& TableWriter::write(const char*);
template TableWriter& TableWriter::write(const bool);
template TableWriter& TableWriter::write(const char);


TableWriter& TableWriter::operator<<(const int t){
  return write(t);
}
TableWriter& TableWriter::operator<<(const unsigned t){
  return write(t);
}
TableWriter& TableWriter::operator<<(const long t){
  return write(t);
}
TableWriter& TableWriter::operator<<(const float t){
  return write(t);
}
TableWriter& TableWriter::operator<<(const double t){
  return write(t);
}
TableWriter& TableWriter::operator<<(const char t){
  return write(t);
}
TableWriter& TableWriter::operator<<(const bool t){
  return write(t);
}
TableWriter& TableWriter::operator<<(const char* t){
  return write(t);
}
TableWriter& TableWriter::operator<<(const std::string t){
  return write(t);
}

END_BCLIB_NAMESPACE
