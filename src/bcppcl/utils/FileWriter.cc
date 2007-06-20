/** 
 *   \file FileWriter.cc
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
#include "bcppcl/FileWriter.h"
#include <iomanip>

//manipulators
void newline(FileWriter& FW){
 FW.needNewLine = true;
}

FileWriter& operator<<(FileWriter& FW, void (*manip)(FileWriter& )){
  manip(FW);
  return FW;
}

FileWriter::FileWriter(){
  needNewLine = false;
}

FileWriter::~FileWriter(){
  if(file.is_open())
    file.close();
}

FileWriter::FileWriter(const char* filename){
  needNewLine = false;
  this->open(filename);
}

void FileWriter::open(const char* filename){
  file.open(filename);
  if(!file.is_open()){
    std::string error_string = "ERROR: could not open ";
    error_string.append(filename);
    throw(error_string);
  }
}

void FileWriter::open(const std::string& filename){
  open(filename.c_str());
}

void FileWriter::close(){
  file.close();
}

bool FileWriter::is_open()const{
  return file.is_open();
}

void FileWriter::setDecimalPrecision(unsigned p){
  file << std::setfill(' ');
  file.setf(std::ios::fixed); 
  file.precision(p);
  file.width(p);
}

template<class T>
FileWriter& FileWriter::write(const T t){
  if(needNewLine){
    file << std::endl;
    needNewLine = false;
  }
  file << t;
  return *this;
}

template FileWriter& FileWriter::write(const unsigned);
template FileWriter& FileWriter::write(const int);
template FileWriter& FileWriter::write(const long);
template FileWriter& FileWriter::write(const float);
template FileWriter& FileWriter::write(const double);
template FileWriter& FileWriter::write(const std::string);
template FileWriter& FileWriter::write(const char*);
template FileWriter& FileWriter::write(const bool);
template FileWriter& FileWriter::write(const char);

FileWriter& FileWriter::operator<<(const int t){
  return write(t);
}
FileWriter& FileWriter::operator<<(const unsigned t){
  return write(t);
}
FileWriter& FileWriter::operator<<(const long t){
  return write(t);
}
FileWriter& FileWriter::operator<<(const float t){
  return write(t);
}
FileWriter& FileWriter::operator<<(const double t){
  return write(t);
}
FileWriter& FileWriter::operator<<(const char t){
  return write(t);
}
FileWriter& FileWriter::operator<<(const bool t){
  return write(t);
}
FileWriter& FileWriter::operator<<(const char* t){
  return write(t);
}
FileWriter& FileWriter::operator<<(const std::string t){
  return write(t);
}

