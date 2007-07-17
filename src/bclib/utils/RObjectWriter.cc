/** 
 *   \file RObjectWriter.cc
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
#include "bclib/RObjectWriter.h"
#include <sstream>
#include <numeric>
#include <iomanip>
#include <cmath>

using namespace::std;
#define COMMENT_CHAR " # "

//manipulators
void newline(RObjectWriter& R){
  R.needNewLine = true;
}

//class members

RObjectWriter::RObjectWriter() {
  needNewLine = false;
  needComma = false;
}

RObjectWriter::RObjectWriter(const char* filename){
  needNewLine = false;
  needComma = false;
  open(filename);
}

RObjectWriter::RObjectWriter(const string& filename){
  needNewLine = false;
  needComma = false;
  open(filename.c_str());
}

void RObjectWriter::open(const char* filename){
  FileWriter::open(filename);

  //start writing R object
  file << "structure(.Data=c(" << endl;
}

void RObjectWriter::open(const string& filename){
  this->open(filename.c_str());
}

void RObjectWriter::close(const vector<int>& dim, const vector<vector<string> >& DimNames){
  WriteDimensions(dim, DimNames);
  file.close();
}

void RObjectWriter::WriteDimensions(const vector<int>& dim, const vector<vector<string> >& DimNames)
{
  //sanity checks
  if(dim.size() < DimNames.size())
    throw string("Error in RobjectWriter::WriteDimensions : dimensions do not match dimnames");

  //write comment for last line if there is one
  if(comment_cache.size()){
    file << COMMENT_CHAR << comment_cache;
    comment_cache.clear();
  }

    //finish defining data 
  file << endl << ")," << endl

  //define dimensions
	<< ".Dim = c(";

  for(unsigned int i = 0; i < dim.size(); i++){
    file << dim[i];
    if(i != dim.size() - 1){
      file << ",";
    }
  }
  file << ")," << endl
  
  //define dimnames
	<< ".Dimnames=list(";
  
  unsigned d = 0;
  for(; d < DimNames.size(); ++d){
    if(d > 0)//seperate dimnames with comma
      file << ", ";    

    if(DimNames[d].size()==0){
      file << "character(0)";
      continue;
    }
    if(dim[d] != (int)DimNames[d].size())
      throw string("Error in RobjectWriter::WriteDimensions : dimnames do not match dimension");

    //write dimnames for labeled dimensions
    file << "c(";
    for(unsigned int i = 0; i < DimNames[d].size(); i++){
      file << "\"" << DimNames[d][i] << "\"";
      if(i != DimNames[d].size() - 1){
	file << ",";
      }
    }
    file << ")";
  }

  //write 'character(0)' for remaining dimensions
  for( ; d < dim.size(); ++d){
    if(d > 0)//seperate dimnames with comma
      file << ", ";    

    file << "character(0)";
  }
  //...and finally
  file << "))" << endl;

}

void RObjectWriter::comment(const char* c){
  if(needNewLine){
    if(needComma)
      file << ',';
    file << COMMENT_CHAR << comment_cache << endl << COMMENT_CHAR << c;
    comment_cache.clear();
    needNewLine = false;
  }else
    comment_cache.append(c);
}
void RObjectWriter::comment(const std::string& s){
  comment(s.c_str()); 
}

template <class T>
RObjectWriter& RObjectWriter::write(const T t){
  if(needComma) file << ",";
  if(needNewLine){
    if(comment_cache.size()){
      file << COMMENT_CHAR << comment_cache;
      comment_cache.clear();
    }
    file << std::endl;
    needNewLine = false;
  }
  file << t;
  needComma = true;
  return *this;
}

template RObjectWriter& RObjectWriter::write(const unsigned);
template RObjectWriter& RObjectWriter::write(const int);
template RObjectWriter& RObjectWriter::write(const long);
template RObjectWriter& RObjectWriter::write(const float);
template RObjectWriter& RObjectWriter::write(const double);
template RObjectWriter& RObjectWriter::write(const std::string);
template RObjectWriter& RObjectWriter::write(const char*);
template RObjectWriter& RObjectWriter::write(const bool);
template RObjectWriter& RObjectWriter::write(const char);

RObjectWriter& RObjectWriter::operator<<(const int t){
  return write(t);
}
RObjectWriter& RObjectWriter::operator<<(const unsigned t){
  return write(t);
}
RObjectWriter& RObjectWriter::operator<<(const long t){
  return write(t);
}
RObjectWriter& RObjectWriter::operator<<(const float t){
  return write(t);
}
RObjectWriter& RObjectWriter::operator<<(const double t){
  return write(t);
}
RObjectWriter& RObjectWriter::operator<<(const char t){
  return write(t);
}
RObjectWriter& RObjectWriter::operator<<(const bool t){
  return write(t);
}
RObjectWriter& RObjectWriter::operator<<(const char* t){
  return write(t);
}
RObjectWriter& RObjectWriter::operator<<(const std::string t){
  return write(t);
}

RObjectWriter& operator<<(RObjectWriter& R, void (*manip)(RObjectWriter& R)){
  manip(R);
  return R;
}
RObjectWriter& operator<<(RObjectWriter& R, const Rcomment& c){
  R.comment(c.str.c_str());
  return R;
}
