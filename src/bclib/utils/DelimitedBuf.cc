/** 
 *   \file DelimitedBuf.cc
 *   Extensions of streambuf filebuf class to write delimited output
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "bclib/DelimitedBuf.h"

BEGIN_BCLIB_NAMESPACE


std::streamsize DelimitedStreamBuf::xsputn (const char* s,
					    std::streamsize num) {
  
  std::streamsize ss = std::streambuf::xsputn(s, num);
  putchar(sep);
  
  return ss; 
}


DelimitedFileBuf::~DelimitedFileBuf(){
  if(needNewLine)
    //putchar('\n');
    std::filebuf::overflow('\n');
  close();
}

void DelimitedFileBuf::open(const char* filename){
  std::filebuf::open(filename, std::ios::out);
  if(!is_open()){
    std::string error_string = "ERROR: could not open ";
    error_string.append(filename);
    throw(error_string);
  }
}

std::streamsize DelimitedFileBuf::xsputn (const char* s,
					  std::streamsize num) {
  if(delim){
    if(needSep)std::filebuf::overflow(sep);
    needSep = true;
  }
  
  if(needNewLine){
    std::filebuf::overflow('\n');
    needNewLine = false;
  }
  
  return std::filebuf::xsputn(s, num);
}

END_BCLIB_NAMESPACE
