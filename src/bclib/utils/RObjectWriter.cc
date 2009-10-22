/** 
 *   \file RObjectWriter.cc
 *   This class encapsulates the details of writing an R object
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "bclib/RObjectWriter.h"
#include <iomanip>

BEGIN_BCLIB_NAMESPACE

void RObjectWriter::WriteFirstLine(){
  //start writing R object
  buf.delimit(false);
  *this << "structure(.Data=c(\n";
  buf.delimit(true);
}

void RObjectWriter::open(const char* filename){
  DelimitedFileWriter::open(filename);
  WriteFirstLine();
}

void RObjectWriter::close(const std::vector<int>& dim, const std::vector<std::vector<std::string> >& DimNames){
  WriteDimensions(dim, DimNames);
  buf.close();
}

//private function
void RObjectWriter::close(){
  buf.close();
}

void RObjectWriter::WriteDimensions(const std::vector<int>& dim, const std::vector<std::vector<std::string> >& DimNames){
    using std::vector;
    using std::string;
    //sanity checks
    if(dim.size() < DimNames.size())
      throw string("Error in RobjectWriter::WriteDimensions : dimensions do not match dimnames");
    
    //turn off delimiting
    buf.delimit( false );
    
    //   //write comment for last line if there is one
    //   if(comment_cache.size()){
    //     file << COMMENT_CHAR << comment_cache;
    //     comment_cache.clear();
    //   }
    
    
    //finish defining data 
  *this << "),\n"
    
    //define dimensions
	<< ".Dim = c(";
  
  for(unsigned int i = 0; i < dim.size(); i++){
    *this << dim[i];
    if(i != dim.size() - 1){
      *this << ",";
    }
  }
  *this << "),\n"// << endl
    
  //define dimnames
	<< ".Dimnames=list(";
  
  unsigned d = 0;
  for(; d < DimNames.size(); ++d){
    if(d > 0)//seperate dimnames with comma
      *this << ", ";    
    
    if(DimNames[d].size()==0){
      *this << "character(0)";
      continue;
    }
    if(dim[d] != (int)DimNames[d].size())
      throw string("Error in RobjectWriter::WriteDimensions : dimnames do not match dimension");

    //write dimnames for labeled dimensions
    *this << "c(";
    for(unsigned int i = 0; i < DimNames[d].size(); i++){
      *this << "\"" << DimNames[d][i] << "\"";
      if(i != DimNames[d].size() - 1){
	*this << ",";
      }
    }
    *this << ")";
  }

  //write 'character(0)' for remaining dimensions
  for( ; d < dim.size(); ++d){
    if(d > 0)//seperate dimnames with comma
      *this << ", ";    

    *this << "character(0)";
  }
  //...and finally
  *this << "))\n";
}

END_BCLIB_NAMESPACE
