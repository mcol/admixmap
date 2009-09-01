// *-*-C++-*-*
/** 
 *   \file RObjectWriter.h
 *   This class encapsulates the details of writing an R object
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

#include "bclib/DelimitedFileWriter.h"
#include "DelimitedBuf.h"
#include <vector>
#include <string>

BEGIN_BCLIB_NAMESPACE

/** \addtogroup bclib
 * @{ */



///extension of DelimitedFileWriter to write an array as an R object
class RObjectWriter : public DelimitedFileWriter{
public:
  /**
     Default constructor.
     Does nothing except set comma as delimiter.
  */
  RObjectWriter():DelimitedFileWriter(','){};

  /**
     Constructor with filename.
     Opens output file and sets comma as delimiter.
  */
  RObjectWriter(const char* filename ):DelimitedFileWriter(filename, ',') {
    WriteFirstLine();
  };

  /**
     opens output file.
  */
  void open(const char* filename);

  ///manipulator to call addNewLine function
  friend void newline(RObjectWriter& R );

  /**
     Close oputput file.
     Writes the .dim and .dimnames arguments of the array then closes the file.
     \param dim vector of dimensions of the array
     \param DimNames 'list' (vector of vectors) of dimension names. May be incomplete. If empty, 'character(0)' is written instead.
  */
  void close(const std::vector<int>& dim, const std::vector<std::vector<std::string> >& DimNames);

private:

  /**
     Private close function to prevent closing without writing dimensions.
  */
  void close();

  ///write first line of object definition
  void WriteFirstLine();
  ///write dimensions and dimension names
  void WriteDimensions(const std::vector<int>& dim, const std::vector<std::vector<std::string> >& DimNames);

  //to prevent undesired copying
  RObjectWriter(const RObjectWriter& );
  void operator=(const RObjectWriter& );
};


/** @} */

END_BCLIB_NAMESPACE

#endif
