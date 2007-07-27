/** 
 *   \file Delimitedostream.cc
 *   An extension std::ostream to write delimited output
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "bclib/Delimitedostream.h"
#include <iomanip>

BEGIN_BCLIB_NAMESPACE

void newline(Delimitedostream& FW){
  FW.putNewLine();
}

Delimitedostream& operator<<(Delimitedostream& FW, void (*manip)(Delimitedostream& )){
  manip(FW);
  return FW;
}

Delimitedostream& operator<<(Delimitedostream& FW, std::ostream& (*manip)(std::ostream& )){
  manip(FW);
  return FW;
}

void Delimitedostream::setDecimalPrecision(unsigned p){
  *this << std::setfill(' ');
  setf(std::ios::fixed); 
  precision(p);
  width(p);
}

END_BCLIB_NAMESPACE
