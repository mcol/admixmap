// *-*-C++-*-*
/** 
 *   AncestryAssocTest.h 
 *   Class to implement score test for association of trait with ancestry
 *   Copyright (c) 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef ANCESTRYASSOCTEST_H
#define ANCESTRYASSOCTEST_H

#include "CopyNumberAssocTest.h"
class Genome;
namespace bclib{
  class LogWriter;
}

class AncestryAssocTest : public CopyNumberAssocTest{

public:
  AncestryAssocTest();
  ~AncestryAssocTest();
  void Initialise(const char* filename, const int NumPopulations, const int NumLoci);
  void Output(const Vector_s& PopLabels, const Genome& Loci);
  void WriteFinalTable(const char* filename, const Vector_s& PopLabels, 
		       const Genome& Loci, bclib::LogWriter& Log);
private:
  unsigned firstpoplabel;

};

#endif
