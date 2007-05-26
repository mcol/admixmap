// *-*-C++-*-*
/** 
 *   HAMIXMAP
 *   HapMixAllelicAssocTest.h 
 *   wrapper for the hapmix model allelic assoc test
 *   Extension of the CopyNumberAssocTest   
 *   
 *   Copyright (c) 2007 David O'Donnell
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#ifndef HAPMIXALLELICASSOCTEST_H
#define HAPMIXALLELICASSOCTEST_H

#include "CopyNumberAssocTest.h"
class InputHapMixData;
class LogWriter;

class HapMixIndividualCollection;
class Regression;

///test using conditional distribution of copies of allele2, rather than sampled values
class HapMixAllelicAssocTest : public CopyNumberAssocTest{
public:
  void Update(const HapMixIndividualCollection* const IC, const Regression* const R, const Genome& Loci);
  void PrintAverageInfo(LogWriter& Log, const InputHapMixData& data, const char* filename);
};

#endif

