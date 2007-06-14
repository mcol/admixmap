// *-*-C++-*-*
/**
   \file HapMapData.h
   functions to write HAPMIXMAP formatted locusfiles and genotypesfiles, based on HapMap data.
   This file is part of FPHD

   This program is free software distributed WITHOUT ANY WARRANTY. 
   You can redistribute it and/or modify it under the terms of the GNU General Public License, 
   version 2 or later, as published by the Free Software Foundation. 
   See the file COPYING for details.

   Copyright (c) David O'Donnell 2007
*/
#ifndef LOCUSFILEWRITER_H
#define LOCUSFILEWRITER_H

#include <string>

class HapMapLegend;

/**
   Write locus file(s).
   \param Legend HapMapLegend object, holding locus info
   \param fileprefix prefix of outputfile
   \param first index of first locus to use
   \param last index of last locus to use
   \beVerbose indicates whether to write information to screen
*/
void WriteLocusFile(HapMapLegend& Legend, const char* fileprefix, 
		    unsigned first, unsigned last, bool beVerbose);

/**
   Write genotypes file(s).
   \param Legend HapMapLegend object, holding locus info
   \param fileprefix prefix of outputfile
   \param dataprefix prefix of HapMap data files (*_phased.txt, *_sample.txt)
   \param first index of first locus to use
   \param last index of last locus to use
   \beVerbose indicates whether to write information to screen
*/
void WriteGenotypesFile(HapMapLegend& Legend, const char* fileprefix, const std::string& dataprefix,  
			unsigned first, unsigned last, bool beVerbose);

///remove monomorphic loci from HapMap data. Call this before reading legend.
void RemoveMonomorphicLoci( const std::string& prefix, bool beVerbose, bool backup);
#endif
