// *-*-C++-*-*
/**
   \file InitialValues.h
   functions to split HAPMIXMAP initial value files
   This file is part of FPHD

   This program is free software distributed WITHOUT ANY WARRANTY. 
   You can redistribute it and/or modify it under the terms of the GNU General Public License, 
   version 2 or later, as published by the Free Software Foundation. 
   See the file COPYING for details.

   Copyright (c) David O'Donnell 2007
*/

#ifndef INITIALVALUES_H
#define INITIALVALUES_H

class HapMapLegend;

void SplitInitialArrivalRateFile(const char* filename, HapMapLegend& Legend);
void SplitInitialMixturePropsFile(const char* filename, HapMapLegend& Legend);
void SplitInitialAlleleFreqFile(const char* filename, HapMapLegend& Legend);
void SplitInitialFreqPriorFile(const char* filename, HapMapLegend& Legend);

#endif
