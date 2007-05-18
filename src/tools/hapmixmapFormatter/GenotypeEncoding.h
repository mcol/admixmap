// *-*-C++-*-*
/**
   GenotypeEncoding.h
   Prototypes of various HAPMIXMAP data formatting functions
   This file is part of FPHD

   This program is free software distributed WITHOUT ANY WARRANTY. 
   You can redistribute it and/or modify it under the terms of the GNU General Public License, 
   version 2 or later, as published by the Free Software Foundation. 
   See the file COPYING for details.

   Copyright (c) 2007 David O'Donnell

*/
#ifndef FORMATTING_H
#define FORMATTING_H
#include <string>
#include "HapMapLegend.h"

///determines if a given string is one of a vector of strings
bool isListedSNP(const std::string HapMapSNP, const std::vector<std::string>SNPs);

///returns the complement of a given base (A/T, C/G)
char getComplement(char a);

///returns the complement of a pair of bases
std::pair<char, char> getComplement(std::pair<char, char> a);

///turns a genotype encoded as a string of 2 bases into a pair of bases
std::pair<char, char> GenotypeString2Pair(const std::string& g);

///Encodes a string containing a pair of bases as a genotype string
std::string getGenotype(const std::string& g, std::pair<char, char> a, char Missing);

///Encodes a pair of bases as a genotype string
std::string getGenotype(const std::pair<char, char>& g, std::pair<char, char> a, char Missing);

//Encodes genotypes given in infilename, sorts and writes to outfilename
unsigned EncodeGenotypes(HapMapLegend& Legend,
                         const char* infilename, const char* outfilename, char missing);

#endif
