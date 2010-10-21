/*
 *   AncestryAssocTest.cc
 *   Class to implement score test for association of trait with ancestry
 *   Copyright (c) 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

//=============================================================================
/// \file AncestryAssocTest.cc
/// Implementation of the AncestryAssocTest class.
//=============================================================================

#include "AncestryAssocTest.h"
#include "Genome.h"
#include "bclib/LogWriter.h"
#include "bclib/DelimitedFileWriter.h"
#include <vector>
#include <string>

AncestryAssocTest::AncestryAssocTest(){
  useprevb = true;
  firstpoplabel = 0;
}

void AncestryAssocTest::Initialise(const char* filename, const int NumStrata, const int NumLoci){
  if(NumStrata == 2){
    firstpoplabel = 1;//skip first pop label if 2 pops
  }
  CopyNumberAssocTest::Initialise(filename, NumStrata, NumLoci);
}

void AncestryAssocTest::Output(const Vector_s& PopLabels, const Genome& Loci){

  for(unsigned int j = 0; j < L; j++ ){
    const string locuslabel = Loci(j)->GetLabel(0);
    for( unsigned k = 0; k < NumOutputStrata; k++ ){//end at 1 for 2pops
      std::string label = locuslabel;
      label += "\",\"" + PopLabels[k+firstpoplabel];//label is output in quotes
      OutputCopyNumberAssocTest(j, k, R, label, false);
    }
  }
  ++numPrintedIterations;
}

void AncestryAssocTest::WriteFinalTable(const char* filename, const Vector_s& PopLabels, 
					const Genome& Loci, bclib::LogWriter& Log){
  bclib::DelimitedFileWriter finaltable(filename);
  Log << bclib::Quiet << "Tests for locus linkage written to " << filename << "\n";
  finaltable <<"Locus" << "Population" << "Score" << "CompleteInfo" << "ObservedInfo" 
	     << "PercentInfo" << "Missing1" << "Missing2" << "StdNormal" << "PValue"
	     << bclib::newline;

  for(unsigned int j = 0; j < L; j++ ){
    const string locuslabel = Loci(j)->GetLabel(0);
    for( unsigned k = 0; k < NumOutputStrata; k++ ){//end at 1 for 2pops
      std::string label = locuslabel;
      label += "\"\t\"" + PopLabels[k+firstpoplabel];//label is output in quotes
      OutputCopyNumberAssocTest(j, k, finaltable, label, true);
    }
  }
  finaltable.close();
}

AncestryAssocTest::~AncestryAssocTest(){
  if(test){
    std::vector<std::vector<std::string> > labels(1);
    labels[0].push_back("Locus");
    labels[0].push_back("Population");  
    labels[0].push_back("minusLog10PValue");

    std::vector<int> dimensions(3,0);
    dimensions[0] = labels[0].size();
    dimensions[1] = L * NumOutputStrata;
    dimensions[2] = (int)(numPrintedIterations);
    
    R.close(dimensions, labels);
  }
}
