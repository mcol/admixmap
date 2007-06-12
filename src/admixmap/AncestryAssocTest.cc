/*
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
#include "AncestryAssocTest.h"
#include <vector>
#include <string>

AncestryAssocTest::AncestryAssocTest(){
  useprevb = true;
}

void AncestryAssocTest::OpenOutputFile(LogWriter &Log, const char* filename){
  OpenFile(Log, &outputfile, filename, "Tests for locus linkage", true);
}

void AncestryAssocTest::ROutput(){
  if(test){
    int KK = K;
    if(KK ==2 )KK = 1;
    
    std::vector<std::string> labels;
    labels.push_back("Locus");
    labels.push_back("Population");  
    labels.push_back("minusLog10PValue");

    std::vector<int> dimensions(3,0);
    dimensions[0] = labels.size();
    dimensions[1] = L * KK;
    dimensions[2] = (int)(numPrintedIterations);
    
    R_output3DarrayDimensions(&outputfile, dimensions, labels);
  }
}
