/** 
 *   HAMIXMAP
 *   HapMixAllelicAssocTest.cc
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
#include "HapMixAllelicAssocTest.h"
#include "HapMixIndividualCollection.h"
#include "InputHapMixData.h"
#include "bcppcl/Regression.h"

HapMixAllelicAssocTest::HapMixAllelicAssocTest(){
  useprevb = false;
  K = 1;
}

//public function
void HapMixAllelicAssocTest::Initialise(const char* filename, const int NumLoci, LogWriter &Log){
  this->Initialise(filename, 1, NumLoci, Log);
}

//private function - calls base function
void HapMixAllelicAssocTest::Initialise(const char* filename, const int , const int NumLoci, LogWriter &Log){
  CopyNumberAssocTest::Initialise(filename, 1, NumLoci, Log);
}

void HapMixAllelicAssocTest::OpenOutputFile(LogWriter &Log, const char* filename){
  OpenFile(Log, &outputfile, filename, "Tests for allelic association", true);
}
void HapMixAllelicAssocTest::Update(const HapMixIndividualCollection* const IC, const Regression* const R, const Genome& Loci){
  this->Reset();
  const double dispersion = R->getDispersion();
  const double* const EY = R->getExpectedOutcome();
  const int NumberOfIndividuals = IC->getNumberOfIndividualsForScoreTests();

  for( int i = 0 ; i < NumberOfIndividuals; i++ ){
    this->UpdateB(R->DerivativeInverseLinkFunction(i), dispersion, 0);
  }
  const int offset = IC->getFirstScoreTestIndividualNumber();
  
  for( int i = 0 ; i < NumberOfIndividuals; i++ ){
    
    const HapMixIndividual* const ind = IC->getHapMixIndividual(i + offset);
    vector<vector<double> > UnorderedProbs(3, vector<double>(1));
    const double YMinusEY = IC->getOutcome(0, i) - EY[i];//individual outcome - its expectation
    double DInvLink = R->DerivativeInverseLinkFunction(i);
    unsigned int numberCompositeLoci = Loci.GetNumberOfCompositeLoci();
    for(unsigned int j = 0; j < numberCompositeLoci; j++ ) {
      // SNPs only
      if (Loci.GetNumberOfStates(j) == 2) {
	UnorderedProbs = ind->getUnorderedProbs(j); 
	CopyNumberAssocTest::Update(j, 0/*no covariates*/, dispersion, YMinusEY, DInvLink, UnorderedProbs);
      }
    }
  }

  this->Accumulate();
}

void HapMixAllelicAssocTest::ROutput(){
  if(test){
    vector<string> labels;
    labels.push_back("Locus");
    labels.push_back("minusLog10PValue");

    vector<int> dimensions(3,0);
    dimensions[0] = labels.size();
    dimensions[1] = L;
    dimensions[2] = (int)(numPrintedIterations);
    
    R_output3DarrayDimensions(&outputfile, dimensions, labels);
  }
}

#define COUNTINFO(S, X) \
  if(S == "NA" || S =="NaN") \
    X +=0.0;\
  else \
    X += atof(S.c_str());

void HapMixAllelicAssocTest::PrintAverageInfo(LogWriter& Log, const InputHapMixData& data, const char* filename){
  //filename is the final table written by this class earlier
  ifstream ScoreTable(filename);
  //TODO?? check this file exists. It should.

  string scrap;
  float sumPI = 0.0, sumMissing1 = 0.0, sumMissing2 = 0.0;
  string sInfo, smissing1, smissing2;
  getline(ScoreTable, scrap);//skip header
  ScoreTable >> scrap >> scrap >> scrap >> scrap >> sInfo >> smissing1 >> smissing2;
  unsigned locus = 0;
  while(getline(ScoreTable, scrap)){
    if(!data.isTypedLocus(locus)){//skip typed loci
      //count NA as 100%

      //count percent info
      COUNTINFO(sInfo, sumPI)
      COUNTINFO(smissing1, sumMissing1)
      COUNTINFO(smissing2, sumMissing2)
    }
    ScoreTable >> scrap >> scrap >> scrap >> scrap >> sInfo >> smissing1 >> smissing2;
    ++locus;
  }
  ScoreTable.close();
  const float size = (float)data.getNumTypedLoci();

  Log << Off << "Average Information extracted across untyped loci in allelic assoc score test: " << sumPI / size << "%"
      << "\nOn average " << sumMissing1/size << "% Missing Info due to uncertainty about genotypes\n and " 
      << sumMissing2/size << "% Missing Info due to uncertainty about model parameters\n\n" ;
}

void HapMixAllelicAssocTest::Update(int locus, const double* Covariates, double phi, double YMinusEY, double DInvLink, 
				    const std::vector<std::vector<double> > Probs) {
  CopyNumberAssocTest::Update(locus, Covariates, phi, YMinusEY, DInvLink, Probs);
}
