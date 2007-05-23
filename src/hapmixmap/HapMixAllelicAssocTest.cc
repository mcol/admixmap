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
#include "Comms.h"
#include "HapMixIndividualCollection.h"
#include "bcppcl/Regression.h"

void HapMixAllelicAssocTest::Update(const HapMixIndividualCollection* const IC, const Regression* const R, const Genome& Loci){
  this->Reset();
  const double dispersion = R->getDispersion();
  const double* const EY = R->getExpectedOutcome();
  const int NumberOfIndividuals = IC->getNumberOfIndividualsForScoreTests();
  unsigned worker_rank = Comms::getWorkerRank();
  if(Comms::isWorker()){
    for( int i = worker_rank ; i < NumberOfIndividuals; i+=Comms::getNumWorkers() ){
      this->UpdateB(R->DerivativeInverseLinkFunction(i), dispersion, 0);
    }
    const int offset = IC->getFirstScoreTestIndividualNumber();
    
    for( int i = worker_rank ; i < NumberOfIndividuals; i+=Comms::getNumWorkers() ){
      
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
  }
  //TODO reduce scores
  if(Comms::isMaster()){
    this->Accumulate();
  }
}

#define COUNTINFO(S, X) \
  if(S == "NA" || S =="NaN") \
    X +=100.0;\
  else \
    X += atof(S.c_str());

void HapMixAllelicAssocTest::PrintAverageInfo(LogWriter& Log, const InputData& data, const char* filename){
  //filename is the final table written by this class earlier
  ifstream ScoreTable(filename);
  //TODO?? check this file exists. It should.

  string scrap;
  float sumPI = 0.0, sumMissing1 = 0.0, sumMissing2 = 0.0;
  string sInfo, smissing1, smissing2;
  ScoreTable >> scrap >> scrap >> scrap >> scrap >> sInfo >> smissing1 >> smissing2;
  unsigned locus = 0;
  while(getline(ScoreTable, scrap)){
    if(!data.isTypedLocus(locus)){//skip typed loci
      //count NA as 100%

      //count percent info
      COUNTINFO(sInfo, sumPI)
      COUNTINFO(smissing1, sumMissing1)
      COUNTINFO(smissing2, sumMissing2)

      ++locus;
    }
    ScoreTable >> scrap >> scrap >> scrap >> scrap >> sInfo >> smissing1 >> smissing2;
  }
  ScoreTable.close();
  const float size = (float)data.getNumTypedLoci();

  Log << Off << "Average Information extracted across untyped loci in allelic assoc score test: " << sumPI / size << "%"
      << "\nOn average " << sumMissing1/size << "% due to ?? and " << sumMissing2/size << "% due to ??\n\n" ;

}
