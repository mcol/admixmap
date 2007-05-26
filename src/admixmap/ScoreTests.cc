/** 
 *   ADMIXMAP
 *   ScoreTests.cc 
 *   Class acts as a container for the following tests:
 *   (1) Score test for allelic association
 *   (2) Score test for within-halpotype association
 *   (3) Score test for admixture association (admixturescoretest)
 *   (4) Score test for linkage with locus ancestry
 *   (5) Affecteds-only score test for linkage with locus ancestry
 *
 *   Copyright (c) 2005 - 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include <numeric>
#include "ScoreTests.h"
#include "bcppcl/Regression.h"
#include "Comms.h"

using namespace std;

ScoreTests::ScoreTests(){
  options = 0;
  individuals = 0;
  rank = Comms::getRank();
  worker_rank = Comms::getWorkerRank();
  NumWorkers = Comms::getNumWorkers();
}

ScoreTests::~ScoreTests(){

}

void ScoreTests::Initialise(AdmixOptions* op, const IndividualCollection* const indiv, const Genome* const Loci, 
			    const Vector_s& PLabels, LogWriter &Log){
  options = op;
  individuals = indiv;
  chrm = Loci->getChromosomes();
  Lociptr = Loci;
  Log.setDisplayMode(Quiet);
  if(worker_rank==-1) worker_rank =indiv->getSize();//to stop master iterating through individuals

  int K = options->getPopulations();
  int L = Lociptr->GetNumberOfCompositeLoci();

  /*----------------------
    | admixture association |
    -----------------------*/
  //TODO check conditions on this test
  if( options->getTestForAdmixtureAssociation() ){
    AdmixtureAssocScoreTest.Initialise(K, indiv->getNumberOfOutcomeVars(), options->getAssocScoreFilename(), PLabels, Log);
  }

  /*------------------------------------
    |affecteds only linkage with ancestry |
    ------------------------------------*/ 
  if( options->getTestForAffectedsOnly() ){
    AffectedsOnlyScoreTest.Initialise(options->getAffectedsOnlyScoreFilename(), K, L, Log);
  }

  /*-----------------------
    | Linkage with ancestry  |
    -----------------------*/
  if( options->getTestForLinkageWithAncestry() ){
    AncestryAssocScoreTest.Initialise(options->getAncestryAssociationScoreFilename(), K, L, Log);
  }
  
  /*----------------------
    | Allelic association  |
    -----------------------*/
  if( options->getTestForAllelicAssociation() ){
    AllelicAssociationTest.Initialise(op, indiv, Loci, Log);
  }  
}

/**
   Resets arrays holding sums of scores and info over individuals to zero.
   Invoked at start of each iteration after burnin.
*/
void ScoreTests::Reset(){
  AdmixtureAssocScoreTest.Reset();  
  AllelicAssociationTest.Reset();
  
}
void ScoreTests::SetAllelicAssociationTest(const std::vector<double> &alpha0){
  AllelicAssociationTest.SetAllelicAssociationTest(alpha0);
}

// ****************************** UPDATES ****************************

void ScoreTests::Update(const vector<Regression* >& R)
{
  Reset();//set sums over individuals to zero
  double DInvLink;
  const int NumberOfIndividuals = individuals->getNumberOfIndividualsForScoreTests();
  //------------------------------------------------------
  // Accumulate Scores over individuals for this iteration
  //------------------------------------------------------
  if(options->getNumberOfOutcomes() > 0){//if regressionmodel
    const double* const EY = R[0]->getExpectedOutcome();
    const double dispersion = R[0]->getDispersion();
    
    //NOTE: in future this loop will be outside score tests classes so the indiv indices can be controlled outside
    const int offset = individuals->getFirstScoreTestIndividualNumber();

    for( int i = worker_rank ; i < NumberOfIndividuals; i+=NumWorkers ){
      
      Individual* ind = individuals->getIndividual(i + offset);
      double YMinusEY = individuals->getOutcome(0, i) - EY[i];//individual outcome - its expectation
      //note that it is the first regression that is used
      DInvLink = R[0]->DerivativeInverseLinkFunction(i);
      bool missingOutcome = individuals->isMissingOutcome(0,i); 
     
      //admixture association
      if( options->getTestForAdmixtureAssociation() && (options->getNumberOfOutcomes() == 1) ){
	AdmixtureAssocScoreTest.UpdateIndividualScore(ind->getAdmixtureProps(), YMinusEY, dispersion, 
						      DInvLink, options->isRandomMatingModel());
      }
      //allelic association
      if( options->getTestForAllelicAssociation() )
	AllelicAssociationTest.Update( ind, YMinusEY,dispersion, DInvLink, (bool)missingOutcome);
    }
  }
  
  if(rank==0){

      //-----------------------------
      //Accumulate Scores, Info etc. over iterations
      //-----------------------------
      
      /*----------------------
	| admixture association |
	-----------------------*/
      if( options->getTestForAdmixtureAssociation() ){
	AdmixtureAssocScoreTest.Accumulate();
      }
      /*-----------------------
	| Linkage with ancestry  |
	-----------------------*/
      if( options->getTestForLinkageWithAncestry() ){
	AncestryAssocScoreTest.Accumulate();
      } 
      /*------------------------------------
	|affecteds-only linkage with ancestry |
	------------------------------------*/ 
      if( options->getTestForAffectedsOnly() ){
	AffectedsOnlyScoreTest.Accumulate();
      }
      if(options->getTestForAllelicAssociation()){
	/*-------------------------------------
	  | Allelic and haplotype association  |
	  -------------------------------------*/
	AllelicAssociationTest.Accumulate();
      }
    
  }//end if rank==0
}

// ********** OUTPUT **********************************************************

void ScoreTests::Output(const Vector_s& PLabels, const Vector_s& LocusLabels, bool final){
  //Allelic association
  if( options->getTestForAllelicAssociation() )    {
    AllelicAssociationTest.Output(LocusLabels, final);
  }//end if allelic assoc test  

  //ancestry association
  if( options->getTestForLinkageWithAncestry() ){
    const char* finalfilename = 0;
    if(final){
      string filename(options->getResultsDir());
      filename.append("/TestsAncestryAssocFinal.txt");
      finalfilename = filename.c_str();
     }
    AncestryAssocScoreTest.Output(PLabels, *Lociptr, final, finalfilename);
  }
  //affectedonly
  if( options->getTestForAffectedsOnly() ){
    const char* finalfilename = 0;
    if(final){
      string filename(options->getResultsDir());
      filename.append("/TestsAffectedsOnlyFinal.txt");
      finalfilename = filename.c_str();
    }
    AffectedsOnlyScoreTest.Output(PLabels, *Lociptr, final, finalfilename);
  }
  
  //admixture association
  if( !final && options->getTestForAdmixtureAssociation() ){
    AdmixtureAssocScoreTest.Output();
  }
}


void ScoreTests::ROutput(){
  /**
   * writes out the dimensions and labels of the 
   * R object previously written to allelicAssocScoreStream
   */
  if(options->getTestForAllelicAssociation()){
      AllelicAssociationTest.ROutput();
  }  
  
  /**
   * writes out the dimensions and labels of the 
   * R-matrix previously written to ancestryAssociationScoreStream
   */
  if (options->getTestForLinkageWithAncestry()){
    AncestryAssocScoreTest.ROutput();
  }
  
  /**
   * writes out the dimensions and labels of the 
   * R-matrix previously written to affectedsOnlyScoreStream
   */
  if (options->getTestForAffectedsOnly()){
    AffectedsOnlyScoreTest.ROutput();
  }
}

AffectedsOnlyTest& ScoreTests::getAffectedsOnlyTest(){
  return AffectedsOnlyScoreTest;
}
CopyNumberAssocTest& ScoreTests::getAncestryAssocTest(){
  return AncestryAssocScoreTest;
}

void ScoreTests::OutputLikelihoodRatios(const char* const filename, const Vector_s& PopLabels){
  AffectedsOnlyScoreTest.OutputLikRatios(filename, PopLabels, *Lociptr);
}
