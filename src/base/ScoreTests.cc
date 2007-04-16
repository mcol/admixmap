/** 
 *   ScoreTests.cc 
 *   Class implements the following score tests:

 *   (1) Score test for allelic association
 *   (2) Score test for within-halpotype association
 *
 *  Also acts as a wrapper for the following tests:
 *   (3) Score test for admixture association (admixturescoretest)
 *   (4) Score test for linkage with locus ancestry
 *   (5) Affecteds-only score test for linkage with locus ancestry
 *   (6) Score test for residual allelic association between adjacent pairs of linked loci
 *   (7) hapmix allelic assoc test (using CopyNumberAssocTest)
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

void ScoreTests::Initialise(Options* op, const IndividualCollection* const indiv, const Genome* const Loci, 
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
    if(options->getHapMixModelIndicator()){//use new test, with conditional distribution
      HapMixAllelicAssocTest.Initialise(options->getAllelicAssociationScoreFilename(), 1, L, Log, false);
    }
    else{//use original allelic assoc test
      AllelicAssociationTest.Initialise(op, indiv, Loci, Log);
    }
  }  

  /*--------------------------------
    | Residual Allelic association  |
    -------------------------------*/
  ResidualAllelicAssocScoreTest.Initialise(op, indiv, Loci, Log);

}

/**
   Resets arrays holding sums of scores and info over individuals to zero.
   Invoked at start of each iteration after burnin.
*/
void ScoreTests::Reset(){
  AdmixtureAssocScoreTest.Reset();  
  AllelicAssociationTest.Reset();
  ResidualAllelicAssocScoreTest.Reset();
  HapMixAllelicAssocTest.Reset();
  
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
    //eg if(hapmixmodel && casecontrolanalysis && i>NumIndividuals)update allelic assoc test
    const int offset = individuals->getFirstScoreTestIndividualNumber();

    if(options->getHapMixModelIndicator() && options->getTestForAllelicAssociation())
      for( int i = worker_rank ; i < NumberOfIndividuals; i+=NumWorkers ){
	HapMixAllelicAssocTest.UpdateB(R[0]->DerivativeInverseLinkFunction(i), dispersion, 0);
      }

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
	if(options->getHapMixModelIndicator()){
	  // vector<double>OrderedProbs(4);
	  vector<vector<double> > UnorderedProbs(3, vector<double>(1));
	  
	  unsigned int numberCompositeLoci = Lociptr->GetNumberOfCompositeLoci();
	  for(unsigned int j = 0; j < numberCompositeLoci; j++ ) {
	    // SNPs only
	    if (Lociptr->GetNumberOfStates(j) == 2) {
	      //        ind->calculateUnorderedProbs(j);
	      UnorderedProbs = ind->getUnorderedProbs(j); 
	      // from this line onward, code does not change - rest of score test update remains same
	      HapMixAllelicAssocTest.Update(j, 0/*no covariates*/, dispersion, YMinusEY, DInvLink, UnorderedProbs);
	    }
	  }
	}
	else
	  AllelicAssociationTest.Update( ind, YMinusEY,dispersion, DInvLink, (bool)missingOutcome);
    }
  }
  
#ifdef PARALLEL
  if( options->getTestForAllelicAssociation() ){
    Comms::ReduceAllelicAssocScores(LocusLinkageAlleleScore, LocusLinkageAlleleInfo, Lociptr->GetNumberOfCompositeLoci(), 
				    dim_,individuals->GetNumCovariates() - individuals->GetNumberOfInputCovariates());
  }

#endif

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
	/*------------------------------------
	  |hapmixmodel allelic assoc test      |
	  ------------------------------------*/ 
	if(options->getHapMixModelIndicator()){
	  HapMixAllelicAssocTest.Accumulate();
	}
      
	/*-------------------------------------
	  | Allelic and haplotype association  |
	  -------------------------------------*/
	else{
	  AllelicAssociationTest.Accumulate();
	}
      }
    
  }//end if rank==0
}

void ScoreTests::UpdateScoresForResidualAllelicAssociation(const FreqArray& AlleleFreqs){
  ResidualAllelicAssocScoreTest.Update(AlleleFreqs, options->getHapMixModelIndicator());
}


// }
// ********** OUTPUT **********************************************************

void ScoreTests::Output(const Vector_s& PLabels, const Vector_s& LocusLabels, bool final){
  //PopLabels = PLabels;

  //Allelic association
  if( options->getTestForAllelicAssociation() )    {
    if(options->getHapMixModelIndicator()){
      string filename; //ignored if !final
      if(final){
	filename = (options->getResultsDir());
	filename.append("/AllelicAssocTestsFinal.txt");
      }
      //else outfile = &allelicAssocScoreStream;
      HapMixAllelicAssocTest.Output(PLabels, *Lociptr, final, filename.c_str());
    }//end hapmixmodel test
    else{
      AllelicAssociationTest.Output(LocusLabels, final);
    }//end if not hapmixmodel

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
  
  //residual allelic association
  ResidualAllelicAssocScoreTest.Output(final, LocusLabels);
  
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
    if(options->getHapMixModelIndicator()){
      HapMixAllelicAssocTest.ROutput();
    }
    else{
      AllelicAssociationTest.ROutput();

    }
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

  /**
   * writes out the dimensions and labels of the 
   * R-matrix previously written to ResAlleleScoreFile
   */
  ResidualAllelicAssocScoreTest.ROutput();

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
