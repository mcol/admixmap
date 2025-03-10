//=============================================================================
//
// Copyright (C) 2005-2007  David O'Donnell, Clive Hoggart and Paul McKeigue
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file ScoreTests.cc
/// Implementation of the ScoreTests class.
//=============================================================================

#include "ScoreTests.h"
#include "AdmixFilenames.h"
#include "bclib/LogWriter.h"

using namespace std;

ScoreTests::ScoreTests() :
  options(0),
  individuals(0) 
{
}

ScoreTests::~ScoreTests(){

}

void ScoreTests::Initialise(AdmixOptions* op, const IndividualCollection* const indiv, const Genome* const Loci, 
			    const Vector_s& PLabels, bclib::LogWriter &Log){
  options = op;
  individuals = indiv;
  Lociptr = Loci;
  Log.setDisplayMode(bclib::Quiet);
 
  int K = options->getPopulations();
  int L = Lociptr->GetNumberOfCompositeLoci();

  const string dirname = options->getResultsDir() + "/";

  /*----------------------
    | admixture association |
    -----------------------*/
  //TODO check conditions on this test
  if( options->getTestForAdmixtureAssociation() ){
    const string filename = dirname + ADMIXTUREASSOCTESTFILE;
    AdmixtureAssocScoreTest.Initialise(K, indiv->getNumberOfOutcomeVars(), filename.c_str(), PLabels, Log);
  }

  /*------------------------------------
    |affecteds only linkage with ancestry |
    ------------------------------------*/ 
  if( options->getTestForAffectedsOnly() ){
    const string filename = dirname + AFFECTEDSONLYTEST_PVALUES;
    AffectedsOnlyScoreTest.Initialise(filename.c_str(), K, L);
  }

  /*-----------------------
    | Linkage with ancestry  |
    -----------------------*/
  if( options->getTestForLinkageWithAncestry() ){
    const string filename = dirname + ANCESTRYASSOCTEST_PVALUES;
    AncestryAssocScoreTest.Initialise(filename.c_str(), K, L);
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
void ScoreTests::MergeRareHaplotypes(const std::vector<double> &alpha0){
  AllelicAssociationTest.MergeRareHaplotypes(alpha0);
}

// ****************************** UPDATES ****************************

void ScoreTests::Update(const vector<bclib::Regression* >& R)
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

    for( int i = 0; i < NumberOfIndividuals; i++ ){
      
      const PedBase & ind = individuals->getElement( i + offset );
      double YMinusEY = individuals->getOutcome(0, i) - EY[i];//individual outcome - its expectation
      //note that it is the first regression that is used
      DInvLink = R[0]->DerivativeInverseLinkFunction(i);
      bool missingOutcome = individuals->isMissingOutcome(0,i); 

      //admixture association
      if( options->getTestForAdmixtureAssociation() && (options->getNumberOfOutcomes() == 1) ){
	AdmixtureAssocScoreTest.UpdateIndividualScore(ind.getAdmixtureProps(), YMinusEY, dispersion, 
						      DInvLink, options->isRandomMatingModel());
      }
      //allelic association
      if( options->getTestForAllelicAssociation() )
        AllelicAssociationTest.Update(&ind, YMinusEY, dispersion,
                                      DInvLink, missingOutcome);
    }
  }
  
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
}

// ********** OUTPUT **********************************************************

void ScoreTests::Output(const Vector_s& PLabels){
  //Allelic association
  if( options->getTestForAllelicAssociation() )    {
    AllelicAssociationTest.Output();
  }//end if allelic assoc test  

  //ancestry association
  if( options->getTestForLinkageWithAncestry() ){
    AncestryAssocScoreTest.Output(PLabels, *Lociptr);
  }
  //affectedsonly
  if( options->getTestForAffectedsOnly() ){
    AffectedsOnlyScoreTest.Output(PLabels, *Lociptr);
  }
  
  //admixture association
  if( options->getTestForAdmixtureAssociation() ){
    AdmixtureAssocScoreTest.Output();
  }
}

void ScoreTests::WriteFinalTables(const Vector_s& PLabels, bclib::LogWriter& Log){

  const string dirname = options->getResultsDir() + "/";

  //Allelic association
  if( options->getTestForAllelicAssociation() )    {
    AllelicAssociationTest.WriteFinalTables(Log);
  }//end if allelic assoc test  

  //ancestry association
  if( options->getTestForLinkageWithAncestry() ){
    const string finalfilename = dirname + ANCESTRYASSOCTEST_FINAL;
    AncestryAssocScoreTest.WriteFinalTable(finalfilename.c_str(), PLabels, *Lociptr, Log);
  }
  //affectedsonly
  if( options->getTestForAffectedsOnly() ){
    string finalfilename = dirname + AFFECTEDSONLYTEST_FINAL;
    AffectedsOnlyScoreTest.WriteFinalTable(finalfilename.c_str(), PLabels, *Lociptr, Log);
  }
  
}

void ScoreTests::OutputLikelihoodRatios(const Vector_s& PopLabels) {
  string filename = options->getResultsDir() + "/" + AFFECTEDSONLY_LIKRATIOFILE;
  AffectedsOnlyScoreTest.OutputLikRatios(filename.c_str(), PopLabels, *Lociptr);
}
