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
#include <gsl/gsl_cdf.h>
#include "ScoreTests.h"
#include "utils/linalg.h"
#include "regression/Regression.h"
#include "Comms.h"

using namespace std;

ScoreTests::ScoreTests(){
  LocusLinkageAlleleScore = 0;
  LocusLinkageAlleleInfo = 0;
  SumLocusLinkageAlleleScore2 = 0;
  SumLocusLinkageAlleleScore = 0;
  SumLocusLinkageAlleleInfo = 0;
  locusObsIndicator = 0;

  ScoreWithinHaplotype = 0;
  InfoWithinHaplotype = 0;
  SumScoreWithinHaplotype = 0;
  SumScore2WithinHaplotype = 0;
  SumInfoWithinHaplotype = 0;

  options = 0;
  individuals = 0;
  rank = Comms::getRank();
  worker_rank = Comms::getWorkerRank();
  NumWorkers = Comms::getNumWorkers();
  dim_ = 0;
  NumOutputs = 0;
  onFirstLineAllelicAssoc = true;
  onFirstLineHapAssoc = true;
}

ScoreTests::~ScoreTests(){
  //delete arrays for allelic assoc score test
  delete[] dim_;
  delete[] locusObsIndicator;
  if(LocusLinkageAlleleScore){//proxy for options->getTestForLinkageWithAncestry()
    for(unsigned i = 0; i < Lociptr->GetNumberOfCompositeLoci(); ++i){
      //it is safe to call this function in Genome as the Loci object will be deleted after ScoreTests
      delete[] LocusLinkageAlleleScore[i];
      delete[] LocusLinkageAlleleInfo[i];
      if(rank==0){
	delete[] SumLocusLinkageAlleleScore[i];
	delete[] SumLocusLinkageAlleleScore2[i];
	delete[] SumLocusLinkageAlleleInfo[i];
      }
    }
    delete[] LocusLinkageAlleleScore;
    delete[] LocusLinkageAlleleInfo;
    if(rank==0){
      delete[] SumLocusLinkageAlleleScore;
      delete[] SumLocusLinkageAlleleScore2;
      delete[] SumLocusLinkageAlleleInfo;
    }
  }

  //delete arrays for haplotype assoc score test
  if(ScoreWithinHaplotype){
    for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); ++j){
      int NumberOfLoci = Lociptr->getNumberOfLoci(j);
      if( ScoreWithinHaplotype[j]){
	//free_matrix(ScoreWithinHaplotype[j], NumberOfLoci);
	//free_matrix(InfoWithinHaplotype[j], NumberOfLoci);

	for(int jj = 0; jj < NumberOfLoci; ++jj){
	  delete[] ScoreWithinHaplotype[j][jj];
	  delete[] InfoWithinHaplotype[j][jj];
	}
	delete[] ScoreWithinHaplotype[j];
	delete[] InfoWithinHaplotype[j];
	if(rank==0){
	  delete[] SumScoreWithinHaplotype[j];
	  delete[] SumScore2WithinHaplotype[j];
	  delete[] SumInfoWithinHaplotype[j];
	}
      }
    }
    delete[] ScoreWithinHaplotype;
    delete[] InfoWithinHaplotype;
    if(rank==0){
      delete[] SumScoreWithinHaplotype;
      delete[] SumScore2WithinHaplotype;
      delete[] SumInfoWithinHaplotype;
    }
  }
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
      NewAllelicAssocTest.Initialise(options->getAllelicAssociationScoreFilename(), 1, L, Log, false);
    }
    else{//use original allelic assoc test
      if(rank==0)OpenFile(Log, &allelicAssocScoreStream, options->getAllelicAssociationScoreFilename(), "Tests for allelic association");
      dim_ = new unsigned[L];//dimensions of arrays
      
      LocusLinkageAlleleScore = new double*[L];
      LocusLinkageAlleleInfo = new double*[L];
      ScoreWithinHaplotype = new double**[ L ];
      InfoWithinHaplotype = new double**[ L ];
      
      if(rank==0){
	SumLocusLinkageAlleleScore2 = new double*[L];
	SumLocusLinkageAlleleScore = new double*[L];
	SumLocusLinkageAlleleInfo = new double*[L];
	
	SumScoreWithinHaplotype = new double*[L];
	SumScore2WithinHaplotype = new double*[L];
	SumInfoWithinHaplotype = new double*[L];
      }
      
      //if(!options->getHapMixModelIndicator()){
#ifdef PARALLEL
      int dimalleleinfo = 0, dimhapinfo = 0;
#endif
      locusObsIndicator = new int[L];
      
      //search for loci with no observed genotypes
      for(int j = 0; j < L; ++j){
	locusObsIndicator[j] = false;
	for(int i = worker_rank + indiv->getFirstScoreTestIndividualNumber(); i < indiv->getNumberOfIndividualsForScoreTests(); i+= NumWorkers){
	  if(!indiv->getIndividual(i)->GenotypeIsMissing(j)){
	    locusObsIndicator[j] = true;
	  }
	}
      }
      //}    
      unsigned NumCovars = indiv->GetNumCovariates() - indiv->GetNumberOfInputCovariates();
      for( int j = 0; j < L; j++ ){
	int NumberOfStates = Lociptr->GetNumberOfStates(j);
	int NumberOfLoci = Lociptr->getNumberOfLoci(j);
	
	
	if(NumberOfLoci > 1 )//haplotype
	  dim_[j] = 1;
	else if(NumberOfStates == 2 )//simple diallelic locus
	  dim_[j] = 1;
	else
	  dim_[j] = NumberOfStates;//simple multiallelic locus
#ifdef PARALLEL
	dimalleleinfo += (dim_[j] + NumCovars) * (dim_[j] + NumCovars);
#endif
	
	//next two lines may not be necessary as these arrays are sized later
	LocusLinkageAlleleScore[j] = new double[ dim_[j] + NumCovars ];
	LocusLinkageAlleleInfo[j] = new double[( dim_[j] + NumCovars) * (dim_[j] + NumCovars )];
	if(rank==0){
	  SumLocusLinkageAlleleScore[j] = new double[ dim_[j] ];
	  SumLocusLinkageAlleleInfo[j] = new double[( dim_[j]) * (dim_[j] )];
	  SumLocusLinkageAlleleScore2[j] = new double[( dim_[j]) * (dim_[j] )];
	  fill(SumLocusLinkageAlleleScore[j], SumLocusLinkageAlleleScore[j]+dim_[j], 0.0);
	  fill(SumLocusLinkageAlleleInfo[j], SumLocusLinkageAlleleInfo[j] + dim_[j]*dim_[j], 0.0);
	  fill(SumLocusLinkageAlleleScore2[j], SumLocusLinkageAlleleScore2[j] + dim_[j]*dim_[j], 0.0);
	}
	
	
	if( NumberOfLoci > 1 ){
#ifdef PARALLEL
	  dimhapinfo += NumberOfLoci * (1 + NumCovars) * (1 + NumCovars);
#endif
	  ScoreWithinHaplotype[ j ] = new double*[NumberOfLoci];
	  InfoWithinHaplotype[ j ] = new double*[NumberOfLoci];
	  
	  for(int jj = 0; jj < NumberOfLoci; ++jj){
	    ScoreWithinHaplotype[j][jj] = new double[1 + NumCovars];
	    InfoWithinHaplotype[j][jj] = new double[(1 + NumCovars)*(1 + NumCovars)];
	  }
	  if(rank==0){
	    SumScoreWithinHaplotype[j] = new double[ NumberOfLoci ];
	    SumScore2WithinHaplotype[j] = new double[ NumberOfLoci ];
	    SumInfoWithinHaplotype[j] = new double[ NumberOfLoci ];
	    fill(SumScoreWithinHaplotype[j], SumScoreWithinHaplotype[j]+NumberOfLoci, 0.0);
	    fill(SumScore2WithinHaplotype[j], SumScore2WithinHaplotype[j]+NumberOfLoci, 0.0);
	    fill(SumInfoWithinHaplotype[j], SumInfoWithinHaplotype[j]+NumberOfLoci, 0.0);
	  }
	  
	}
	else{
	  ScoreWithinHaplotype[ j ] = 0;
	  InfoWithinHaplotype[ j ] = 0;
	}
      }//end loop over loci
#ifdef PARALLEL
      Comms::SetDoubleWorkspace(max(dimalleleinfo, dimhapinfo), rank==0);
      Comms::SetIntegerWorkspace(L, (rank==0));
      //accumulate locusObsIndicator over individuals (processes) 
      Comms::AllReduce_int(locusObsIndicator, L);
#endif
    }
  }  

  /*----------------------
    | haplotype association |
    ----------------------*/  
  if( strlen( options->getHaplotypeAssociationScoreFilename() ) ){
    if(Lociptr->GetTotalNumberOfLoci() > Lociptr->GetNumberOfCompositeLoci()){//cannot test for SNPs in Haplotype if only simple loci
      if(rank==0)OpenFile(Log, &HaplotypeAssocScoreStream, options->getHaplotypeAssociationScoreFilename(), "Tests for haplotype associations");
    }
    else {
      op->setTestForHaplotypeAssociation(false);
      if(rank==0)Log << "ERROR: Cannot test for haplotype associations if all loci are simple\n" << "This option will be ignored\n";
    }
  }

  /*--------------------------------
    | Residual Allelic association  |
    -------------------------------*/
  ResidualAllelicAssocScoreTest.Initialise(op, indiv, Loci, Log);

}

void ScoreTests::OpenFile(LogWriter &Log, std::ofstream* outputstream, const char* filename, std::string testname){
  outputstream->open(filename, ios::out);
  if(!outputstream->is_open()){
    string error_string = "ERROR: could not open ";
    error_string.append(filename);
    throw(error_string);
  }
  Log << testname << " written to " << filename << "\n";
  //start writing R object
  *outputstream << "structure(.Data=c(";

}

void ScoreTests::Reset(){
  //resets arrays holding sums of scores and info over individuals to zero; invoked at start of each iteration after burnin.

  if( options->getTestForAllelicAssociation() && !options->getHapMixModelIndicator()){
    int K = individuals->GetNumCovariates() - individuals->GetNumberOfInputCovariates();
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      fill(LocusLinkageAlleleScore[j], LocusLinkageAlleleScore[j]+dim_[j]+K, 0.0);
      fill(LocusLinkageAlleleInfo[j], LocusLinkageAlleleInfo[j]+(dim_[j]+K)*(dim_[j]+K), 0.0);
	    
#ifndef PARALLEL
      if( Lociptr->getNumberOfLoci(j) > 1 ){
	for(int jj = 0; jj < (* Lociptr)(j)->GetNumberOfLoci(); ++jj){
	  fill(ScoreWithinHaplotype[j][jj], ScoreWithinHaplotype[j][jj]+K+1, 0.0);
	  fill(InfoWithinHaplotype[j][jj], InfoWithinHaplotype[j][jj]+(K+1)*(K+1), 0.0);
	}
      }
#endif
    }
  }

  AdmixtureAssocScoreTest.Reset();  
  ResidualAllelicAssocScoreTest.Reset();
  NewAllelicAssocTest.Reset();
  
}

void ScoreTests::SetAllelicAssociationTest(const std::vector<double> &alpha0){
  /*
    Invokes merging of haplotypes in CompositeLocus and resizes arrays for allelic association test accordingly  
    alpha0 = alpha[0], pop admixture dirichlet params, from Latent
  */

  //first scale alphas so they sum to 1
  const unsigned K = options->getPopulations();
  double* alphaScaled = new double[K];
  double sum  = accumulate(alpha0.begin(), alpha0.end(), 0.0, std::plus<double>());//sum of alpha0 over pops
  for( int k = 0; k < options->getPopulations(); k++ )
    alphaScaled[k] = alpha0[k] / sum;

  unsigned NumCovars = individuals->GetNumCovariates() - individuals->GetNumberOfInputCovariates();
#ifdef PARALLEL
  int dimalleleinfo = 0;
#endif
  //merge rare haplotypes
  for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){//skip simple loci
      (*Lociptr)(j)->SetDefaultMergeHaplotypes( alphaScaled);
      // dim is set to 1 if single diallelic locus, to number of
      // merged haplotypes if compound locus, to number of states
      // if >2 alleles
      dim_[j] = (*Lociptr)(j)->GetNumberOfMergedHaplotypes();
      
      //resize arrays for Allelic association test
      delete[] LocusLinkageAlleleScore[j];
      delete[] LocusLinkageAlleleInfo[j];
      LocusLinkageAlleleScore[j] = new double[ dim_[j] + NumCovars];
      LocusLinkageAlleleInfo[j] = new double[( dim_[j] + NumCovars) * (dim_[j] + NumCovars )];

      if(rank==0){
	delete[] SumLocusLinkageAlleleScore[j];
	delete[] SumLocusLinkageAlleleInfo[j];
	delete[] SumLocusLinkageAlleleScore2[j];
	SumLocusLinkageAlleleScore[j] = new double[ dim_[j] ];
	SumLocusLinkageAlleleInfo[j] = new double[( dim_[j]) * (dim_[j] )];
	SumLocusLinkageAlleleScore2[j] = new double[( dim_[j]) * (dim_[j] )];
	fill(SumLocusLinkageAlleleScore[j], SumLocusLinkageAlleleScore[j]+dim_[j], 0.0);
	fill(SumLocusLinkageAlleleInfo[j], SumLocusLinkageAlleleInfo[j] + dim_[j]*dim_[j], 0.0);
	fill(SumLocusLinkageAlleleScore2[j], SumLocusLinkageAlleleScore2[j] + dim_[j]*dim_[j], 0.0);
      }

    }
#ifdef PARALLEL
    dimalleleinfo += ( dim_[j] + NumCovars) * (dim_[j] + NumCovars );
#endif
  }
#ifdef PARALLEL
  Comms::SetDoubleWorkspace(dimalleleinfo, (rank==0));
#endif

  delete[] alphaScaled;
}

// ****************************** UPDATES ****************************

void ScoreTests::Update(const vector<Regression* >& R)
{
  Reset();//set sums over individuals to zero
  double DInvLink;
  int NumberOfIndividuals = individuals->getNumberOfIndividualsForScoreTests();
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
	NewAllelicAssocTest.UpdateB(R[0]->DerivativeInverseLinkFunction(i), dispersion, 0);
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
	      NewAllelicAssocTest.Update(j, 0/*no covariates*/, dispersion, YMinusEY, DInvLink, UnorderedProbs);
	    }
	  }
	}
	else
	  UpdateScoreForAllelicAssociation( ind, YMinusEY,dispersion, DInvLink, (bool)missingOutcome);
    }
  }
  
#ifdef PARALLEL
  if( options->getTestForAllelicAssociation() ){
    Comms::ReduceAllelicAssocScores(LocusLinkageAlleleScore, LocusLinkageAlleleInfo, Lociptr->GetNumberOfCompositeLoci(), 
				    dim_,individuals->GetNumCovariates() - individuals->GetNumberOfInputCovariates());
  }

#endif

  if(rank==0){
    try{
      //-----------------------------
      //Accumulate Scores, Info etc. over iterations
      //-----------------------------
      
      /*----------------------
	| admixture association |
	-----------------------*/
      if( options->getTestForAdmixtureAssociation() ){
	AdmixtureAssocScoreTest.Update();
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
      /*------------------------------------
	|hapmixmodel allelic assoc test      |
	------------------------------------*/ 
      if(options->getHapMixModelIndicator() && options->getTestForAllelicAssociation()){
	NewAllelicAssocTest.Accumulate();
      }
      
      /*-------------------------------------
	| Allelic and haplotype association  |
	-------------------------------------*/      
      for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
	const int NumberOfLoci = Lociptr->getNumberOfLoci(j);
	
	if(!options->getHapMixModelIndicator()){
	  if( options->getTestForAllelicAssociation() ){
	    
	    //loop over simple loci within haplotypes
	    if( NumberOfLoci > 1 ){
	      for( int l = 0; l < NumberOfLoci; l++ ){
		CentreAndSum(1, ScoreWithinHaplotype[j][l], InfoWithinHaplotype[j][l], &(SumScoreWithinHaplotype[ j ][ l ]),
			     &(SumScore2WithinHaplotype[j][l]), &(SumInfoWithinHaplotype[ j ][ l ]));
	      }
	    }
	  }
	  
	  if((options->getTestForAllelicAssociation() && NumberOfLoci == 1) || options->getTestForHaplotypeAssociation()){
	    if(options->getHapMixModelIndicator() || locusObsIndicator[j]){//skip loci with no observed genotypes
	      CentreAndSum(dim_[j], LocusLinkageAlleleScore[j], LocusLinkageAlleleInfo[j],SumLocusLinkageAlleleScore[j],
			   SumLocusLinkageAlleleScore2[j],SumLocusLinkageAlleleInfo[j]); 
	    }
	  }
	}
	
      }//end comp locus loop
    }
    catch(string s){
      string error_string = "Error accumulating scores for allelicassociation or haplotype association scoretest\n";
      error_string.append(s);
      throw(error_string);
    }
    
  }//end if rank==0
}

void ScoreTests::UpdateScoreForAllelicAssociation( const Individual* const ind, double YMinusEY, double phi, double DInvLink, bool missingOutcome)
{
  int locus = 0;

  for(unsigned int j = 0; j < Lociptr->GetNumberOfChromosomes(); j++ ){
    for(unsigned int jj = 0; jj < Lociptr->GetSizeOfChromosome(j); jj++ ){
      //skip loci with missing genotypes as hap pairs have not been sampled for these
      bool condition = true;

      if(options->getHapMixModelIndicator())condition = !missingOutcome;
      else condition = !ind->GenotypeIsMissing(locus);
      //in hapmixmodel, skip individuals with observed outcome
      if(condition){
	//retrieve sampled hap pair from Individual
	const int* happair = ind->getSampledHapPair(locus);
	const unsigned numStates = Lociptr->GetNumberOfStates(locus);
	const unsigned numLoci = Lociptr->getNumberOfLoci(locus);
	
	// count alleles / haplotypes      
	vector<int> counts;
	
	// if diallelic, evaluate score for allele 2 only, otherwise evaluate score for all alleles or haplotypes
	
	// special case for SNP (counts has size 1)
	if( numStates == 2 ){
	  unsigned allele2counts = (happair[0]==1) + (happair[1]==1);
	  // (*Lociptr)(locus)->getAlleleCounts(2, happair)[0];
	  counts.push_back( allele2counts );
	  
	  // general case for simple locus (counts has size nStates)
	} else if(numLoci == 1 ){
	  for( unsigned k = 0; k < numStates; k++ ){
	    counts.push_back((*Lociptr)(locus)->getAlleleCounts(k+1, happair)[0] );
	  }
	}
	// general case for composite locus (counts has size nMergedHaplotypes)
	else {
	  //count copies of allele2
	  const vector<int> allele2Counts = (*Lociptr)(locus)->getAlleleCounts(2, happair);
	  UpdateScoreForWithinHaplotypeAssociation(ind, allele2Counts, locus, YMinusEY,phi , DInvLink);
	  
	  //update score and info for each simple locus within a compound locus
	  // 	for( int l = 0; l < (*Lociptr)(j)->GetNumberOfLoci(); l++ ){
	  // 	  vector<int> a(1, allele2Counts[0]);
	  // 	  UpdateAlleleScores(ScoreWithinHaplotype[locus][l], InfoWithinHaplotype[locus][l], ind->getAdmixtureProps(), a, 
	  // 			     YMinusEY, phi, DInvLink);
	  // 	}
	  
	  
	  if( options->getTestForHaplotypeAssociation() ){
	    //count numbers of each haplotype
	    counts = (*Lociptr)(locus)->getHaplotypeCounts(happair);
	  }
	}
	UpdateAlleleScores(LocusLinkageAlleleScore[locus],LocusLinkageAlleleInfo[locus], ind->getAdmixtureProps(), counts, 
			   YMinusEY, phi, DInvLink);
      }
      locus++;
    }
  }
}

// This function calculates score for allelic association at each simple locus within a compound locus
void ScoreTests::UpdateScoreForWithinHaplotypeAssociation( const Individual* const ind, const vector<int> allele2Counts, 
							   int j, double YMinusEY, double phi, double DInvLink)
{
  int K = individuals->GetNumCovariates() - individuals->GetNumberOfInputCovariates();
  double* x = new double[ K + 1 ];


  x[ K ] = 1.0;
  for( int k = 0; k < K - 1; k++ )
    x[ k + 1 ] = ind->getAdmixtureProps()[k];//?? should be [k+1]

  for( int l = 0; l < (*Lociptr)(j)->GetNumberOfLoci(); l++ ){//loop over simple loci within compound locus
    x[0] = (double)allele2Counts[l];

    for( int k = 0; k < K + 1; k++ ){
      ScoreWithinHaplotype[j][l][ k ] += phi * x[k] * YMinusEY;
      for( int kk = 0; kk < K + 1; kk++ )
	InfoWithinHaplotype[j][l][ k*(K+1) + kk ] += phi*DInvLink * x[ k ] * x[ kk ];
    }
  }
  delete[] x;
}
void ScoreTests::UpdateAlleleScores( double* score, double* info, const double* admixtureProps, const vector<int> Counts, 
				     double YMinusEY, double phi, double DInvLink)
{
  int K = individuals->GetNumCovariates() - individuals->GetNumberOfInputCovariates();
  unsigned dim = Counts.size();
  double* x = new double[ K + dim ];

  // ** Set x co-ordinate for regression parameter under test
  copy(Counts.begin(), Counts.end(), x);
 
  // ** Set x-co-ordinates of covariates in model
  x[ dim ] = 1.0;//intercept
  for( int k = 1; k < K; k++ )
    x[ dim+k ] = admixtureProps[k];
  //if(!options->getHapMixModelIndicator())
  //copy(admixtureProps+1, admixtureProps+K, x+dim+1);

  //accumulate score and info
  for( unsigned k = 0; k < K + dim; k++ ){
    score[ k ] += phi * x[k] * YMinusEY;
    for( unsigned kk = 0; kk < K + dim; kk++ )
      info[ k*(K+dim) + kk ] += x[ k ] * x[ kk ] * phi*DInvLink;
  }
  delete[] x;
}


// corrects for covariance between score and regression parameters
// then accumulates sumscore, sumscoresq and suminfo
void ScoreTests::CentreAndSum(unsigned dim, double *score, double* info, 
			      double *sumscore, double* sumscoresq, double* suminfo)
{
  double *cscore = new double[dim];//centred score
  double *cinfo = new double[dim*dim];//centred info

  unsigned NumCovars = individuals->GetNumCovariates() - individuals->GetNumberOfInputCovariates();
  CentredGaussianConditional( dim, score, info, cscore, cinfo, NumCovars+dim );
  for(unsigned d = 0; d < dim; ++d){
    *(sumscore + d) += cscore[d];
    for(unsigned dd = 0; dd < dim; ++dd){
      *(sumscoresq + d*dim +dd) += cscore[d] * cscore[dd];
      *(suminfo + d*dim +dd) += cinfo[d*dim+dd];
    }
  }
  delete[] cscore;
  delete[] cinfo;
}

void ScoreTests::UpdateScoresForResidualAllelicAssociation(const FreqArray& AlleleFreqs){
  ResidualAllelicAssocScoreTest.Update(AlleleFreqs, options->getHapMixModelIndicator());
}


// }
// ********** OUTPUT **********************************************************

void ScoreTests::Output(int iterations, const Vector_s& PLabels, const Vector_s& LocusLabels, bool final){
  //PopLabels = PLabels;
  if(!final)++NumOutputs;
  string sep = final ? "\t" : ",";//separator
  ofstream* outfile;

  //Allelic association
  if( options->getTestForAllelicAssociation() )    {
    if(options->getHapMixModelIndicator()){
      string filename; //ignored if !final
      if(final){
	filename = (options->getResultsDir());
	filename.append("/AllelicAssocTestsFinal.txt");
      }
      else outfile = &allelicAssocScoreStream;
      NewAllelicAssocTest.Output(PLabels, *Lociptr, final, filename.c_str());
    }
    else{
    if(final){
      string filename(options->getResultsDir());
      filename.append("/AllelicAssocTestsFinal.txt");
      outfile = new ofstream(filename.c_str(), ios::out);
      *outfile << "Locus\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tStdNormal\tPValue\tChiSquare\n";
    }
    else outfile = &allelicAssocScoreStream;
    for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ )
      if(options->getHapMixModelIndicator() || locusObsIndicator[j]){
	const int NumberOfLoci = Lociptr->getNumberOfLoci(j);
	const int NumberOfStates = Lociptr->GetNumberOfStates(j);
	const string locuslabel = LocusLabels[j];

	//case of simple locus
	if( NumberOfLoci == 1 ){
	  if(NumberOfStates==2)
	    {//SNP
	    OutputScalarScoreTest(iterations, outfile, locuslabel,
				  SumLocusLinkageAlleleScore[j][0], SumLocusLinkageAlleleScore2[j][0], 
				  SumLocusLinkageAlleleInfo[j][0], final, onFirstLineAllelicAssoc);
	  }
	  else {//multiallelic simple locus
	    vector<string> labels;
	    for(unsigned i = 0; i < dim_[j]; ++i){
	      stringstream ss;
	      ss << "\"" << locuslabel << "(" << i+1 << ")\"";
	      labels.push_back(ss.str());
	    }
	    OutputScoreTest(iterations, outfile, dim_[j], labels, 
			    SumLocusLinkageAlleleScore[j], SumLocusLinkageAlleleScore2[j], 
			    SumLocusLinkageAlleleInfo[j], final, onFirstLineAllelicAssoc, dim_[j]);
	  }
	}//end if simple locus
	//case of haplotype
	else{
	  //for(int i = 0; i < (*Lociptr)(j)->GetNumberOfLoci(); ++i)labels.push_back("\""+(*Lociptr)(j)->GetLabel(i)+"\"");
	  for(int simplelocus = 0; simplelocus < NumberOfLoci; ++simplelocus){
	    string simplelocuslabel = (*Lociptr)(j)->GetLabel(simplelocus);
	    OutputScalarScoreTest(iterations, outfile, simplelocuslabel, SumScoreWithinHaplotype[ j ][simplelocus], 
				  SumScore2WithinHaplotype[ j ][simplelocus], SumInfoWithinHaplotype[ j ][simplelocus], final, onFirstLineAllelicAssoc);
	  }
	}
	onFirstLineAllelicAssoc = false;
      }//end j loop over comp loci
    if(final)delete outfile;
    }//end if not hapmixmodel

  }//end if allelic assoc test  
  //haplotype association
  if( options->getTestForHaplotypeAssociation() ){
    if(final){
      string filename(options->getResultsDir());
      filename.append("/HaplotypeAssocTestsFinal.txt");
      outfile = new ofstream(filename.c_str(), ios::out);
      *outfile << "Locus\tHaplotype\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tStdNormal\tPValue\tChiSquare\n";
    }
    else outfile = &HaplotypeAssocScoreStream;
    vector<string> labels;
    const int *hap = 0;
    for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ) if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
      //here, dim_[j] = NumberOfMergedHaplotypes
      //create labels as "locuslabel","haplabel"
      labels.clear();
      for(unsigned i = 0; i < dim_[j]; ++i){
	stringstream ss;
	ss  << "\"" << (*Lociptr)(j)->GetLabel(0) << "\""<< sep;
	if( i < dim_[j] - 1 ){
	  hap = (*Lociptr)(j)->GetHapLabels(i);
	  ss  << "\"";
	  for( int kk = 0; kk < (*Lociptr)(j)->GetNumberOfLoci() - 1; kk++ ){
	    ss  << hap[kk] << "-";
	  }
	  ss  << hap[(*Lociptr)(j)->GetNumberOfLoci() - 1] << "\"";
	}
	else
	  ss  << "\"others\"";
	labels.push_back(ss.str());
      }
      OutputScoreTest(iterations, outfile, dim_[j], labels, SumLocusLinkageAlleleScore[j], 
		      SumLocusLinkageAlleleScore2[j], SumLocusLinkageAlleleInfo[j], final, onFirstLineHapAssoc, dim_[j]-1);
      onFirstLineHapAssoc = false;
    }
    if(final)delete outfile;
  }//end if hap assoc test
  
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
  ResidualAllelicAssocScoreTest.Output(iterations, final, LocusLabels);
  
  //admixture association
  if( !final && options->getTestForAdmixtureAssociation() ){
    AdmixtureAssocScoreTest.Output(iterations);
  }
}

// next few function calculate score tests from the cumulative sums of
// the score, score squared, and information score and info can be
// scalars, or respectively a vector and a matrix

// generic scalar score test
//TODO: move output of NA in chisq column outside as it is only required if along with vector tests
void ScoreTests::OutputScalarScoreTest( int iterations, ofstream* outputstream, string label, const double score, const double scoresq, const double info, bool final, bool firstline)
{
  string sep = final? "\t" : ",";
  double Score = score / (double) iterations;
  double CompleteInfo = info / (double) iterations;
  double MissingInfo = scoresq / ( double)iterations - Score * Score;
  double ObservedInfo = CompleteInfo - MissingInfo;

  if(!firstline){
    //if not the first line, output a separator after previous output
    *outputstream << sep;
  }
  //now start a new line (for ease of human reading)
  *outputstream << endl;


  //output label
  *outputstream << "\"" << label << "\"" << sep;
  if(final)
    *outputstream << double2R(Score, 3)        << sep
		  << double2R(CompleteInfo, 3) << sep
		  << double2R(ObservedInfo, 3) << sep;

  if( (MissingInfo < CompleteInfo) && (CompleteInfo > 0.0) ) {
    double PercentInfo = 100.0*ObservedInfo / CompleteInfo;
    double zscore = Score / sqrt( ObservedInfo );
    double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
    if(final)
      *outputstream << double2R(PercentInfo, 2) << sep
		    << double2R(zscore,3)   << sep 
		    << double2R(pvalue);
    else
      *outputstream << double2R(-log10(pvalue));
  }
  else{
    if(final)*outputstream << "NA" << sep << "NA" << sep;
    *outputstream << "NA" ;
  }
  if(final)*outputstream << sep << "NA";//NA in chisquare column in final table 
  //*outputstream << endl;
}

//generic vector score test
void ScoreTests::OutputScoreTest( int iterations, ofstream* outputstream, unsigned dim, vector<string> labels,
				  const double* score, const double* scoresq, const double* info, bool final, bool firstline, unsigned dim2)
{
  //given cumulative scores, square of scores and info, of dimension dim, over iterations, computes expectation of score, complete info and observed info and outputs to output stream along with a summary chi-square statistic and p-value. Also performs scalar test for each element.
  //if final=false, only the log(-pvalue)'s are printed

  double *ScoreVector = 0, *CompleteInfo = 0, *ObservedInfo = 0;
  string sep = final? "\t" : ",";

  if(!firstline){
    //if not the first line, output a separator after previous output
    *outputstream << sep;
  }
  //now start a new line (for ease of human reading)
  *outputstream << endl;

  
  ScoreVector = new double[dim];
  copy(score, score+dim, ScoreVector);
  scale_matrix(ScoreVector, 1.0/( iterations), dim, 1);
  
  CompleteInfo = new double[dim*dim];
  copy(info, info + dim*dim, CompleteInfo);
  scale_matrix(CompleteInfo, 1.0/(double) iterations, dim, dim);
  
  ObservedInfo = new double[dim*dim];
  for(unsigned d1 = 0; d1 < dim; ++d1)for(unsigned d2 = 0; d2 < dim; ++d2)
    ObservedInfo[d1*dim + d2] = CompleteInfo[d1*dim+d2] + ScoreVector[d1]*ScoreVector[d2] -
      scoresq[d1*dim+d2]/(double) iterations;
  for( unsigned k = 0; k < dim; k++ ){
    if(k > 0)
      *outputstream << sep << endl;
    // ** output labels
    *outputstream << labels[k] << sep;
    if(final){
      *outputstream  << double2R(ScoreVector[k], 3) << sep
		     << double2R(CompleteInfo[k*dim+k], 3) << sep//prints diagonal of CI matrix
		     << double2R(ObservedInfo[k*dim+k], 3) << sep;//   "      "     "  MI   "
      if(CompleteInfo[k*dim+k]>0.0)
	*outputstream<< double2R(100.0*ObservedInfo[k*dim+k] / CompleteInfo[k*dim+k], 2) << sep;//%Observed Info
      else 
	*outputstream << "NA" << sep;
    }
    double zscore = ScoreVector[ k ] / sqrt( ObservedInfo[k*dim+k] );
    if(final)*outputstream  << double2R(zscore, 3) << sep;//z-score
    double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
    if(final)*outputstream << double2R(pvalue) << sep;
    else *outputstream << double2R(-log10(pvalue));
    // if not last allele at locus, output unquoted "NA" in chi-square column
    if( final && k != dim - 1 ){
      *outputstream  << "NA" ;
    }
  }//end loop over alleles
  if(final){
    double chisq=0.0;
    try{
      if(dim2==dim) chisq = GaussianQuadraticForm(ScoreVector, ObservedInfo, dim);
      else chisq = GaussianMarginalQuadraticForm( dim2, ScoreVector, ObservedInfo, dim );//marginalise over first dim2 elements
      if(chisq < 0.0)
	*outputstream << "NA" ;
      else *outputstream << double2R(chisq) ;
    }
    catch(...){//in case ObservedInfo is rank deficient
      *outputstream  << "NA";
    }
  }
  //TODO:?? output p-value for chisq
	
  delete[] ScoreVector;
  delete[] CompleteInfo;
  delete[] ObservedInfo;
}

void ScoreTests::ROutput(){
  const int numPrintedIterations = NumOutputs;//(options->getTotalSamples() - options->getBurnIn()) / (options->getSampleEvery() * 10);
  /**
   * writes out the dimensions and labels of the 
   * R object previously written to allelicAssocScoreStream
   */
  int count;
  if(options->getTestForAllelicAssociation()){
    if(options->getHapMixModelIndicator()){
      NewAllelicAssocTest.ROutput();
    }
    else{

    count = 0;
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ )if(options->getHapMixModelIndicator() || locusObsIndicator[j]){
#ifdef PARALLEL
      const int NumberOfLoci = 1;
#else
      const int NumberOfLoci = (* Lociptr)(j)->GetNumberOfLoci();
#endif
      if( NumberOfLoci == 1 ) count += dim_[j];
      else count += NumberOfLoci;
    }         
    vector<int> dimensions(3,0);
    dimensions[0] = 2;
    dimensions[1] = count;
    dimensions[2] = (int)(numPrintedIterations);
    vector<string> labels(dimensions[0],"");
    labels[0] = "Locus";
    labels[1] = "log10Pvalue";
    R_output3DarrayDimensions(&allelicAssocScoreStream,dimensions,labels);
  }
  }  
  /** 
   * writes out the dimensions and labels of the        
   * R-matrix for score tests of SNPs in haplotypes.        
   */ 
  if( options->getTestForHaplotypeAssociation()  ){      
    count = 0;
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
	count += (*Lociptr)(j)->GetNumberOfMergedHaplotypes();
      }
    }         
    vector<int> dimensions(3,0);
    dimensions[0] = 3;
    dimensions[1] = count;
    dimensions[2] = (int)(numPrintedIterations);
    vector<string> labels(dimensions[0],"");
    labels[0] = "Locus";
    labels[1] = "Haplotype";
    labels[2] = "log10PValue";
    R_output3DarrayDimensions(&HaplotypeAssocScoreStream,dimensions,labels);
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

//finishes writing scoretest output as R object
void ScoreTests::R_output3DarrayDimensions(ofstream* stream, const vector<int> dim, const vector<string> labels)
{
  *stream << endl << ")," << endl;
  *stream << ".Dim = c(";
  for(unsigned int i=0;i<dim.size();i++){
    *stream << dim[i];
    if(i != dim.size() - 1){
      *stream << ",";
    }
  }
  *stream << ")," << endl;
  *stream << ".Dimnames=list(c(";
  for(unsigned int i=0;i<labels.size();i++){
    *stream << "\"" << labels[i] << "\"";
    if(i != labels.size() - 1){
      *stream << ",";
    }
  }
  *stream << "), character(0), character(0)))" << endl;
}

//converts a double to a string for R to read
//useful only for converting infs and nans to "NaN"
string ScoreTests::double2R( double x )
{
  if( isnan(x) || isinf(x) )
    return "NaN";
  else{
    stringstream ret;
    ret << x;
    return( ret.str() );
  }
}
string ScoreTests::double2R( double x, int precision )
{
  if( isnan(x) || isinf(x) )
    return "NaN";
  else{
    stringstream ret;

    if(x < numeric_limits<float>::max( ))//in case x too large to represent as floating point number
      ret << setiosflags(ios::fixed) << setprecision(precision);
    ret << x;
    return( ret.str() );
  }
}
AffectedsOnlyTest& ScoreTests::getAffectedsOnlyTest(){
  return AffectedsOnlyScoreTest;
}
CopyNumberAssocTest& ScoreTests::getAncestryAssocTest(){
  return AncestryAssocScoreTest;
}

void ScoreTests::OutputLikelihoodRatios(const char* const filename, int iterations, const Vector_s& PopLabels){
  AffectedsOnlyScoreTest.OutputLikRatios(filename, iterations, PopLabels, *Lociptr);
}
