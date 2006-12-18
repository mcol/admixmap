/** 
 *   ADMIXMAP
 *   Scoretests.cc 
 *   Class implements the following score tests:
 *   (1) Score test for admixture association (admixturescoretest)
 *   (2) Score test for allelic association
 *   (3) Score test for within-halpotype association
 *   (4) Score test for linkage with locus ancestry
 *   (5) Affecteds-only score test for linkage with locus ancestry
 *   Copyright (c) 2005, 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "ScoreTests.h"
#include <numeric>
#include "utils/linalg.h"
#include "regression/Regression.h"
#include "Comms.h"

using namespace std;

ScoreTests::ScoreTests(){
  SumAncestryScore = 0;
  SumAncestryInfo = 0;
  SumAncestryVarScore = 0;
  SumAncestryScore2 = 0;

//   SumAffectedsScore = 0;
//   SumAffectedsScore2 = 0;
//   SumAffectedsInfo = 0;
//   SumAffectedsVarScore = 0;

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
}

ScoreTests::~ScoreTests(){
  //delete arrays for ancestry assoc score test
  delete[] SumAncestryScore;
  delete[] SumAncestryInfo;
  delete[] SumAncestryVarScore;
  delete[] SumAncestryScore2;


//   //delete arrays for affecteds-only score test
//   delete[] SumAffectedsScore;
//   delete[] SumAffectedsScore2;
//   delete[] SumAffectedsInfo;
//   delete[] SumAffectedsVarScore;

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
 //    OpenFile(Log, &affectedsOnlyScoreStream, options->getAffectedsOnlyScoreFilename(), "Affected-only tests for association");
//     //KK is the number of populations for which to perform test. For 2way admixture, we only want 2nd population.
//     int KK = K;
//     if(K == 2)KK = 1;
    
//     SumAffectedsScore2 = new double[L * KK];
//     SumAffectedsScore = new double[L * KK];
//     SumAffectedsVarScore = new double[L * KK];
//     SumAffectedsInfo = new double[L * KK];
//     fill(SumAffectedsScore, SumAffectedsScore +L*KK, 0.0);
//     fill(SumAffectedsScore2, SumAffectedsScore2 +L*KK, 0.0);
//     fill(SumAffectedsInfo, SumAffectedsInfo + L*KK, 0.0);
//     fill(SumAffectedsVarScore, SumAffectedsVarScore + L*KK, 0.0);
    AffectedsOnlyScoreTest.Initialise(options->getAffectedsOnlyScoreFilename(), K, L, Log);
  }

  /*-----------------------
    | Linkage with ancestry  |
    -----------------------*/
  if( options->getTestForLinkageWithAncestry() ){
    OpenFile(Log, &ancestryAssociationScoreStream, options->getAncestryAssociationScoreFilename(), "Tests for locus linkage");
	
    int KK = K;
    if(K == 2)KK = 1;//only keeping scores for second population when 2 populations
    SumAncestryScore = new double[L * KK];
    SumAncestryInfo = new double[L * KK];
    SumAncestryScore2 = new double[L * KK];
    SumAncestryVarScore = new double[L * KK];
    fill(SumAncestryScore, SumAncestryScore +L*KK, 0.0);
    fill(SumAncestryScore2, SumAncestryScore2 +L*KK, 0.0);
    fill(SumAncestryInfo, SumAncestryInfo + L*KK, 0.0);
    fill(SumAncestryVarScore, SumAncestryVarScore + L*KK, 0.0);
  }
  
  /*----------------------
    | Allelic association  |
    -----------------------*/
  if( options->getTestForAllelicAssociation() ){
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
  *outputstream << "structure(.Data=c(" << endl;

}

void ScoreTests::Reset(){
  //resets arrays holding sums of scores and info over individuals to zero; invoked at start of each iteration after burnin.

  if( options->getTestForAllelicAssociation() ){
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
    //-----------------------------
    //Accumulate Scores, Info etc. over iterations
    //-----------------------------
  
    /*----------------------
      | admixture association |
      -----------------------*/
    if( options->getTestForAdmixtureAssociation() ){
      AdmixtureAssocScoreTest.Update();
    }
  
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      const int NumberOfLoci = Lociptr->getNumberOfLoci(j);

      /*-------------------------------------
	| Allelic and haplotype association  |
	-------------------------------------*/
      try{
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
      catch(string s){
	string error_string = "Error accumulating scores for allelicassociation or haplotype association scoretest\n";
	error_string.append(s);
	throw(error_string);
      }
    
      /*-----------------------
	| Linkage with ancestry  |
	-----------------------*/
      if( options->getTestForLinkageWithAncestry() ){
	AdmixedIndividual::SumScoresForAncestry(j, SumAncestryScore, SumAncestryInfo, SumAncestryScore2, SumAncestryVarScore);
      } 
      /*------------------------------------
	|affecteds-only linkage with ancestry |
	------------------------------------*/ 
      if( options->getTestForAffectedsOnly() ){
	//Individual::SumScoresForLinkageAffectedsOnly(j, SumAffectedsScore, SumAffectedsVarScore, SumAffectedsScore2, SumAffectedsInfo);
	AffectedsOnlyScoreTest.Accumulate(j);
      }
    }//end comp locus loop
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
#ifdef PARALLEL
	const unsigned numStates = 2;
	const unsigned numLoci = 1;
#else
	const unsigned numStates = (*Lociptr)(locus)->GetNumberOfStates();
	const unsigned numLoci = (*Lociptr)(locus)->GetNumberOfLoci();
#endif
	
	// count alleles / haplotypes      
	vector<int> counts;
	
	// if diallelic, evaluate score for allele 2 only, otherwise evaluate score for all alleles or haplotypes
	
	// special case for SNP (counts has size 1)
	if( numStates == 2 ){
	  //sets counts to -1, 0 or 1 according to whether happair is 00, 01 or 11 (ie alleles 11, 12 or 22)
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

void ScoreTests::UpdateScoresForResidualAllelicAssociation(const array_of_allelefreqs& AlleleFreqs){
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
				  SumLocusLinkageAlleleInfo[j][0], final);
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
			    SumLocusLinkageAlleleInfo[j], final, dim_[j]);
	  }
	}//end if simple locus
	//case of haplotype
	else{
	  //for(int i = 0; i < (*Lociptr)(j)->GetNumberOfLoci(); ++i)labels.push_back("\""+(*Lociptr)(j)->GetLabel(i)+"\"");
	  for(int simplelocus = 0; simplelocus < NumberOfLoci; ++simplelocus){
	    string simplelocuslabel = (*Lociptr)(j)->GetLabel(simplelocus);
	    OutputScalarScoreTest(iterations, outfile, simplelocuslabel, SumScoreWithinHaplotype[ j ][simplelocus], 
				  SumScore2WithinHaplotype[ j ][simplelocus], SumInfoWithinHaplotype[ j ][simplelocus], final);
	  }
	}
      }//end j loop over comp loci
    if(final)delete outfile;
  }
  
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
		      SumLocusLinkageAlleleScore2[j], SumLocusLinkageAlleleInfo[j], final, dim_[j]-1);
    }
    if(final)delete outfile;
  }
  
  //ancestry association
  if( options->getTestForLinkageWithAncestry() ){
    if(final){
      string filename(options->getResultsDir());
      filename.append("/TestsAncestryAssocFinal.txt");
      outfile = new ofstream(filename.c_str(), ios::out);
      *outfile <<"Locus\tPopulation\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tMissing1\tMissing2\tStdNormal\tPValue\n";
    }
    else outfile = &ancestryAssociationScoreStream;
    OutputTestsForLocusLinkage( iterations, outfile, PLabels,
				SumAncestryScore, SumAncestryVarScore,
				SumAncestryScore2, SumAncestryInfo, sep );
    if(final)delete outfile;
  }
  //affectedonly
  if( options->getTestForAffectedsOnly() ){
    const char* finalfilename = 0;
     if(final){
       string filename(options->getResultsDir());
       filename.append("/TestsAffectedsOnlyFinal.txt");
       finalfilename = filename.c_str();
       //       outfile = new ofstream(filename.c_str(), ios::out);
       //*outfile <<"Locus\tPopulation\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tMissing1\tMissing2\tStdNormal\tPValue\n";
     }//else outfile = &affectedsOnlyScoreStream;
//     OutputTestsForLocusLinkage( iterations, outfile, PLabels, 
// 				SumAffectedsScore, SumAffectedsVarScore,
// 				SumAffectedsScore2, SumAffectedsInfo, sep );
    //if(final)delete outfile;
    AffectedsOnlyScoreTest.Output(iterations, PLabels, *Lociptr, final, finalfilename);
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
void ScoreTests::OutputScalarScoreTest( int iterations, ofstream* outputstream, string label, const double score, const double scoresq, const double info, bool final)
{
  string sep = final? "\t" : ",";
  double Score = score / ( iterations );
  double CompleteInfo = info / ( iterations );
  double MissingInfo = scoresq / ( iterations ) - Score * Score;
  double ObservedInfo = CompleteInfo - MissingInfo;
  //output label
  *outputstream << "\"" << label << "\"" << sep;
  if(final)
    *outputstream << double2R(Score, 3)        << sep
		  << double2R(CompleteInfo, 3) << sep
		  << double2R(ObservedInfo, 3) << sep;
  if( (MissingInfo < CompleteInfo) && (CompleteInfo > 0.0) ) {
    double PercentInfo = 100*ObservedInfo / CompleteInfo;
    double zscore = Score / sqrt( ObservedInfo );
    double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
    if(final)
      *outputstream << double2R(PercentInfo, 2) << sep
		    << double2R(zscore,3)   << sep 
		    << double2R(pvalue) << sep;// << endl;
    else
      *outputstream << double2R(-log10(pvalue)) << sep;// << endl;
  }
  else{
    if(final)*outputstream << "NaN" << sep << "NaN" << sep;
    *outputstream << "NaN" << sep;// << endl;
  }
  if(final)*outputstream << "NA";//NA in chisquare column in final table 
  *outputstream << endl;
}

//generic vector score test
void ScoreTests::OutputScoreTest( int iterations, ofstream* outputstream, unsigned dim, vector<string> labels,
				  const double* score, const double* scoresq, const double* info, bool final, unsigned dim2)
{
  //given cumulative scores, square of scores and info, of dimension dim, over iterations, computes expectation of score, complete info and observed info and outputs to output stream along with a summary chi-square statistic and p-value. Also performs scalar test for each element.
  //if final=false, only the log(-pvalue)'s are printed

  double *ScoreVector = 0, *CompleteInfo = 0, *ObservedInfo = 0;
  string sep = final? "\t" : ",";
  
  ScoreVector = new double[dim];
  copy(score, score+dim, ScoreVector);
  scale_matrix(ScoreVector, 1.0/( iterations), dim, 1);
  
  CompleteInfo = new double[dim*dim];
  copy(info, info + dim*dim, CompleteInfo);
  scale_matrix(CompleteInfo, 1.0/( iterations), dim, dim);
  
  ObservedInfo = new double[dim*dim];
  for(unsigned d1 = 0; d1 < dim; ++d1)for(unsigned d2 = 0; d2 < dim; ++d2)
    ObservedInfo[d1*dim + d2] = CompleteInfo[d1*dim+d2] + ScoreVector[d1]*ScoreVector[d2] -
      scoresq[d1*dim+d2]/( iterations );
  for( unsigned k = 0; k < dim; k++ ){
    // ** output labels
    *outputstream << labels[k] << sep;
    if(final){
      *outputstream  << double2R(ScoreVector[k], 3) << sep
		     << double2R(CompleteInfo[k*dim+k], 3) << sep//prints diagonal of CI matrix
		     << double2R(ObservedInfo[k*dim+k], 3) << sep;//   "      "     "  MI   "
      if(CompleteInfo[k*dim+k]>0.0)
	*outputstream<< double2R(100*ObservedInfo[k*dim+k] / CompleteInfo[k*dim+k], 2) << sep;//%Observed Info
      else 
	*outputstream << "NA" << sep;
    }
    double zscore = ScoreVector[ k ] / sqrt( ObservedInfo[k*dim+k] );
    if(final)*outputstream  << double2R(zscore, 3) << sep;//z-score
    double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
    if(final)*outputstream << double2R(pvalue) << sep;
    else *outputstream << double2R(-log10(pvalue)) << sep << endl;
    // if not last allele at locus, output unquoted "NA" in chi-square column
    if( final && k != dim - 1 ){
      *outputstream  << "NA" << sep << endl;
    }
  }//end loop over alleles
  if(final){
    double chisq=0.0;
    try{
      if(dim2==dim) chisq = GaussianQuadraticForm(ScoreVector, ObservedInfo, dim);
      else chisq = GaussianMarginalQuadraticForm( dim2, ScoreVector, ObservedInfo, dim );//marginalise over first dim2 elements
      if(chisq < 0.0)
	*outputstream << "NA" << sep << endl;
      else *outputstream << double2R(chisq) << sep << endl;
    }
    catch(...){//in case ObservedInfo is rank deficient
      *outputstream  << "NA" << sep << endl;
    }
  }
  //TODO:?? output p-value for chisq
	
  delete[] ScoreVector;
  delete[] CompleteInfo;
  delete[] ObservedInfo;
}

// the next function is for ancestry assoc tests
void ScoreTests::OutputTestsForLocusLinkage( int iterations, ofstream* outputstream, const Vector_s& PopLabels,
					     const double* Score, const double* VarScore,
					     const double* Score2, const double* Info, string separator )
//used for affectedsonly test and ancestry association test
//Score2 = Score^2
{
  //when there are two populations: 
  //we only compute scores for the 2nd population.
  
  int KK = options->getPopulations(), k1 = 0;
  
  if(KK == 2 ){
    KK = 1;k1 = 1;
  }
  
  double VU = 0.0, EU = 0.0, missing = 0.0, complete = 0.0;
  for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    for( int k = 0; k < KK; k++ ){//end at 1 for 2pops
      *outputstream << "\"" << (*Lociptr)(j)->GetLabel(0) << "\"" << separator;
      *outputstream << "\"" << PopLabels[k+k1] << "\"" << separator; //need offset to get second poplabel for 2pops
      
      EU = Score[ j*KK + k] / ( iterations );
      VU = VarScore[ j*KK + k ] / ( iterations );
      missing = Score2[ j*KK + k ] / ( iterations ) - EU * EU + VU;
      complete =  Info[ j*KK + k ] / ( iterations );

      *outputstream << double2R(EU, 3)                                << separator//score
		    << double2R(complete, 3)                          << separator//complete info
		    << double2R(complete - missing, 3)                << separator//observed info
		    << double2R(100*(complete - missing)/complete, 2) << separator;//%observed info
      if(complete > 0.0){
	*outputstream << double2R(100*(VU/complete), 2)                 << separator//%missing info attributable to locus ancestry
		      << double2R(100*(missing-VU)/complete, 2)         << separator;//%remainder of missing info      
	if(complete - missing > 0.0){
	  double zscore = EU / sqrt( complete - missing );
	  double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
	  *outputstream << double2R(zscore,3)  << separator << double2R(pvalue) << separator << endl;
	}
	else *outputstream << "NaN" << separator << "NaN" << separator << endl;
      }
      else{
	*outputstream << "NaN" << separator << "NaN" << separator << "NaN" << separator << "NaN" << separator << endl; 
      }
    }
  }
}

void ScoreTests::ROutput(){
  const int numPrintedIterations = NumOutputs;//(options->getTotalSamples() - options->getBurnIn()) / (options->getSampleEvery() * 10);
  /**
   * writes out the dimensions and labels of the 
   * R object previously written to allelicAssocScoreStream
   */
  int count;
  if(options->getTestForAllelicAssociation()){
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
    int KK = options->getPopulations();
    if(KK ==2 )KK = 1;
    vector<int> dimensions(3,0);
    dimensions[0] = 10;
    dimensions[1] = Lociptr->GetNumberOfCompositeLoci() * KK;
    dimensions[2] = (int)(numPrintedIterations);
    
    vector<string> labels(dimensions[0],"");
    labels[0] = "Locus";
    labels[1] = "Population";
    labels[2] = "Score";
    labels[3] = "CompleteInfo";
    labels[4] = "ObservedInfo";
    labels[5] = "PercentInfo";
    labels[6] = "Missing1";
    labels[7] = "Missing2";
    labels[8] = "StdNormal";
    labels[9] = "PValue";
    
    R_output3DarrayDimensions(&ancestryAssociationScoreStream,dimensions,labels);
  }
  
  /**
   * writes out the dimensions and labels of the 
   * R-matrix previously written to affectedsOnlyScoreStream
   */
  if (options->getTestForAffectedsOnly()){
//     int KK = options->getPopulations();
//     if(KK ==2 )KK = 1;
//     vector<int> dimensions(3,0);
//     dimensions[0] = 10;
//     dimensions[1] = Lociptr->GetNumberOfCompositeLoci() * KK;
//     dimensions[2] = (int)(numPrintedIterations);
    
//     vector<string> labels(dimensions[0],"");
//     labels[0] = "Locus";
//     labels[1] = "Population";
//     labels[2] = "Score";
//     labels[3] = "CompleteInfo";
//     labels[4] = "ObservedInfo";
//     labels[5] = "PercentInfo";
//     labels[6] = "Missing1";
//     labels[7] = "Missing2"; 
//     labels[8] = "StdNormal";
//     labels[9] = "PValue";
//    R_output3DarrayDimensions(&affectedsOnlyScoreStream,dimensions,labels);
    AffectedsOnlyScoreTest.ROutput(numPrintedIterations);
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
  *stream << ")," << endl;
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
  return AffectedsOnlyScoreTest;}


void ScoreTests::OutputLikelihoodRatios(const char* const filename, int iterations, const Vector_s& PopLabels){
  AffectedsOnlyScoreTest.OutputLikRatios(filename, iterations, PopLabels, *Lociptr);
}
