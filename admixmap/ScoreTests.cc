/** 
 *   ADMIXMAP
 *   Scoretests.cc 
 *   Class implements the following score tests:
 *   (1) Score test for admixture association (admixturescoretest)
 *   (2) Score test for allelic association
 *   (3) Score test for within-halpotype association
 *   (4) Score test for linkage with locus ancestry
 *   (5) Affecteds-only score test for linkage with locus ancestry
 *   (6) Score test for residual allelic association between adjacent pairs of linked loci
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
#include "gsl/gsl_cdf.h" //for residual allelic assoc test
#include "StringConvertor.h"

using namespace std;

ScoreTests::ScoreTests(){
  SumAncestryScore = 0;
  SumAncestryInfo = 0;
  SumAncestryVarScore = 0;
  SumAncestryScore2 = 0;

  SumAffectedsScore = 0;
  SumAffectedsScore2 = 0;
  SumAffectedsInfo = 0;
  SumAffectedsVarScore = 0;

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

  AdmixtureScore = 0; 
  AdmixtureInfo = 0; 
  SumAdmixtureScore = 0; 
  SumAdmixtureScore2 = 0;
  SumAdmixtureInfo = 0;

  SumAlleleScore= 0;
  SumAlleleScore2 = 0;
  SumAlleleInfo = 0;

  options = 0;
  individuals = 0;

}

ScoreTests::~ScoreTests(){
  delete[] SumAncestryScore;
  delete[] SumAncestryInfo;
  delete[] SumAncestryVarScore;
  delete[] SumAncestryScore2;

  delete[] SumAffectedsScore;
  delete[] SumAffectedsScore2;
  delete[] SumAffectedsInfo;
  delete[] SumAffectedsVarScore;

  delete[] AdmixtureScore;
  delete[] AdmixtureInfo;
  delete[] SumAdmixtureScore;
  delete[] SumAdmixtureInfo;
  delete[] SumAdmixtureScore2;

  //TODO: delete these properly
  delete[] ScoreWithinHaplotype;
  delete[] InfoWithinHaplotype;
  delete[] LocusLinkageAlleleScore;
  delete[] LocusLinkageAlleleInfo;
  delete[] SumLocusLinkageAlleleScore;
  delete[] SumLocusLinkageAlleleScore2;
  delete[] SumLocusLinkageAlleleInfo;
  delete[] locusObsIndicator;

  delete[] SumAlleleScore;
  delete[] SumAlleleScore2;
  delete[] SumAlleleInfo;

}

void ScoreTests::Initialise(AdmixOptions* op, const IndividualCollection* const indiv, const Genome* const Loci, 
			    const std::string *PLabels,
			    LogWriter &Log){
  options = op;
  individuals = indiv;
  chrm = Loci->getChromosomes();
  Lociptr = Loci;
  Log.setDisplayMode(Quiet);

  int K = options->getPopulations();
  int L = Lociptr->GetNumberOfCompositeLoci();

  /*----------------------
    | admixture association |
    -----------------------*/
  //TODO check conditions on this test
  if( options->getTestForAdmixtureAssociation() ){
    if ( strlen( options->getAssocScoreFilename() ) ){
      //can't use this yet as this test doesn't write an R object
      //OpenFile(Log, &assocscorestream, options->getAssocScoreFilename(), "Tests for admixture association");
      assocscorestream.open( options->getAssocScoreFilename(), ios::out );
      if( !assocscorestream ){
	Log.setDisplayMode(On);
	Log << "ERROR: Couldn't open admixturescorefile\n";
	exit( 1 );}
      else {
	Log << "Writing tests for admixture association to: " << options->getAssocScoreFilename() << "\n";
	assocscorestream << setiosflags( ios::fixed );

	int NumOutcomeVars = indiv->getNumberOfOutcomeVars();
	AdmixtureScore = new double[K * NumOutcomeVars];
	SumAdmixtureScore = new double[K * NumOutcomeVars];
	SumAdmixtureScore2 = new double[K * NumOutcomeVars];
	AdmixtureInfo = new double[K * NumOutcomeVars];
	SumAdmixtureInfo = new double[K * NumOutcomeVars];
	fill(AdmixtureScore, AdmixtureScore + K * NumOutcomeVars, 0.0);
	fill(AdmixtureInfo, AdmixtureInfo + K * NumOutcomeVars, 0.0);
	fill(SumAdmixtureScore, SumAdmixtureScore + K * NumOutcomeVars, 0.0);
	fill(SumAdmixtureScore2, SumAdmixtureScore2 + K * NumOutcomeVars, 0.0);
	fill(SumAdmixtureInfo, SumAdmixtureInfo + K * NumOutcomeVars, 0.0);
      }
    }
    else{
      Log.setDisplayMode(On);
      Log << "No admixturescorefile given\n";
      exit(1);}
  }

  /*------------------------------------
    |affecteds only linkage with ancestry |
    ------------------------------------*/ 
  if( options->getTestForAffectedsOnly() ){
    OpenFile(Log, &affectedsOnlyScoreStream, options->getAffectedsOnlyScoreFilename(), "Affected-only tests for association");
    //KK is the number of populations for which to perform test. For 2way admixture, we only want 2nd population.
    int KK = K;
    if(K == 2)KK = 1;
    
    SumAffectedsScore2 = new double[L * KK];
    SumAffectedsScore = new double[L * KK];
    SumAffectedsVarScore = new double[L * KK];
    SumAffectedsInfo = new double[L * KK];
    fill(SumAffectedsScore, SumAffectedsScore +L*KK, 0.0);
    fill(SumAffectedsScore2, SumAffectedsScore2 +L*KK, 0.0);
    fill(SumAffectedsInfo, SumAffectedsInfo + L*KK, 0.0);
    fill(SumAffectedsVarScore, SumAffectedsVarScore + L*KK, 0.0);
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
    OpenFile(Log, &allelicAssocScoreStream, options->getAllelicAssociationScoreFilename(), "Test for allelic association");
    
    dim_ = new unsigned[L];//dimensions of arrays
    LocusLinkageAlleleScore = new double*[L];
    LocusLinkageAlleleInfo = new double*[L];
    SumLocusLinkageAlleleScore2 = new double*[L];
    SumLocusLinkageAlleleScore = new double*[L];
    SumLocusLinkageAlleleInfo = new double*[L];
    
    ScoreWithinHaplotype = new double**[ L ];
    InfoWithinHaplotype = new double**[ L ];
    SumScoreWithinHaplotype = new double*[L];
    SumScore2WithinHaplotype = new double*[L];
    SumInfoWithinHaplotype = new double*[L];
    locusObsIndicator = new bool[L];
    
    
    //search for loci with no observed genotypes
    for(int j = 0; j < L; ++j){
      locusObsIndicator[j] = false;
      for(int i = 0; i < indiv->getSize(); ++i){
	if(!indiv->getIndividual(i)->GenotypeIsMissing(j)){
	  locusObsIndicator[j] = true;
	}
      }
    }
    
    for( int j = 0; j < L; j++ ){
      int NumberOfLoci = (*Lociptr)(j)->GetNumberOfLoci();
      
      if(NumberOfLoci > 1 )//haplotype
	dim_[j] = 1;
      else if((*Lociptr)(j)->GetNumberOfStates() == 2 )//simple diallelic locus
	dim_[j] = 1;
      else
	dim_[j] = (*Lociptr)(j)->GetNumberOfStates();//simple multiallelic locus

      //next two lines may not be necessary as these arrays are sized later
      LocusLinkageAlleleScore[j] = new double[ dim_[j] + K ];
      LocusLinkageAlleleInfo[j] = new double[( dim_[j] + K) * (dim_[j] + K )];

      SumLocusLinkageAlleleScore[j] = new double[ dim_[j] ];
      SumLocusLinkageAlleleInfo[j] = new double[( dim_[j]) * (dim_[j] )];
      SumLocusLinkageAlleleScore2[j] = new double[( dim_[j]) * (dim_[j] )];
      fill(SumLocusLinkageAlleleScore[j], SumLocusLinkageAlleleScore[j]+dim_[j], 0.0);
      fill(SumLocusLinkageAlleleInfo[j], SumLocusLinkageAlleleInfo[j] + dim_[j]*dim_[j], 0.0);
      fill(SumLocusLinkageAlleleScore2[j], SumLocusLinkageAlleleScore2[j] + dim_[j]*dim_[j], 0.0);
      
      if( NumberOfLoci > 1 ){
	ScoreWithinHaplotype[ j ] = new double*[NumberOfLoci];
	InfoWithinHaplotype[ j ] = new double*[NumberOfLoci];
	for(int jj = 0; jj < NumberOfLoci; ++jj){
	  ScoreWithinHaplotype[j][jj] = new double[1 + K];
	  InfoWithinHaplotype[j][jj] = new double[(1 + K)*(1 + K)];
	}
	SumScoreWithinHaplotype[j] = new double[ NumberOfLoci ];
	SumScore2WithinHaplotype[j] = new double[ NumberOfLoci ];
	SumInfoWithinHaplotype[j] = new double[ NumberOfLoci ];
	fill(SumScoreWithinHaplotype[j], SumScoreWithinHaplotype[j]+NumberOfLoci, 0.0);
	fill(SumScore2WithinHaplotype[j], SumScore2WithinHaplotype[j]+NumberOfLoci, 0.0);
	fill(SumInfoWithinHaplotype[j], SumInfoWithinHaplotype[j]+NumberOfLoci, 0.0);
      }
    }
  }  

  /*----------------------
    | haplotype association |
    ----------------------*/  
  if( strlen( options->getHaplotypeAssociationScoreFilename() ) ){
    if(Lociptr->GetTotalNumberOfLoci() >= Lociptr->GetNumberOfCompositeLoci()){//cannot test for SNPs in Haplotype if only simple loci
      OpenFile(Log, &HaplotypeAssocScoreStream, options->getHaplotypeAssociationScoreFilename(), "Tests for haplotype associations");
    }
    else {
      op->setTestForHaplotypeAssociation(false);
      Log << "ERROR: Cannot test for haplotype associations if all loci are simple\n" << "This option will be ignored\n";
    }
  }

  /*--------------------------------
    | Residual Allelic association  |
    -------------------------------*/
  if(options->getTestForResidualAllelicAssoc()){
    //open output file
    OpenFile(Log, &ResAlleleScoreFile, options->getResidualAllelicAssocScoreFilename(), "Tests for residual allelic association");

    SumAlleleScore = new double**[Lociptr->GetNumberOfChromosomes()];
    SumAlleleScore2 = new double**[Lociptr->GetNumberOfChromosomes()];
    SumAlleleInfo = new double**[Lociptr->GetNumberOfChromosomes()];
    for(unsigned j = 0; j < Lociptr->GetNumberOfChromosomes(); ++j){
      SumAlleleScore[j] = new double*[chrm[j]->GetSize()-1];
      SumAlleleScore2[j] = new double*[chrm[j]->GetSize()-1];
      SumAlleleInfo[j] = new double*[chrm[j]->GetSize()-1];
      for(unsigned k = 0; k < chrm[j]->GetSize()-1; ++k){
	int locus = chrm[j]->GetLocus(k);
	unsigned dim = ((*Lociptr)(locus)->GetNumberOfStates()-1) * ((*Lociptr)(locus+1)->GetNumberOfStates()-1);
	SumAlleleScore[j][k] = new double[dim];
	fill(SumAlleleScore[j][k], SumAlleleScore[j][k]+dim, 0.0);
	SumAlleleScore2[j][k] = new double[dim*dim];
	fill(SumAlleleScore2[j][k], SumAlleleScore2[j][k]+dim*dim, 0.0);
	SumAlleleInfo[j][k] = new double[dim*dim];
	fill(SumAlleleInfo[j][k], SumAlleleInfo[j][k]+dim*dim, 0.0);
      }

    }
  }
  
  InitialiseAssocScoreFile(PLabels);
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
//Initialise ergodic average score file
void ScoreTests::InitialiseAssocScoreFile(const std::string *PLabels){
  if( options->getTestForAdmixtureAssociation() ){
    PopLabels = PLabels;
    assocscorestream << "Ergodic averages of score statistic for populations:\n";
    for( int i = 0; i < options->getPopulations(); i++ ){
      assocscorestream << PopLabels[i] << " ";
      if( !i )
	assocscorestream << " ";
    }
    assocscorestream << "\ncomplete  missing   statistic  ";
    for( int i = 1; i < options->getPopulations(); i++ )
      assocscorestream << "complete  missing   statistic ";
    assocscorestream << endl;
  }
  
}

void ScoreTests::Reset(){
  //resets arrays holding sums of scores and info over individuals to zero; invoked at start of each iteration after burnin.

  if( options->getTestForAdmixtureAssociation() ){
    fill(AdmixtureScore, AdmixtureScore + options->getPopulations() * individuals->getNumberOfOutcomeVars(), 0.0);
    fill(AdmixtureInfo, AdmixtureInfo + options->getPopulations() * individuals->getNumberOfOutcomeVars(), 0.0);
  }

  if( options->getTestForAllelicAssociation() ){
    int K = options->getPopulations();
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      fill(LocusLinkageAlleleScore[j], LocusLinkageAlleleScore[j]+dim_[j]+K, 0.0);
      fill(LocusLinkageAlleleInfo[j], LocusLinkageAlleleInfo[j]+(dim_[j]+K)*(dim_[j]+K), 0.0);

      if((* Lociptr)(j)->GetNumberOfLoci() > 1 ){
	for(int jj = 0; jj < (* Lociptr)(j)->GetNumberOfLoci(); ++jj){
	  fill(ScoreWithinHaplotype[j][jj], ScoreWithinHaplotype[j][jj]+K+1, 0.0);
	  fill(InfoWithinHaplotype[j][jj], InfoWithinHaplotype[j][jj]+(K+1)*(K+1), 0.0);
	}
      }
    }
  }
}

void ScoreTests::SetAllelicAssociationTest(const std::vector<double> &alpha0){
  /*
    Invokes merging of haplotypes in CompositeLocus and resizes arrays for allelic association test accordingly  
    alpha0 = alpha[0], pop admixture dirichlet params, from Latent
  */

  //first scale alphas so they sum to 1
  double alphaScaled[options->getPopulations()];
  double sum  = accumulate(alpha0.begin(), alpha0.end(), 0.0, std::plus<double>());//sum of alpha0 over pops
  for( int k = 0; k < options->getPopulations(); k++ )
    alphaScaled[k] = alpha0[k] / sum;

  //merge rare haplotypes
  for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){//skip simple loci
      (*Lociptr)(j)->SetDefaultMergeHaplotypes( alphaScaled);
      
      //resize arrays for Allelic association test
      delete[] LocusLinkageAlleleScore[j];
      delete[] LocusLinkageAlleleInfo[j];
      delete[] SumLocusLinkageAlleleScore[j];
      delete[] SumLocusLinkageAlleleInfo[j];
      delete[] SumLocusLinkageAlleleScore2[j];

      // dim is set to 1 if single diallelic locus, to number of
      // merged haplotypes if compound locus, to number of states
      // if >2 alleles
      dim_[j] = (*Lociptr)(j)->GetNumberOfMergedHaplotypes();
      LocusLinkageAlleleScore[j] = new double[ dim_[j] + options->getPopulations()];
      LocusLinkageAlleleInfo[j] = new double[( dim_[j] + options->getPopulations()) * (dim_[j] + options->getPopulations() )];

      SumLocusLinkageAlleleScore[j] = new double[ dim_[j] ];
      SumLocusLinkageAlleleInfo[j] = new double[( dim_[j]) * (dim_[j] )];
      SumLocusLinkageAlleleScore2[j] = new double[( dim_[j]) * (dim_[j] )];
      fill(SumLocusLinkageAlleleScore[j], SumLocusLinkageAlleleScore[j]+dim_[j], 0.0);
      fill(SumLocusLinkageAlleleInfo[j], SumLocusLinkageAlleleInfo[j] + dim_[j]*dim_[j], 0.0);
      fill(SumLocusLinkageAlleleScore2[j], SumLocusLinkageAlleleScore2[j] + dim_[j]*dim_[j], 0.0);
    }
  }
}

// ****************************** UPDATES ****************************

void ScoreTests::Update(double dispersion)
//Note: dispersion = dispersion in regression model
//                 = precision for linear reg, 1.0 for logistic
{
  Reset();//set sums over individuals to zero
  
  double DInvLink;
  
  //----------------------------------
  // Update Scores for each individual
  //----------------------------------
  for( int i = 0; i < individuals->getSize(); i++ ){
    if(options->getNumberOfOutcomes() > 0){//if regressionmodel
      Individual* ind = individuals->getIndividual(i);
      double YMinusEY = individuals->getOutcome(0, i) - individuals->getExpectedY(i);//individual outcome - its expectation
      DInvLink = individuals->DerivativeInverseLinkFunction(i);
      
      //admixture association
      if( options->getTestForAdmixtureAssociation() && (options->getNumberOfOutcomes() == 1) ){
	UpdateScoreForAdmixtureAssociation(ind->getAdmixtureProps(), YMinusEY,dispersion, DInvLink);
      }
      //allelic association
      if( options->getTestForAllelicAssociation() )
	UpdateScoreForAllelicAssociation( ind, YMinusEY,dispersion, DInvLink);
    }
  }
  
  //-----------------------------
  //Accumulate Scores, Info etc. over iterations
  //-----------------------------
  
  /*----------------------
    | admixture association |
    -----------------------*/
  if( options->getTestForAdmixtureAssociation() ){
    int K = options->getPopulations();
    int NumOutcomeVars = individuals->getNumberOfOutcomeVars();
    //SumAdmixtureScore += AdmixtureScore;
    transform(AdmixtureScore, AdmixtureScore + K*NumOutcomeVars, SumAdmixtureScore, SumAdmixtureScore, std::plus<double>());
    //SumAdmixtureInfo += AdmixtureInfo;
    transform(AdmixtureInfo, AdmixtureInfo + K*NumOutcomeVars, SumAdmixtureInfo, SumAdmixtureScore, std::plus<double>());
    
    for( int k = 0; k < K; k++ )
      for( int kk = 0; kk < NumOutcomeVars; kk++ )
	SumAdmixtureScore2[ k*NumOutcomeVars + kk ] += AdmixtureScore[ k*NumOutcomeVars + kk ] * AdmixtureScore[ k*NumOutcomeVars + kk ];
  }
  
  for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    
    /*-------------------------------------
      | Allelic and haplotype association  |
      -------------------------------------*/
    try{
      if( options->getTestForAllelicAssociation() ){
	//loop over simple loci within haplotypes
	if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
	  for( int l = 0; l < (*Lociptr)(j)->GetNumberOfLoci(); l++ ){
	    CentreAndSum(1, ScoreWithinHaplotype[j][l], InfoWithinHaplotype[j][l], &(SumScoreWithinHaplotype[ j ][ l ]),
			 &(SumScore2WithinHaplotype[j][l]), &(SumInfoWithinHaplotype[ j ][ l ]));
	  }
	}
      }
      
      if( (options->getTestForAllelicAssociation()  && (*Lociptr)(j)->GetNumberOfLoci() == 1) || options->getTestForHaplotypeAssociation() ){
	if(locusObsIndicator[j]){//skip loci with no observed genotypes
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
      Individual::SumScoresForAncestry(j, SumAncestryScore, SumAncestryInfo, SumAncestryScore2, SumAncestryVarScore);
    } 
    /*------------------------------------
      |affecteds-only linkage with ancestry |
      ------------------------------------*/ 
    if( options->getTestForAffectedsOnly() ){
      Individual::SumScoresForLinkageAffectedsOnly(j, SumAffectedsScore, SumAffectedsVarScore, SumAffectedsScore2, SumAffectedsInfo);
    }
  }
}

void ScoreTests::UpdateScoreForAdmixtureAssociation( const double* const Theta, double YMinusEY, double phi, double DInvLink)
{
  //Updates score and info for score test for admixture association
  double x;
  int NumOutcomeVars = individuals->getNumberOfOutcomeVars();
  for( int k = 0; k < options->getPopulations(); k++ ){
    if( options->isRandomMatingModel() )
      x = 0.5 * ( Theta[k] + Theta[ options->getPopulations() + k ] );
    else
      x = Theta[k];
    AdmixtureScore[ k*NumOutcomeVars] += phi * x * YMinusEY;// only for 1st outcomevar 
    AdmixtureInfo[ k*NumOutcomeVars] += phi * x * x *DInvLink;
  }
}

void ScoreTests::UpdateScoreForAllelicAssociation( const Individual* const ind, double YMinusEY, double phi, double DInvLink)
{
  int locus = 0;
  //int K = options->getPopulations();

  for(unsigned int j = 0; j < Lociptr->GetNumberOfChromosomes(); j++ ){
    for(unsigned int jj = 0; jj < chrm[j]->GetSize(); jj++ ){
      if(!ind->GenotypeIsMissing(locus)){//skip loci with missing genotypes as hap pairs have not been sampled for these
	//retrieve sampled hap pair from Individual
	const int* happair = ind->getSampledHapPair(locus);
	const unsigned numStates = (*Lociptr)(locus)->GetNumberOfStates();
	const unsigned numLoci = (*Lociptr)(locus)->GetNumberOfLoci();
	
	// count alleles / haplotypes      
	vector<int> counts;
	
	// if diallelic, evaluate score for allele 2 only, otherwise evaluate score for all alleles or haplotypes
	
	// special case for SNP (counts has size 1)
	if( numStates == 2 ){
	  //sets counts to -1, 0 or 1 according to whether happair is 11, 12 or 22
	  counts.push_back((*Lociptr)(locus)->getAlleleCounts(2, happair)[0] - 1 );
	  
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
  int K = options->getPopulations();
  double* x = new double[ K + 1 ];


  x[ K ] = 1.0;
  for( int k = 0; k < K - 1; k++ )
    x[ k + 1 ] = ind->getAdmixtureProps()[k];

  for( int l = 0; l < (*Lociptr)(j)->GetNumberOfLoci(); l++ ){
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
  int K = options->getPopulations();
  unsigned dim = Counts.size();
  double* x = new double[ K + dim ];

  // ** Set x co-ordinate for regression parameter under test
  copy(Counts.begin(), Counts.end(), x);
 
  // ** Set x-co-ordinates of covariates in model
  x[ dim ] = 1.0;//intercept
  //for( int k = 1; k < K; k++ )
  //x[ dim+k ] = admixtureProps[k];
  copy(admixtureProps+1, admixtureProps+K, x+dim+1);


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
  double *cscore = new double[dim];
  double *cinfo = new double[dim*dim];

  CentredGaussianConditional( dim, score, info, cscore, cinfo, options->getPopulations()+dim );
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

void ScoreTests::UpdateScoresForResidualAllelicAssociation(const double* const* AlleleFreqs){
  for(unsigned c = 0; c < Lociptr->GetNumberOfChromosomes(); ++c)
    for(unsigned j = 0; j < chrm[c]->GetSize()-1; ++j){
      int abslocus = chrm[c]->GetLocus(j);
      UpdateScoresForResidualAllelicAssociation(c, j, AlleleFreqs[abslocus], AlleleFreqs[abslocus+1]);
    }
}


int delta(int i, int j){
  return (i == j);
}

void ScoreTests::UpdateScoresForResidualAllelicAssociation(int c, int locus,  
							   const double* const AlleleFreqsA, const double* const AlleleFreqsB){
  int abslocus = chrm[c]->GetLocus(locus);//number of this locus
  int M = (*Lociptr)(abslocus)->GetNumberOfStates()-1;
  int N = (*Lociptr)(abslocus+1)->GetNumberOfStates()-1;
  int dim = M*N;
  if(dim == 1)UpdateScoresForResidualAllelicAssociation_1D(c, locus, AlleleFreqsA, AlleleFreqsB);
  else{
    
    int Populations = options->getPopulations();
    //compute frequencies of last haplotype as AlleleFreqs omits it
    double lastFreqsA[Populations];fill(lastFreqsA, lastFreqsA+Populations, 1.0);
    double lastFreqsB[Populations];fill(lastFreqsB, lastFreqsB+Populations, 1.0);
    for(int k = 0; k < Populations; ++k){
      for(int i = 0; i < M; ++i)lastFreqsA[k] -= AlleleFreqsA[i*Populations+k];
      for(int j = 0; j < N; ++j)lastFreqsB[k] -= AlleleFreqsB[j*Populations+k];
    }

    vector<double> AlleleScore(dim);//temporary array to hold sums of scores over individuals. 
    vector<double> AlleleInfo(dim*dim);
    //Note that the dimensions vary between calls
    fill(AlleleScore.begin(), AlleleScore.end(), 0.0);
    fill(AlleleInfo.begin(), AlleleInfo.end(), 0.0);
    int count = 0;
    
    int ancA[2];//ancestry at A
    int ancB[2];//ancestry at B

    for(int i = 0; i < individuals->getSize(); ++i){
      Individual* ind = individuals->getIndividual(i);
      if(!ind->GenotypeIsMissing(abslocus)){//skip loci with missing genotypes as hap pairs have not been sampled for these
	ind->GetLocusAncestry(c, locus, ancA);
	ind->GetLocusAncestry(c, locus+1, ancB);
	const int* hA = ind->getSampledHapPair(abslocus);//realized hap pair at locus A
	const int* hB = ind->getSampledHapPair(abslocus+1);//realized hap pair at locus B
	for(int g = 0; g < 2; ++g){
	  //if(ancA[g] == ancB[g]){ 
	  ++count;//count number of gametes with ancestry states the same at both loci
	  
	  for(int i = 0; i < M; ++i)
	    for(int j = 0; j < N; ++j){//i and j index rows
	      //update score
	      AlleleScore[i*N +j] += ( delta(hA[g], i) - delta(hA[g], M) ) * ( delta(hB[g], j) - delta(hB[g], N) ) //observed
		- ( AlleleFreqsA[i*Populations + ancA[g]] - lastFreqsA[ancA[g]] ) 
		* ( AlleleFreqsB[j*Populations + ancB[g]] - lastFreqsB[ancB[g]]);// - expected
	      
	      for(int m = 0; m < M ; ++m)for(int n = 0; n < N; ++n){//m and n index columns
		//update info
		AlleleInfo[(i*N+j) * dim + (m*N+n)] += 
		  ( delta(i,m)*AlleleFreqsA[i*Populations + ancA[g]] + lastFreqsA[ancA[g]] ) 
		  * ( delta(j,n)*AlleleFreqsB[j*Populations + ancB[g]] + lastFreqsB[ancB[g]] );
	      }
	    }
	  
	}//end gamete loop
	//}//end condition (ancestry equal)
      }
    }//end ind loop
    
    //accumulate score and score squared
    for(int j = 0; j < dim; ++j){
      SumAlleleScore[c][locus][j] += AlleleScore[j];
      for(int jj = 0; jj < dim; ++jj){
	SumAlleleScore2[c][locus][j*dim +jj] += AlleleScore[j]*AlleleScore[jj];
	SumAlleleInfo[c][locus][j*dim +jj] = AlleleInfo[j*dim + jj] * (double)count;
      }
    }
  }//end else
}
void ScoreTests::UpdateScoresForResidualAllelicAssociation_1D(int c, int locus,  
							   const double* const AlleleFreqsA, const double* const AlleleFreqsB){
  int ancA[2];
  int ancB[2];
  //int Populations = options->getPopulations();
  int count = 0;
  double AlleleScore = 0.0;

  int abslocus = chrm[c]->GetLocus(locus);
  for(int i = 0; i < individuals->getSize(); ++i){
    Individual* ind = individuals->getIndividual(i);
    if(!ind->GenotypeIsMissing(abslocus)){//skip loci with missing genotypes as hap pairs have not been sampled for these
      ind->GetLocusAncestry(c, locus, ancA);
      ind->GetLocusAncestry(c, locus+1, ancB);
      const int* hA = ind->getSampledHapPair(abslocus);//realized hap pair at locus A
      const int* hB = ind->getSampledHapPair(abslocus+1);//realized hap pair at locus B
      
      for(int g = 0; g < 2; ++g){
	//if(ancA[g] == ancB[g]){ 
	  ++count;//count number of gametes with ancestry states the same at both loci
	  int h = (hA[g] == hB[g]);//indicator for homozygosity
	  double phiA = AlleleFreqsA[ancA[g]];//frequency of first allele in this gamete's ancestry state at locus A
	  double phiB = AlleleFreqsB[ancB[g]];
	  
	  //h - (1-h) evaluates to 1 if allele states equal, -1 if not
	  //expected counts factorises to  --\/
	  AlleleScore += h - (1-h) - (2*phiA - 1.0)*(2*phiB - 1.0);
	  
	  //}//end condition on equal ancestry states
      }//end gamete loop
    }
  }//end individual loop
  //accumulate score and score squared

  SumAlleleScore[c][locus][0] += AlleleScore;
  SumAlleleScore2[c][locus][0] += AlleleScore*AlleleScore;
  SumAlleleInfo[c][locus][0] += (double)count;

}
// ********** OUTPUT **********************************************************

void ScoreTests::Output(int iteration, const std::string * PLabels){
  PopLabels = PLabels;
  int iterations = iteration - options->getBurnIn();
  //Allelic association
  if( options->getTestForAllelicAssociation() )    {
    for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
       //case of simple locus
      if((* Lociptr)(j)->GetNumberOfLoci() == 1 )
	OutputTestsForAllelicAssociation(iterations, &allelicAssocScoreStream, j, dim_[j], 
					 SumLocusLinkageAlleleScore[j], SumLocusLinkageAlleleScore2[j], 
					 SumLocusLinkageAlleleInfo[j], ",");
      //case of haplotype
      else
	OutputTestsForAllelicAssociation(iterations, &allelicAssocScoreStream, j, (*Lociptr)(j)->GetNumberOfLoci(), 
					 SumScoreWithinHaplotype[ j ], SumScore2WithinHaplotype[ j ], 
					 SumInfoWithinHaplotype[ j ], ",");
    }//end j loop over comp loci
  }
  
  //haplotype association
  if( options->getTestForHaplotypeAssociation() ){
    OutputTestsForHaplotypeAssociation( iterations, &HaplotypeAssocScoreStream, "," );
  }
  
  //ancestry association
  if( options->getTestForLinkageWithAncestry() ){
    
    OutputTestsForLocusLinkage( iterations, &ancestryAssociationScoreStream,
				SumAncestryScore, SumAncestryVarScore,
				SumAncestryScore2, SumAncestryInfo, "," );
  }
  //affectedonly
  if( options->getTestForAffectedsOnly() ){
    OutputTestsForLocusLinkage( iterations, &affectedsOnlyScoreStream,
				SumAffectedsScore, SumAffectedsVarScore,
				SumAffectedsScore2, SumAffectedsInfo, "," );
  }
  
  //admixture association
  if( options->getTestForAdmixtureAssociation() ){
    OutputAdmixtureScoreTest( iterations );
  }
  //residual allelic association
  if(options->getTestForResidualAllelicAssoc()){

    OutputTestsForResidualAllelicAssociation(iterations, &ResAlleleScoreFile, ",");
  }
}

void ScoreTests::WriteFinalTables(){
  ofstream finalTable;
  int iterations = options->getTotalSamples() - options->getBurnIn();
  //ancestry association
  if( options->getTestForLinkageWithAncestry() ){
    string filename(options->getResultsDir());
    filename.append("/TestsAncestryAssocFinal.txt");
    finalTable.open(filename.c_str(), ios::out);
    finalTable <<"Locus\tPopulation\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tMissing1\tMissing2\tStdNormal\tPValue\n";
    OutputTestsForLocusLinkage( iterations, &finalTable,
				SumAncestryScore, SumAncestryVarScore,
				SumAncestryScore2, SumAncestryInfo, "\t" );
    finalTable.close();
  }
  //affectedonly
  if( options->getTestForAffectedsOnly() ){
    string filename(options->getResultsDir());
    filename.append("/TestsAffectedsOnlyFinal.txt");
    finalTable.open(filename.c_str(), ios::out);
    finalTable <<"Locus\tPopulation\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tMissing1\tMissing2\tStdNormal\tPValue\n";
    OutputTestsForLocusLinkage( iterations, &finalTable,
				SumAffectedsScore, SumAffectedsVarScore,
				SumAffectedsScore2, SumAffectedsInfo, "\t" );
    finalTable.close();
  }

  //residual allelic association
  if(options->getTestForResidualAllelicAssoc()){
    string filename(options->getResultsDir());
    filename.append("/ResidualLDTestFinal.txt");
    finalTable.open(filename.c_str(), ios::out);
    finalTable << "Loci\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tdf\tChiSquared\tPValue\n";
    OutputTestsForResidualAllelicAssociation(iterations, &finalTable, "\t");
    finalTable.close();
  }
  //haplotype association
  if( options->getTestForHaplotypeAssociation() ){
    string filename(options->getResultsDir());
    filename.append("/HaplotypeAssocTestsFinal.txt");
    finalTable.open(filename.c_str(), ios::out);
    finalTable << "Locus\tHaplotype\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tStdNormal\tPValue\tChiSquare\n";
    OutputTestsForHaplotypeAssociation( iterations, &finalTable, "\t" );
    finalTable.close();
  }
  //Allelic association
  if( options->getTestForAllelicAssociation() )    {
    string filename(options->getResultsDir());
    filename.append("/AllelicAssocTestsFinal.txt");
    finalTable.open(filename.c_str(), ios::out);
    finalTable << "Locus\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tStdNormal\tPValue\n";
    for(unsigned j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
       //case of simple locus
      if((* Lociptr)(j)->GetNumberOfLoci() == 1 )
	OutputTestsForAllelicAssociation(iterations, &finalTable, j, dim_[j], 
					 SumLocusLinkageAlleleScore[j], SumLocusLinkageAlleleScore2[j], 
					 SumLocusLinkageAlleleInfo[j], "\t");
      //case of haplotype
      else
	OutputTestsForAllelicAssociation(iterations, &finalTable, j, (*Lociptr)(j)->GetNumberOfLoci(), 
					 SumScoreWithinHaplotype[ j ], SumScore2WithinHaplotype[ j ], 
					 SumInfoWithinHaplotype[ j ], "\t");
    }//end j loop over comp loci
    finalTable.close();
  }
}

// next few function calculate score tests from the cumulative sums of
// the score, score squared, and information score and info can be
// scalars, or respectively a vector and a matrix should have just one
// method or class to do this
void ScoreTests::OutputAdmixtureScoreTest(int iterations)
{
  int NumOutcomeVars = individuals->getNumberOfOutcomeVars();
  for( int j = 0; j < options->getPopulations(); j++ ){
    for( int jj = 0; jj < NumOutcomeVars; jj++ ){
      double EU = SumAdmixtureScore[ j*NumOutcomeVars + jj ] / ( iterations );
      double complete = SumAdmixtureInfo[ j*NumOutcomeVars + jj ] / ( iterations );
      double missing = SumAdmixtureScore2[ j*NumOutcomeVars + jj ] / ( iterations ) - EU * EU;
      assocscorestream.width(9);
      assocscorestream << setprecision(6) << double2R(complete) << " ";
      assocscorestream.width(9);
      assocscorestream << setprecision(6) << double2R(missing) << " ";
      assocscorestream.width(9);
      assocscorestream << setprecision(6) << double2R(EU / sqrt( complete - missing )) << " ";
    }
  }
  assocscorestream << endl;
}

void ScoreTests::OutputTestsForHaplotypeAssociation( int iterations, ofstream* outputstream, string sep )
// loops over composite loci that have 2 or more simple loci, and calculates score tests for each 
// haplotype separately, together with a summary chi-square
{
  int NumberOfMergedHaplotypes;
  const int *hap;
  double *ScoreVector, *CompleteMatrix, *ObservedMatrix;
  
  for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){ // compound locus contains 2 or more simple loci
      
      ScoreVector = new double[dim_[j]];
      copy(SumLocusLinkageAlleleScore[j], SumLocusLinkageAlleleScore[j]+dim_[j], ScoreVector);
      scale_matrix(ScoreVector, 1.0/( iterations), dim_[j], 1);
      
      CompleteMatrix = new double[dim_[j]*dim_[j]];
      copy(SumLocusLinkageAlleleInfo[j], SumLocusLinkageAlleleInfo[j]+ dim_[j]*dim_[j], CompleteMatrix);
      scale_matrix(CompleteMatrix, 1.0/( iterations), dim_[j], dim_[j]);
      
      ObservedMatrix = new double[dim_[j]*dim_[j]];
      for(unsigned d1 = 0; d1 < dim_[j]; ++d1)for(unsigned d2 = 0; d2 < dim_[j]; ++d2)
	ObservedMatrix[d1*dim_[j] + d2] = CompleteMatrix[d1*dim_[j]+d2] + ScoreVector[d1]*ScoreVector[d2] -
	  SumLocusLinkageAlleleScore2[j][d1*dim_[j]+d2]/( iterations );
      
      NumberOfMergedHaplotypes = dim_[j];
      for( int k = 0; k < NumberOfMergedHaplotypes; k++ ){
	*outputstream  << "\"" << (*Lociptr)(j)->GetLabel(0) << "\"" << sep;
	if( k < NumberOfMergedHaplotypes - 1 ){
	  hap = (*Lociptr)(j)->GetHapLabels(k);
	  *outputstream  << "\"";
	  for( int kk = 0; kk < (*Lociptr)(j)->GetNumberOfLoci() - 1; kk++ ){
	    *outputstream  << hap[kk] << "-";
	  }
	  *outputstream  << hap[(*Lociptr)(j)->GetNumberOfLoci() - 1] << "\"" << sep;
	}
	else
	  *outputstream  << "\"others\"" << sep;
	*outputstream  << double2R(ScoreVector[k], 3) << sep
		       << double2R(CompleteMatrix[k*dim_[j]+k], 3) << sep
		       << double2R(ObservedMatrix[k*dim_[j]+k], 3) << sep
		       << double2R(100*ObservedMatrix[k*dim_[j]+k] / CompleteMatrix[k*dim_[j]+k], 2) << sep;//%Observed Info
	double zscore = ScoreVector[ k ] / sqrt( ObservedMatrix[k*dim_[j]+k] );
	*outputstream  << double2R(zscore, 2) << sep;//z-score
	double pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
	*outputstream << double2R(pvalue, 2) << sep;
	// if not last allele at locus, output unquoted "NA" in chi-square column
	if( k != NumberOfMergedHaplotypes - 1 ){
	  *outputstream  << "NA" << sep << endl;
	}
      }//end loop over haplotypes
      // calculate summary chi-square statistic
      double chisq = 0.0;
      try{
	chisq = GaussianConditionalQuadraticForm( NumberOfMergedHaplotypes - 1, ScoreVector, ObservedMatrix, dim_[j] );
	*outputstream  << double2R(chisq, 2) << sep << endl;
      }
      catch(...){
	*outputstream  << "NaN" << sep << endl;// if ObservedMatrix is rank deficient
      }

      delete[] ScoreVector;
      delete[] CompleteMatrix;
      delete[] ObservedMatrix;
    }//end if
  }//end loop over loci
}

void ScoreTests::OutputTestsForAllelicAssociation( int iterations, ofstream* outputstream, int locus, unsigned dim, 
						   const double* score, const double* scoresq, const double* info, string sep)
{
  double Score, CompleteInfo, MissingInfo, ObservedInfo, PercentInfo, zscore, pvalue;
  for(unsigned a = 0; a < dim; ++a){
    Score = score[a] / ( iterations );
    CompleteInfo = info[a] / ( iterations );
    MissingInfo = scoresq[a] / ( iterations ) - Score * Score;
    ObservedInfo = CompleteInfo - MissingInfo;
    
    string locuslabel = (*Lociptr)(locus)->GetLabel(a);
    if(dim==1 || (*Lociptr)(locus)->GetNumberOfLoci()>1) *outputstream << "\"" << locuslabel << "\"" << sep;
    else *outputstream << "\"" << locuslabel<< "("<<a+1<<")\""<< sep;
    *outputstream << double2R(Score, 3)        << sep
		  << double2R(CompleteInfo, 3) << sep
		  << double2R(ObservedInfo, 3) << sep;
    if(CompleteInfo > 0.0) {
      PercentInfo = 100*ObservedInfo / CompleteInfo;
      zscore = Score / sqrt( ObservedInfo );
      pvalue = 2.0 * gsl_cdf_ugaussian_P(-fabs(zscore));
      *outputstream << double2R(PercentInfo, 2) << sep
		    << double2R(zscore, 2)    << sep 
		    << double2R(pvalue, 2) << sep << endl;
    }
    else{
      *outputstream << "NaN" << sep << "NaN" << sep << endl;
    }
    
  }
}

void ScoreTests::OutputTestsForLocusLinkage( int iterations, ofstream* outputstream,
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
  
  double VU, EU, missing, complete;
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
	  *outputstream << double2R(zscore, 2)  << separator << double2R(pvalue, 2) << separator << endl;
	}
	else *outputstream << "NaN" << separator << "NaN" << separator << endl;
	}
      else{
	*outputstream << "NaN" << separator << "NaN" << separator << "NaN" << separator << "NaN" << separator << endl; 
      }
    }
  }
}

void ScoreTests::OutputTestsForResidualAllelicAssociation(int iterations, ofstream* outputstream, string separator){
  double *Score, *ObservedInfo;
  *outputstream << setiosflags(ios::fixed) << setprecision(3) ;
  for(unsigned int c = 0; c < Lociptr->GetNumberOfChromosomes(); c++ )
    for(unsigned j = 0; j < chrm[c]->GetSize()-1; ++j){
      int abslocus = chrm[c]->GetLocus(j);
      int M = (*Lociptr)(abslocus)->GetNumberOfStates()-1;
      int N = (*Lociptr)(abslocus+1)->GetNumberOfStates()-1;
      int dim = M*N;
      Score = new double[dim];
      ObservedInfo = new double[dim*dim];
      double obsinfo = 0.0;
      double compinfo = 0.0;

      for(int i = 0; i < dim; ++i){
	Score[i] = SumAlleleScore[c][j][i]/(double) iterations; //score
	for(int ii = 0; ii < dim; ++ii){
	  ObservedInfo[i*dim +ii] = SumAlleleInfo[c][j][i*dim+ii]/ (double) iterations//complete info
	    + Score[i]*Score[ii] - SumAlleleScore2[c][j][i*dim+ii]/(double)iterations;//-missing info
	}
	obsinfo += ObservedInfo[i*dim +i];//trace of Observed Info
	compinfo += SumAlleleInfo[c][j][i*dim+i]/ (double) iterations;//trace of Complete Info
      }

      //output labels
      *outputstream << "\"" << (*Lociptr)(abslocus)->GetLabel(0) << "/" << (*Lociptr)(abslocus+1)->GetLabel(0) << "\""<< separator
			 << Score[0] << separator << compinfo << separator <<obsinfo << separator;

      //compute chi-squared statistic
      double VinvU[dim];
      try{
	HH_solve(dim, ObservedInfo, Score, VinvU);
	double chisq = 0.0;
	for(int i = 0; i < dim; ++i){
	  chisq += Score[i] * VinvU[i];
	}
	
	*outputstream << (obsinfo*100)/compinfo << separator<< dim << separator;

	if(chisq < 0.0){
	  *outputstream << "NA" << separator << "NA" << separator << endl;
	}
	else {
	  //compute p-value
	  double pvalue = gsl_cdf_chisq_Q (chisq, dim);
	  *outputstream << chisq << separator << pvalue << separator << endl;
	}
      }
      catch(...){//in case ObservedInfo is rank deficient
	*outputstream  << "NA" << separator << dim 
			    << separator << "NA" << separator << "NA" << separator << endl;
      }
      delete[] Score;
      delete[] ObservedInfo;
    }//end loop over loci on chromosome
}

void ScoreTests::ROutput(){
  int numPrintedIterations = (options->getTotalSamples() - options->getBurnIn()) / (options->getSampleEvery() * 10);
  /**
   * writes out the dimensions and labels of the 
   * R object previously written to allelicAssocScoreStream
   */
  int count;
  if(options->getTestForAllelicAssociation()){
    count = 0;
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      if((* Lociptr)(j)->GetNumberOfLoci() == 1 ) count += dim_[j];
      else count += (*Lociptr)(j)->GetNumberOfLoci();
    }         
    vector<int> dimensions(3,0);
    dimensions[0] = 7;
    dimensions[1] = count;
    dimensions[2] = (int)(numPrintedIterations);
    vector<string> labels(dimensions[0],"");
    labels[0] = "Locus";
    labels[1] = "Score";
    labels[2] = "CompleteInfo";
    labels[3] = "ObservedInfo";
    labels[4] = "PercentInfo";
    labels[5] = "StdNormal";
    labels[6] = "Pvalue";
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
    dimensions[0] = 9;
    dimensions[1] = count;
    dimensions[2] = (int)(numPrintedIterations);
    vector<string> labels(dimensions[0],"");
    labels[0] = "Locus";
    labels[1] = "Haplotype";
    labels[2] = "Score";
    labels[3] = "CompleteInfo";
    labels[4] = "ObservedInfo";
    labels[5] = "PercentInfo";
    labels[6] = "StdNormal";
    labels[7] = "PValue";
    labels[8] = "ChiSquare";
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
    R_output3DarrayDimensions(&affectedsOnlyScoreStream,dimensions,labels);
  }

  /**
   * writes out the dimensions and labels of the 
   * R-matrix previously written to ResAlleleScoreFile
   */
  if(options->getTestForResidualAllelicAssoc()){
    vector<int> dimensions(3,0);
    dimensions[0] = 8;
    dimensions[1] = Lociptr->GetNumberOfCompositeLoci() - Lociptr->GetNumberOfChromosomes();
    dimensions[2] = (int)(numPrintedIterations);
    
    vector<string> labels(dimensions[0],"");
    labels[0] = "Loci";
    labels[1] = "Score";
    labels[2] = "CompleteInfo";
    labels[3] = "ObservedInfo";
    labels[4] = "PercentInfo";
    labels[5] = "df";
    labels[6] = "ChiSquare";
    labels[7] = "P-value";

    R_output3DarrayDimensions(&ResAlleleScoreFile, dimensions, labels);
  }

}

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
string ScoreTests::double2R( double x )
{
  if( isnan(x) )
    return "NaN";
  else{
    stringstream ret;
    ret << x;
    return( ret.str() );
  }
}
string ScoreTests::double2R( double x, int precision )
{
  if( isnan(x) )
    return "NaN";
  else{
    stringstream ret;

    ret << setiosflags(ios::fixed) << setprecision(precision)<< x;
    return( ret.str() );
  }
}
