/** 
 *   ADMIXMAP
 *   Scoretests.cc 
 *   Class implements the following score tests:
 *   (1) Score test for admixture association (admixturescoretest)
 *   (2) Score test for allelic association
 *   (3) Score test for within-halpotype association
 *   (4) Score test for linkage with locus ancestry
 *   (5) Affecteds-only score test for linkage with locus ancestry
 *   Copyright (c) 2005 LSHTM
 *  
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include "ScoreTests.h"
#include <numeric>


using namespace std;

ScoreTests::ScoreTests(){
  Matrix_d null_Matrix_d(1,1);

  //   LocusLinkageScore = 0;
  //   LocusLinkageInfo = 0;
  //   SumLocusLinkageScore = 0;
  //   SumLocusLinkageScore2 = 0;
  //   SumLocusLinkageInfo = 0;

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

  ancestryAssociationScoreStream = 0;
  SNPsAssociationScoreStream = 0;
  affectedsOnlyScoreStream = 0;

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
}

void ScoreTests::Initialise(AdmixOptions * op, IndividualCollection *indiv, Genome *Loci, Chromosome **c,std::string *PLabels,
			    LogWriter *log){
  options = op;
  individuals = indiv;
  chrm = c;
  Lociptr=Loci;
  Logptr = log;

  int K = options->getPopulations();
  int L = Lociptr->GetNumberOfCompositeLoci();

  /*----------------------
    | admixture association |
    -----------------------*/
  //TODO check conditions on this test
  if( options->getTestForAdmixtureAssociation() ){
    if ( strlen( options->getAssocScoreFilename() ) ){
      assocscorestream.open( options->getAssocScoreFilename(), ios::out );
      if( !assocscorestream ){
	Logptr->logmsg(true,"ERROR: Couldn't open admixturescorefile\n");
	exit( 1 );}
      else {
	Logptr->logmsg(true,"Writing tests for admixture association to: ");    
	Logptr->logmsg(true,options->getAssocScoreFilename());
	Logptr->logmsg(true,"\n");
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
      Logptr->logmsg(true,"No admixturescorefile given\n");
      exit(1);}
  }

  /*------------------------------------
    |affecteds only linkage with ancestry |
    ------------------------------------*/ 
  if( options->getTestForAffectedsOnly() ){
    affectedsOnlyScoreStream = new ofstream(options->getAffectedsOnlyScoreFilename());
    if( !affectedsOnlyScoreStream ){
      Logptr->logmsg(true,"ERROR: Couldn't open affectedsonlyscorefile");
      Logptr->logmsg(true,options->getAffectedsOnlyScoreFilename());
      Logptr->logmsg(true,"\n");
      exit( 1 );//remove?
    }
    else{
      Logptr->logmsg(true,"Writing affected only tests for association to ");
      Logptr->logmsg(true,options->getAffectedsOnlyScoreFilename());
      Logptr->logmsg(true,"\n");
      *affectedsOnlyScoreStream << "structure(.Data=c(" << endl;
	
      //KK is the number of populations for which to perform test. For 2way admixture, we only want 2nd population.
      int KK = K;
      if(K == 2)KK = 1;
      //Individual::InitialiseAffectedsOnlyScores(L, KK);
      SumAffectedsScore2 = new double[L * KK];
      SumAffectedsScore = new double[L * KK];
      SumAffectedsVarScore = new double[L * KK];
      SumAffectedsInfo = new double[L * KK];
      fill(SumAffectedsScore, SumAffectedsScore +L*KK, 0.0);
      fill(SumAffectedsScore2, SumAffectedsScore2 +L*KK, 0.0);
      fill(SumAffectedsInfo, SumAffectedsInfo + L*KK, 0.0);
      fill(SumAffectedsVarScore, SumAffectedsVarScore + L*KK, 0.0);
    }
  }

  /*-----------------------
    | Linkage with ancestry  |
    -----------------------*/
  if( options->getTestForLinkageWithAncestry() ){
    ancestryAssociationScoreStream = new ofstream(options->getAncestryAssociationScoreFilename());
    if( !ancestryAssociationScoreStream ){
      Logptr->logmsg(true,"ERROR: Couldn't open ancestry association scorefile\n");
      exit( 1 );
    }
    else{
      Logptr->logmsg(true,"Writing tests for locus linkage to ");
      Logptr->logmsg(true,options->getAncestryAssociationScoreFilename());
      Logptr->logmsg(true,"\n");
      *ancestryAssociationScoreStream << "structure(.Data=c(" << endl;
    }
	
    //code for old method
    //     LocusLinkageScore = new Matrix_d[L];
    //     LocusLinkageInfo = new Matrix_d[L];
    //     SumLocusLinkageScore = new Matrix_d[L];
    //     SumLocusLinkageScore2 = new Matrix_d[L];
    //     SumLocusLinkageInfo = new Matrix_d[L];
    //     for( int j = 0; j < L; j++ ){
    //       LocusLinkageScore[j].SetNumberOfElements(2 * K, 1 );
    //       LocusLinkageInfo[j].SetNumberOfElements( 2 * K, 2 * K );
    //       SumLocusLinkageScore2[j].SetNumberOfElements(K, 1 );
    //       SumLocusLinkageScore[j].SetNumberOfElements(K, 1 );
    //       SumLocusLinkageInfo[j].SetNumberOfElements(K, 1 );
    //     }

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
    genescorestream.open( options->getAllelicAssociationScoreFilename(), ios::out );
    if( !genescorestream ){
      Logptr->logmsg(true,"ERROR: Couldn't open locusscorefile\n");
      exit( 1 );
    }
    else{
      Logptr->logmsg(false,"Test for allelic association written to ");
      Logptr->logmsg(false,options->getAllelicAssociationScoreFilename());
      Logptr->logmsg(false,"\n");
      genescorestream << "structure(.Data=c(" << endl;
    }
    
    int kk;
    LocusLinkageAlleleScore = new Matrix_d[L];
    LocusLinkageAlleleInfo = new Matrix_d[L];
    SumLocusLinkageAlleleScore2 = new Matrix_d[L];
    SumLocusLinkageAlleleScore = new Matrix_d[L];
    SumLocusLinkageAlleleInfo = new Matrix_d[L];
    
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
	if(!indiv->getIndividual(i)->IsMissing(j)){
	  locusObsIndicator[j] = true;
	}
      }
    }
    
    for( int j = 0; j < L; j++ ){
      int NumberOfLoci = (*Lociptr)(j)->GetNumberOfLoci();
      
      if(NumberOfLoci > 1 )
	kk = 1;
      else if((*Lociptr)(j)->GetNumberOfStates() == 2 )
	kk = 1;
      else
	kk = (*Lociptr)(j)->GetNumberOfStates();

      LocusLinkageAlleleScore[j].SetNumberOfElements( kk + K, 1 );
      LocusLinkageAlleleInfo[j].SetNumberOfElements( kk + K, kk + K );
      SumLocusLinkageAlleleScore2[j].SetNumberOfElements( kk, kk );
      SumLocusLinkageAlleleScore[j].SetNumberOfElements( kk, 1 );
      SumLocusLinkageAlleleInfo[j].SetNumberOfElements( kk, kk );
      
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

  /*---------------------------------
    | Misspecified allele frequencies  |
    ---------------------------------*/
  //moved to new class

  /*-------------------
    | SNPs in haplotype  |
    -------------------*/  
  if( strlen( options->getTestsForSNPsInHaplotypeOutputFilename() ) ){
    if( !options->getTestForAllelicAssociation() ){
      Logptr->logmsg(true,"Can't Test For SNPs In Haplotype if allelicassociationscorefile is not specified\n");
      exit(1);
    }
    if(Lociptr->GetTotalNumberOfLoci() > Lociptr->GetNumberOfCompositeLoci()){//cannot test for SNPs in Haplotype if only simple loci
      SNPsAssociationScoreStream = new ofstream( options->getTestsForSNPsInHaplotypeOutputFilename(), ios::out );
      if( !SNPsAssociationScoreStream ){
	Logptr->logmsg(true,"ERROR: Couldn't open haplotypeassociationscorefile\n");
	exit( 1 );
      }
      else{
	Logptr->logmsg(true,"Tests for associations of SNPs in haplotypes written to\n");
	Logptr->logmsg(true,options->getTestsForSNPsInHaplotypeOutputFilename());
	Logptr->logmsg(true,"\n");
	*SNPsAssociationScoreStream << "structure(.Data=c(" << endl;
      }
    }
    else {
      options->setTestForSNPsInHaplotype(false);
      Logptr->logmsg(true, "ERROR: Cannot test for SNPs in Haplotype if all loci are simple\n");
      Logptr->logmsg(true, "This option will be ignored\n");
    }

  }


  if( options->getTextIndicator() )InitialiseAssocScoreFile(PLabels);
}

//Initialise ergodic average score file
void ScoreTests::InitialiseAssocScoreFile(std::string *PLabels){
  if( options->getTestForAdmixtureAssociation() ){
    PopLabels = PLabels;
    assocscorestream << "Ergodic averages of score statistic for populations:\n";
    for( int i = 0; i < options->getPopulations(); i++ ){
      if(options->IsPedFile())
	assocscorestream << "\"" << PopLabels[i] << "\"" << " ";
      else
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
  if( options->getTestForAdmixtureAssociation() ){
    fill(AdmixtureScore, AdmixtureScore + options->getPopulations() * individuals->getNumberOfOutcomeVars(), 0.0);
    fill(AdmixtureInfo, AdmixtureInfo + options->getPopulations() * individuals->getNumberOfOutcomeVars(), 0.0);
  }

  if( options->getTestForAllelicAssociation() ){
    int K = options->getPopulations();
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      LocusLinkageAlleleScore[j].SetElements(0);
      LocusLinkageAlleleInfo[j].SetElements(0);

      if((* Lociptr)(j)->GetNumberOfLoci() > 1 ){
	for(int jj = 0; jj < (* Lociptr)(j)->GetNumberOfLoci(); ++jj){
	  fill(ScoreWithinHaplotype[j][jj], ScoreWithinHaplotype[j][jj]+K+1, 0.0);
	  fill(InfoWithinHaplotype[j][jj], InfoWithinHaplotype[j][jj]+(K+1)*(K+1), 0.0);
	}
      }
    }
  }

  //if( options->getTestForLinkageWithAncestry() ){
  //     LocusLinkageScore.SetElements(0);
  //     LocusLinkageInfo.SetElements(0);
  //}

}

void ScoreTests::SetAllelicAssociationTest(std::vector<double> &alpha0){
  //alpha0 = alpha[0], pop admixture dirichlet params, from Latent
  double pp[options->getPopulations()];
  double sum  = accumulate(alpha0.begin(), alpha0.end(), 0.0, std::plus<double>());//sum of alpha0 over pops
  for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
      for( int k = 0; k < options->getPopulations(); k++ )
	pp[k] = alpha0[k] / sum;
      (*Lociptr)(j)->SetDefaultMergeHaplotypes( pp);
      
      //possible memory leaks here
      int kk = (*Lociptr)(j)->GetNumberOfMergedHaplotypes();
      LocusLinkageAlleleScore[j].SetNumberOfElements( kk + options->getPopulations(),1 );
      LocusLinkageAlleleInfo[j].SetNumberOfElements( kk + options->getPopulations(), kk + options->getPopulations() );
      SumLocusLinkageAlleleScore2[j].SetNumberOfElements( kk, kk );
      SumLocusLinkageAlleleScore[j].SetNumberOfElements( kk, 1 );
      SumLocusLinkageAlleleInfo[j].SetNumberOfElements( kk, kk );
    }
  }
}

void ScoreTests::Update(double dispersion)
//Note: dispersion = dispersion in regression model
//                 = precision for linear reg, 1.0 for logistic
{
  Reset();

  double DInvLink;

  //----------------------------------
  // Update Scores for each individual
  //----------------------------------
  for( int i = 0; i < individuals->getSize(); i++ ){
    Individual* ind = individuals->getIndividual(i);
    double YMinusEY = individuals->getOutcome(0, i) - individuals->getExpectedY(i);
    DInvLink = individuals->DerivativeInverseLinkFunction(options->getAnalysisTypeIndicator(),i);

    //admixture association
    if( options->getTestForAdmixtureAssociation() && (options->getAnalysisTypeIndicator()>1 && options->getAnalysisTypeIndicator()<5) ){
      UpdateScoreForAdmixtureAssociation(ind->getAdmixtureProps(), YMinusEY,dispersion, DInvLink);
    }
    //allelic association
    if( options->getTestForAllelicAssociation() )
      UpdateScoreForAllelicAssociation( ind, YMinusEY,dispersion, DInvLink);
  }

  //the rest formerly was the function UpdateScoreStats
  //-----------------------------
  //Accumulate Scores, Info etc. over iterations
  //-----------------------------
  int kk;
  Matrix_d score, info; //used to transform score and info for regression-based tests

  /*---------------------------------
    | Misspecified allele frequencies  |
    ---------------------------------*/
  //moved to new class

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

  if( options->getAnalysisTypeIndicator() != 1 ){


    if( options->getTestForAllelicAssociation() ) SumScoreForWithinHaplotypeAssociation();


    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){

      /*----------------------
	| Allelic association  |
	-----------------------*/
      if( (options->getTestForAllelicAssociation()  && (*Lociptr)(j)->GetNumberOfLoci() == 1) || options->getTestForSNPsInHaplotype() ){
	if(locusObsIndicator[j]){//skip loci with no observed genotypes
	  // kk is set to 1 if single diallelic locus, to number of
	  // merged haplotypes if compound locus, to number of states
	  // if >2 alleles
	  if( (*Lociptr)(j)->GetNumberOfStates() == 2 )
	    kk = 1;
	  else if( (*Lociptr)(j)->GetNumberOfLoci() > 1 )
	    kk = (*Lociptr)(j)->GetNumberOfMergedHaplotypes();
	  else
	    kk = (*Lociptr)(j)->GetNumberOfStates();
	  // correct for posterior covariance between regression parameters
	  CentredGaussianConditional( kk, LocusLinkageAlleleScore[j], LocusLinkageAlleleInfo[j], &score, &info );
	  SumLocusLinkageAlleleScore[j] += score;
	  SumLocusLinkageAlleleScore2[j] += score * score.Transpose();
	  SumLocusLinkageAlleleInfo[j] += info;
	}
      }

      /*-----------------------
	| Linkage with ancestry  |
	-----------------------*/
      if( options->getTestForLinkageWithAncestry() ){
	//code for old method
	// 	CentredGaussianConditional( options->getPopulations(),
	// 				  LocusLinkageScore[j], LocusLinkageInfo[j],
	// 				  &score, &info );

	// 	SumLocusLinkageScore[j] += score;
	// 	for( int kk = 0; kk < options->getPopulations(); kk++ )
	// 	  SumLocusLinkageScore2[j]( kk, 0 ) += score( kk, 0 ) * score( kk, 0 );
	// 	SumLocusLinkageInfo[j] += info.GetDiagonal().ColumnMatrix();

	Individual::SumScoresForAncestry(j,  
					 SumAncestryScore, SumAncestryInfo, SumAncestryScore2, SumAncestryVarScore);

      } 
      /*------------------------------------
	|affecteds-only linkage with ancestry |
	------------------------------------*/ 
      if( options->getTestForAffectedsOnly() ){
	Individual::SumScoresForLinkageAffectedsOnly(j, SumAffectedsScore, SumAffectedsVarScore, SumAffectedsScore2, SumAffectedsInfo);
      }
    }
  }
}

void ScoreTests::UpdateScoreForAllelicAssociation( Individual* ind,double YMinusEY, double phi, double DInvLink)
//Note: EY0 = ExpectedY(0)(i,0), lambda0 = lambda(0)
{

  int dim = 0; //? this removes warning, but isn't necessarily the Right Thing
  int locus = 0;
  Vector_i ancestry;

  //double IndividualOutcomeMinusP = individuals->getOutcome(0)( i, 0 ) - EY0;

  for(unsigned int j = 0; j < Lociptr->GetNumberOfChromosomes(); j++ ){
    for(unsigned int jj = 0; jj < chrm[j]->GetSize(); jj++ ){

      if( (*Lociptr)(locus)->GetNumberOfStates() == 2 ){
	dim = options->getPopulations() + 1;
      } else if( (*Lociptr)(locus)->GetNumberOfLoci() == 1 ){
	dim = (*Lociptr)(locus)->GetNumberOfStates() + options->getPopulations();
      } else {
	dim = (*Lociptr)(locus)->GetNumberOfMergedHaplotypes() + options->getPopulations();
      }


      // Set x-co-ordinates of covariates in model
      Matrix_d cov_x_coord(dim,1);
      cov_x_coord( dim - 1, 0 ) = 1;
      for( int k = 0; k < options->getPopulations() - 1; k++ ){
	cov_x_coord( k + dim - options->getPopulations(), 0 ) = ind->getAdmixtureProps()[k];
      }
     
      // Set x co-ordinate for regression parameter under test

      // special case for SNP
      if( (*Lociptr)(locus)->GetNumberOfStates() == 2 ){
	//? next line relies on encoding system
	cov_x_coord( 0, 0 ) = (int)(ind->getGenotype(locus)[0][0]) + (int)(ind->getGenotype(locus)[0][1] - 3.0);
	// general case for simple locus
      } else if((*Lociptr)(locus)->GetNumberOfLoci() == 1 ){
	for( int k = 0; k < (*Lociptr)(locus)->GetNumberOfStates(); k++ ){
	  if( (int)(ind->getGenotype(locus)[0][0]) == k + 1 )
	    cov_x_coord( k, 0 )++;
	  if( (int)(ind->getGenotype(locus)[0][1]) == k + 1 )
	    cov_x_coord( k, 0 )++;
	}
	// general case for composite locus
      } else {
	UpdateScoreForWithinHaplotypeAssociation(ind, locus, YMinusEY,phi , DInvLink);
	if( options->getTestForSNPsInHaplotype() ){
	  ancestry = ind->GetLocusAncestry( j, jj );
	  int hap[2] = {0,0};
	  (*Lociptr)(locus)->SampleHapPair(hap, ind->getPossibleHapPairs(locus), ancestry);
	  for( int k = 0; k <(*Lociptr)(locus)->GetNumberOfStates(); k++ ){
	    if( hap[0] == k )
	      cov_x_coord( (*Lociptr)(locus)->GetMergedHaplotype(k), 0 )++;
	    if( hap[1] == k )
	      cov_x_coord((*Lociptr)(locus)->GetMergedHaplotype(k), 0 )++;
	  }
	}
      }

      // Check if missing genotype
      if( !(ind->IsMissing(locus)) ){
	LocusLinkageAlleleScore[ locus ] += cov_x_coord * YMinusEY * phi;
	LocusLinkageAlleleInfo[ locus ] += (cov_x_coord * cov_x_coord.Transpose()) * phi * DInvLink;
      }
      locus++;
    }
  }
}

// need to change usage of chrm to make this work
// void ScoreTests::UpdateScoreForAncestryOld( Individual* ind, double Y,  int regressonindicator, double EY, double lambda0)
// //This function is obsolete, replaced by UpdateScoreForAncestry, and is preserved here in case it is needed for reference.
// //Note: EY = ExpectedY(0)(i,0), lambda0 = lambda(0)
// {
//   // regressonindicator = 0 => linear regression model
//   // regressonindicator = 1 => logistic regression model

//   int locus = 0,dim = 2 * options->getPopulations();
//   Vector_i ancestry;
//   Matrix_d X;

//   X.SetNumberOfElements(dim, 1);

//   for( int j = 0; j < chrm->size(); j++ ){
//     for( int jj = 0; jj < (*chrm)(j)->GetSize(); jj++ ){
//       X.SetElements(0); 
        
//       // Set x-co-ordinates of covariates in model
//       X( dim - 1, 0 ) = 1;
//       for( int k = 0; k < options->getPopulations() - 1; k++ ){
// 	X( k + options->getPopulations(), 0 ) = ind->getAdmixtureprops()( k, 0 );
// 
//       }
     
//       // Set x co-ordinate for regression parameter under test
// 	ancestry = ind->GetLocusAncestry( j, jj );
// 	X( ancestry(0), 0 )++;
// 	X( ancestry(1), 0 )++;
 
// 	//should be possible to make this independent of regression type as for ancestry association score test
// 	if( regressonindicator ){//logistic
// 	  LocusLinkageScore( locus ) += X * (Y - EY);
// 	  LocusLinkageInfo( locus ) += (X * X.Transpose()) * (EY * ( 1 - EY ));
// 	} else {//linear
// 	  LocusLinkageScore( locus ) += lambda0 * (Y - EY) * X;
// 	  LocusLinkageInfo( locus ) += (X * X.Transpose()) * lambda0;
// 	}
//       locus++;
//     }
//   }
// }


void ScoreTests::UpdateScoreForAdmixtureAssociation( double  *Theta, double YMinusEY, double phi, double DInvLink)
{
  //Updates score and info for score test for admixture association
  double x;
  int NumOutcomeVars = individuals->getNumberOfOutcomeVars();
  for( int k = 0; k < options->getPopulations(); k++ ){
    if( options->isRandomMatingModel() )
      x = 0.5 * ( Theta[k] + Theta[ options->getPopulations() + k ] );
    else
      x = Theta[k];
    AdmixtureScore[ k*NumOutcomeVars] += phi * x * YMinusEY;//only for 1st outcomevar ??
    AdmixtureInfo[ k*NumOutcomeVars] += phi * x * x *DInvLink;
  }
}

// This method calculates score for allelic association at each simple locus within a compound locus
void ScoreTests::UpdateScoreForWithinHaplotypeAssociation( Individual *ind, int j, double YMinusEY,double phi, double DInvLink)
{
  int K = options->getPopulations();
  double x[ K + 1 ];
 
  double* info = new double[ (K + 1) * (K + 1 )];

  x[ K ] = 1.0;
  for( int k = 0; k < K - 1; k++ )
    x[ k + 1 ] = ind->getAdmixtureProps()[k];

  for( int l = 0; l < (*Lociptr)(j)->GetNumberOfLoci(); l++ ){
    x[0] = GetAllele2CountsInHaplotype(l, ind->getGenotype(j));
    for( int k = 0; k < K + 1; k++ ){
      //if( x(0) != 99 )ScoreWithinHaplotype[j](l)( k, 0 ) += phi * x(k) * ( individuals->getOutcome(0)( i, 0 ) - p );
      if( x[0] != 99 )ScoreWithinHaplotype[j][l][ k ] += phi * x[k] * YMinusEY;
      for( int kk = 0; kk < K + 1; kk++ )
	info[ k*(K+1) + kk ] = x[ k ] * x[ kk ];
    }
    if( x[0] != 99 ){
      scale_matrix(info, phi*DInvLink, K+1, K+1);//info *= phi*DInvLink
      add_matrix(InfoWithinHaplotype[j][l], info, K+1, K+1);//INfoWithinHaplotype += info
    }
  }
  delete[] info;
}

/**
 * Called only by UpdateScoresForSNPsWithinHaplotype in ScoreTests
 * Given an unordered genotype, returns number of copies of allele
 * 2 a simple locus in a composite locus.
 * Used to test individual loci in haplotype for association.
 * 
 * genotype - a 2way array in which each element of
 * alleles coded as unsigned integers numbered.  
 * 
 *returns:
 * the number of copies of allele 2 at each locus.
 *
 * n.b. this function is only useful in composite loci composed of diallelic simple loci
 * should be generalized to deal with multi-allelic loci
 */

int ScoreTests::GetAllele2CountsInHaplotype(int locus, unsigned short **genotype)
{
  /**
   * AlleleCounts is the number of 2 alleles at a
   * locus in haplotype.  Used to test individual loci in haplotype
   * for association.  Only use for haplotypes made up of SNPs.
   */

  int AlleleCounts = 0;

  if(genotype[locus][0]!=0 && genotype[locus][1]!=0){
    if(genotype[locus][0] == 2){
      AlleleCounts++;
    }
    if(genotype[locus][1] == 2){
      AlleleCounts++;
    }
  } else {
    AlleleCounts = 99;
  }
  return AlleleCounts;
}


void ScoreTests::SumScoreForWithinHaplotypeAssociation()
{
  int kk = 1;
  //Matrix_d score1, score2, info, Vbb, Vab, Vaa;
  double *score = new double[1];
  double *info = new double[1];

  for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
      for( int l = 0; l < (*Lociptr)(j)->GetNumberOfLoci(); l++ ){
	//score1 = ScoreWithinHaplotype[j][l].SubMatrix( 0, kk - 1, 0, 0 );
	//score2 = ScoreWithinHaplotype[j][l].SubMatrix( kk, kk + options->getPopulations() - 1, 0, 0 );
	//Vaa = InfoWithinHaplotype[j][l].SubMatrix( 0, kk - 1, 0, kk - 1 );
	//Vbb = InfoWithinHaplotype[j][l].SubMatrix( kk, kk + options->getPopulations() - 1,
	//					   kk, kk + options->getPopulations() - 1 );
	//Vab = InfoWithinHaplotype[j][l].SubMatrix( 0, kk - 1, kk, kk + options->getPopulations() - 1 );
	//Vbb.InvertUsingLUDecomposition();
	//score1 = score1 - Vab * Vbb * score2;
	//info = Vaa - Vab * Vbb * Vab.Transpose();
	CentredGaussianConditional( kk, ScoreWithinHaplotype[j][l], InfoWithinHaplotype[j][l],
				    score, info, options->getPopulations()+1 );
	SumScoreWithinHaplotype[ j ][ l ] += score[0];
	SumScore2WithinHaplotype[ j ][ l ] += score[0] * score[0];
	SumInfoWithinHaplotype[ j ][ l ] += info[0];
      }
    }
  }
  delete[] score;
  delete[] info;
}

void ScoreTests::Output(int iteration,std::string * PLabels){
  PopLabels = PLabels;
  //Allelic association
  if( options->getTestForAllelicAssociation() ){
    OutputTestsForAllelicAssociation( iteration );
  }
  //haplotype association
  if( options->getTestForSNPsInHaplotype() ){
    OutputTestsForSNPsInHaplotype( iteration );
  }

  //ancestry association
  if( options->getTestForLinkageWithAncestry() ){
    //old version
    //    OutputTestsForLocusLinkage( iteration, ancestryAssociationScoreStream,
    // 				SumLocusLinkageScore,
    // 				SumLocusLinkageScore2, SumLocusLinkageInfo );
    
    OutputTestsForLocusLinkage( iteration, ancestryAssociationScoreStream,
				SumAncestryScore, SumAncestryVarScore,
				SumAncestryScore2, SumAncestryInfo );
  }
  //affectedonly
  if( options->getTestForAffectedsOnly() ){
    OutputTestsForLocusLinkage( iteration, affectedsOnlyScoreStream,
				SumAffectedsScore, SumAffectedsVarScore,
				SumAffectedsScore2, SumAffectedsInfo );
    
  }
  //misspecified allele freqs
  //moved to new class


  //admixture association
  if( options->getTestForAdmixtureAssociation() ){
    OutputAdmixtureScoreTest( iteration );
  }

}

// next few methods calculate score tests from the cumulative sums of
// the score, score squared, and information score and info can be
// scalars, or respectively a vector and a matrix should have just one
// method or class to do this
void ScoreTests::OutputTestsForSNPsInHaplotype( int iteration )
{
  int NumberOfMergedHaplotypes;
  const int *hap;
  Matrix_d ScoreMatrix, CompleteMatrix, ObservedMatrix;
  for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
      ScoreMatrix = SumLocusLinkageAlleleScore[j] / ( iteration - options->getBurnIn() );
      CompleteMatrix =  SumLocusLinkageAlleleInfo[j] / ( iteration - options->getBurnIn() );
      ObservedMatrix = CompleteMatrix + ScoreMatrix * ScoreMatrix.Transpose() - SumLocusLinkageAlleleScore2[j] / ( iteration - options->getBurnIn() );
      NumberOfMergedHaplotypes = ScoreMatrix.GetNumberOfRows();
      for( int k = 0; k < NumberOfMergedHaplotypes; k++ ){
	if(options->IsPedFile())
	  *SNPsAssociationScoreStream  << "\"" << (*Lociptr)(j)->GetLabel(0) << "\"" << ",";
	else
	  *SNPsAssociationScoreStream  << (*Lociptr)(j)->GetLabel(0) << ",";
	if( k < NumberOfMergedHaplotypes - 1 ){
	  hap = (*Lociptr)(j)->GetHapLabels(k);
	  *SNPsAssociationScoreStream  << "\"";
	  for( int kk = 0; kk < (*Lociptr)(j)->GetNumberOfLoci() - 1; kk++ ){
	    *SNPsAssociationScoreStream  << hap[kk] << "-";
	  }
	  *SNPsAssociationScoreStream  << hap[(*Lociptr)(j)->GetNumberOfLoci() - 1] << "\",";
	}
	else
	  *SNPsAssociationScoreStream  << "\"others\",";
	*SNPsAssociationScoreStream  << double2R(ScoreMatrix(k,0)) << ",";
	*SNPsAssociationScoreStream  << double2R(CompleteMatrix(k,k)) << ",";
	*SNPsAssociationScoreStream  << double2R(ObservedMatrix(k,k)) << ",";
	*SNPsAssociationScoreStream  << double2R(100*ObservedMatrix(k,k) / CompleteMatrix(k,k)) << ",";
	*SNPsAssociationScoreStream  << double2R(ScoreMatrix( k, 0 ) / sqrt( ObservedMatrix( k, k ) ));
	*SNPsAssociationScoreStream  << ",";
	// if not last allele at locus, output NA in chi-square column
	if( k != NumberOfMergedHaplotypes - 1 ){
	  *SNPsAssociationScoreStream  << "\"NA\"," << endl;
	}
      }
      // calculate summary chi-square statistic
      ObservedMatrix = ObservedMatrix.SubMatrix( 0, NumberOfMergedHaplotypes - 2,
						 0, NumberOfMergedHaplotypes - 2 );
      ScoreMatrix = ScoreMatrix.SubMatrix( 0, NumberOfMergedHaplotypes - 2, 0, 0 );
      ObservedMatrix.InvertUsingLUDecomposition();
      *SNPsAssociationScoreStream  << double2R((ScoreMatrix.Transpose() * ObservedMatrix * ScoreMatrix)(0,0)) << "," << endl;
    }
  }
}

void ScoreTests::OutputAdmixtureScoreTest(int iteration)
{
  int NumOutcomeVars = individuals->getNumberOfOutcomeVars();
  for( int j = 0; j < options->getPopulations(); j++ ){
    for( int jj = 0; jj < NumOutcomeVars; jj++ ){
      double EU = SumAdmixtureScore[ j*NumOutcomeVars + jj ] / ( iteration - options->getBurnIn() );
      double complete = SumAdmixtureInfo[ j*NumOutcomeVars + jj ] / ( iteration - options->getBurnIn() );
      double missing = SumAdmixtureScore2[ j*NumOutcomeVars + jj ] / ( iteration - options->getBurnIn() ) - EU * EU;
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

void ScoreTests::OutputTestsForAllelicAssociation( int iteration )
{
  double Score, CompleteInfo, MissingInfo, ObservedInfo, PercentInfo, zscore;
  for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
    for( int l = 0; l < (*Lociptr)(j)->GetNumberOfLoci(); l++ ){
      if((* Lociptr)(j)->GetNumberOfLoci() == 1 ){
	Score = SumLocusLinkageAlleleScore[j](0,0) / ( iteration - options->getBurnIn() );
	CompleteInfo =  SumLocusLinkageAlleleInfo[j](0,0) / ( iteration - options->getBurnIn() );
	MissingInfo = SumLocusLinkageAlleleScore2[j](0,0) / ( iteration - options->getBurnIn() ) - Score * Score;
	ObservedInfo = CompleteInfo - MissingInfo;
      }
      else{
	Score = SumScoreWithinHaplotype[ j ][ l ] / ( iteration - options->getBurnIn() );
	CompleteInfo = SumInfoWithinHaplotype[ j ][ l ] / ( iteration - options->getBurnIn() );
	MissingInfo = SumScore2WithinHaplotype[ j ][ l ] / ( iteration - options->getBurnIn() ) - Score * Score;
	ObservedInfo = CompleteInfo - MissingInfo;
      }
      if(CompleteInfo > 0.0) {
	PercentInfo = 100*ObservedInfo / CompleteInfo;
	zscore = Score / sqrt( ObservedInfo );
      }
      else{
	PercentInfo = 0.0;

      }
      if(options->IsPedFile())
	genescorestream << "\"" << (*Lociptr)(j)->GetLabel(l) << "\"" << ",";
      else
	genescorestream << (*Lociptr)(j)->GetLabel(l) << ",";
      genescorestream << double2R(Score)        << ","
		      << double2R(CompleteInfo) << ","
		      << double2R(ObservedInfo) << ",";
      if(CompleteInfo > 0.0) {
	genescorestream << double2R(100*ObservedInfo / CompleteInfo) << ",";
	genescorestream << double2R(Score / sqrt( ObservedInfo ))    << "," << endl;
      }
      else{
	genescorestream << "NaN" << ","
			<< "NaN" << ",";
      }
    }
  }
}

//this function is obsolete, was used for old ancestry score test
// void ScoreTests::OutputTestsForLocusLinkage( int iteration, ofstream* outputstream,
// 					Matrix_d *Score, Matrix_d *Score2,
// 					Matrix_d *Info )
// {
//   double EU, missing, complete;
//   for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
//     for( int k = 0; k < options->getPopulations(); k++ ){
//       if(options->IsPedFile())
// 	*outputstream << "\"" << (*Lociptr)(j)->GetLabel(0) << "\"" << ",";
//       else
// 	*outputstream << (*Lociptr)(j)->GetLabel(0) << ",";
//       *outputstream << PopLabels[k] << ",";
      
//       EU = Score[ j ]( k, 0 ) / ( iteration - options->getBurnIn() );
//       missing = Score2[ j ]( k, 0 ) / ( iteration - options->getBurnIn() ) - EU * EU;
//       complete =  Info[ j ]( k, 0 ) / ( iteration - options->getBurnIn() );
      
//       *outputstream << double2R(EU)                                << ",";
//       *outputstream << double2R(complete)                          << ",";
//       *outputstream << double2R(complete - missing)                << ",";
//       *outputstream << double2R(100*(complete - missing)/complete) << ",";
//       *outputstream << double2R(EU / sqrt( complete - missing ))   << "," << endl;
//     }
//   }
// }

void ScoreTests::OutputTestsForLocusLinkage( int iteration, ofstream* outputstream,
					     double* Score, double* VarScore,
					     double* Score2, double* Info )
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
      if(options->IsPedFile())
	*outputstream << "\"" << (*Lociptr)(j)->GetLabel(0) << "\"" << ",";
      else
	*outputstream << (*Lociptr)(j)->GetLabel(0) << ",";
      *outputstream << PopLabels[k+k1] << ","; //need offset to get second poplabel for 2pops
      
      EU = Score[ j*KK + k] / ( iteration - options->getBurnIn() );
      VU = VarScore[ j*KK + k ] / ( iteration - options->getBurnIn() );
      missing = Score2[ j*KK + k ] / ( iteration - options->getBurnIn() ) - EU * EU + VU;
      complete =  Info[ j*KK + k ] / ( iteration - options->getBurnIn() );
      
      *outputstream << double2R(EU)                                << ","//score
		    << double2R(complete)                          << ","//complete info
		    << double2R(complete - missing)                << ","//observed info
		    << double2R(100*(complete - missing)/complete) << ",";//%observed info
      if(complete > 0.0){
	*outputstream << double2R(100*(VU/complete))                 << ","//%missing info attributable to locus ancestry
		      << double2R(100*(missing-VU)/complete)         << ",";//%remainder of missing info      
	if(complete - missing > 0.0)
	  *outputstream << double2R(EU / sqrt( complete - missing ))   << "," << endl;
	else *outputstream << "NaN,"  << endl;
      }
      else{
	*outputstream << "NaN,NaN,NaN," << endl; 
      }
    }
  }
}
void ScoreTests::ROutput(){
  int numPrintedIterations = options->getTotalSamples()/ options->getSampleEvery() / 10  -  options->getBurnIn() / options->getSampleEvery() / 10;

  /**
   * writes out the dimensions and labels of the 
   * R-matrix previously written to genescorestream
   */
  int count;
  if(options->getTestForAllelicAssociation()){
    count = 0;
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      count += (*Lociptr)(j)->GetNumberOfLoci();
    }         
    vector<int> dimensions(3,0);
    dimensions[0] = 6;
    dimensions[1] = count;
    dimensions[2] = (int)(numPrintedIterations);
    vector<string> labels(dimensions[0],"");
    labels[0] = "Locus";
    labels[1] = "Score";
    labels[2] = "CompleteInfo";
    labels[3] = "ObservedInfo";
    labels[4] = "PercentInfo";
    labels[5] = "StdNormal";
    R_output3DarrayDimensions(&genescorestream,dimensions,labels);
  }

  /** 
   * writes out the dimensions and labels of the        
   * R-matrix for score tests of SNPs in haplotypes.        
   */ 
  if( options->getTestForSNPsInHaplotype()  ){      
    count = 0;
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
	count += (*Lociptr)(j)->GetNumberOfMergedHaplotypes();
      }
    }         
    vector<int> dimensions(3,0);
    dimensions[0] = 8;
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
    labels[7] = "ChiSquare";
    R_output3DarrayDimensions(SNPsAssociationScoreStream,dimensions,labels);
  }
   


  //old ancestry score test code
  //   /**
  //    * writes out the dimensions and labels of the 
  //    * R-matrix previously written to ancestryAssociationScoreStream
  //    */
  //   if (options->getTestForLinkageWithAncestry()){
  //     vector<int> dimensions(3,0);
  //     dimensions[0] = 7;
  //     dimensions[1] = Lociptr->GetNumberOfCompositeLoci() * options->getPopulations();
  //     dimensions[2] = (int)(options->getTotalSamples()/ options->getSampleEvery() / 10  -  options->getBurnIn() / options->getSampleEvery() / 10);

  //     vector<string> labels(7,"");
  //     labels[0] = "Locus";
  //     labels[1] = "Population";
  //     labels[2] = "Score";
  //     labels[3] = "CompleteInfo";
  //     labels[4] = "ObservedInfo";
  //     labels[5] = "PercentInfo";
  //     labels[6] = "StdNormal";

  //     R_output3DarrayDimensions(ancestryAssociationScoreStream,dimensions,labels);
  //   }

  /**
   * writes out the dimensions and labels of the 
   * R-matrix previously written to ancestryAssociationScoreStream
   */
  if (options->getTestForLinkageWithAncestry()){
    int KK = options->getPopulations();
    if(KK ==2 )KK = 1;
    vector<int> dimensions(3,0);
    dimensions[0] = 9;
    dimensions[1] = Lociptr->GetNumberOfCompositeLoci() * KK;
    dimensions[2] = (int)(numPrintedIterations);
     
    vector<string> labels(9,"");
    labels[0] = "Locus";
    labels[1] = "Population";
    labels[2] = "Score";
    labels[3] = "CompleteInfo";
    labels[4] = "ObservedInfo";
    labels[5] = "PercentInfo";
    labels[6] = "Missing1";
    labels[7] = "Missing2";
    labels[8] = "StdNormal";

    R_output3DarrayDimensions(ancestryAssociationScoreStream,dimensions,labels);
  }

  /**
   * writes out the dimensions and labels of the 
   * R-matrix previously written to affectedsOnlyScoreStream
   */
  if (options->getTestForAffectedsOnly()){
    int KK = options->getPopulations();
    if(KK ==2 )KK = 1;
    vector<int> dimensions(3,0);
    dimensions[0] = 9;
    dimensions[1] = Lociptr->GetNumberOfCompositeLoci() * KK;
    dimensions[2] = (int)(numPrintedIterations);
     
    vector<string> labels(9,"");
    labels[0] = "Locus";
    labels[1] = "Population";
    labels[2] = "Score";
    labels[3] = "CompleteInfo";
    labels[4] = "ObservedInfo";
    labels[5] = "PercentInfo";
    labels[6] = "Missing1";
    labels[7] = "Missing2"; 
    labels[8] = "StdNormal";

    R_output3DarrayDimensions(affectedsOnlyScoreStream,dimensions,labels);
  }
}

void
ScoreTests::R_output3DarrayDimensions(ofstream* stream,vector<int> dim,vector<string> labels)
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
string
ScoreTests::double2R( double x )
{
  if( isnan(x) )
    return "NaN";
  else{
    stringstream ret;
    ret << x;
    return( ret.str() );
  }
}
