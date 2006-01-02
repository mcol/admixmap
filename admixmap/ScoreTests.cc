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

void ScoreTests::Initialise(AdmixOptions* op, const IndividualCollection* const indiv,const Genome* const Loci, 
			    const Chromosome* const* c, const std::string *PLabels,
			    LogWriter &Log){
  options = op;
  individuals = indiv;
  chrm = c;
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
    affectedsOnlyScoreStream = new ofstream(options->getAffectedsOnlyScoreFilename());
    if( !affectedsOnlyScoreStream ){
      Log << "ERROR: Couldn't open affectedsonlyscorefile" << options->getAffectedsOnlyScoreFilename() << "\n";
      exit( 1 );//remove?
    }
    else{
      Log << "Writing affected only tests for association to " << options->getAffectedsOnlyScoreFilename() << "\n";
      *affectedsOnlyScoreStream << "structure(.Data=c(" << endl;
	
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
  }

  /*-----------------------
    | Linkage with ancestry  |
    -----------------------*/
  if( options->getTestForLinkageWithAncestry() ){
    ancestryAssociationScoreStream = new ofstream(options->getAncestryAssociationScoreFilename());
    if( !ancestryAssociationScoreStream ){
      Log << "ERROR: Couldn't open ancestry association scorefile\n";
      exit( 1 );
    }
    else{
      Log << "Writing tests for locus linkage to " << options->getAncestryAssociationScoreFilename() << "\n";
      *ancestryAssociationScoreStream << "structure(.Data=c(" << endl;
    }
	
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
      Log << "ERROR: Couldn't open locusscorefile\n";
      exit( 1 );
    }
    else{
      Log << "Test for allelic association written to " << options->getAllelicAssociationScoreFilename() << "\n";
      genescorestream << "structure(.Data=c(" << endl;
    }
    
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
	if(!indiv->getIndividual(i)->IsMissing(j)){
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
  if( strlen( options->getTestsForSNPsInHaplotypeOutputFilename() ) ){
    if(Lociptr->GetTotalNumberOfLoci() > Lociptr->GetNumberOfCompositeLoci()){//cannot test for SNPs in Haplotype if only simple loci
      SNPsAssociationScoreStream = new ofstream( options->getTestsForSNPsInHaplotypeOutputFilename(), ios::out );
      if( !SNPsAssociationScoreStream ){
	Log.setDisplayMode(On);
	Log << "ERROR: Couldn't open haplotypeassociationscorefile\n";
	exit( 1 );
      }
      else{
	Log << "Tests for haplotype associations written to\n" << options->getTestsForSNPsInHaplotypeOutputFilename() << "\n";
	*SNPsAssociationScoreStream << "structure(.Data=c(" << endl;
      }
    }
    else {
      op->setTestForSNPsInHaplotype(false);
      Log << "ERROR: Cannot test for haplotype associations if all loci are simple\n" << "This option will be ignored\n";
    }
  }

  InitialiseAssocScoreFile(PLabels);
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
  
  if( options->getNumberOfOutcomes() > 0 ){//no updates if no regression model
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){

      //-----------------------
      // haplotype association 
      //-----------------------  
      if( options->getTestForAllelicAssociation() ){
	//loop over simple loci within haplotypes
	if( (*Lociptr)(j)->GetNumberOfLoci() > 1 ){
	  for( int l = 0; l < (*Lociptr)(j)->GetNumberOfLoci(); l++ ){
	    CentreAndSum(1, ScoreWithinHaplotype[j][l], InfoWithinHaplotype[j][l], &(SumScoreWithinHaplotype[ j ][ l ]),
			 &(SumScore2WithinHaplotype[j][l]), &(SumInfoWithinHaplotype[ j ][ l ]));
	  }
	}
      }
      
      /*----------------------
	| Allelic association  |
	-----------------------*/
      if( (options->getTestForAllelicAssociation()  && (*Lociptr)(j)->GetNumberOfLoci() == 1) || options->getTestForSNPsInHaplotype() ){
	//if(locusObsIndicator[j]){//skip loci with no observed genotypes
	  CentreAndSum(dim_[j], LocusLinkageAlleleScore[j], LocusLinkageAlleleInfo[j],SumLocusLinkageAlleleScore[j],
		       SumLocusLinkageAlleleScore2[j],SumLocusLinkageAlleleInfo[j]); 
	  //}
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
      
      //retrieve sampled hap pair from Individual
      const int* happair = ind->getSampledHapPair(locus);
      const unsigned numStates = (*Lociptr)(locus)->GetNumberOfStates();
      const unsigned numLoci = (*Lociptr)(locus)->GetNumberOfLoci();

      // count alleles / haplotypes      
      vector<int> counts;

      // if diallelic, evaluate score for allele 2 only, otherwise evaluate score for all alleles or haplotypes

      // special case for SNP (X has size K+1)
      if( numStates == 2 ){
 	//sets X[0] to -1, 0 or 1 according to whether genotype is 11, 12 or 22
	counts.push_back((*Lociptr)(locus)->getAlleleCounts(2, happair)[0] - 1 );
	
	// general case for simple locus (X has size K + nStates)
      } else if(numLoci == 1 ){
 	for( unsigned k = 0; k < numStates; k++ ){
	  counts.push_back((*Lociptr)(locus)->getAlleleCounts(k+1, happair)[0] );
 	}
      }
      // general case for composite locus (X has size (K + nMergedHaplotypes))
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

	
	if( options->getTestForSNPsInHaplotype() ){
	  //count numbers of each haplotype
	  counts = (*Lociptr)(locus)->getHaplotypeCounts(happair);
	}
      }

      UpdateAlleleScores(LocusLinkageAlleleScore[locus],LocusLinkageAlleleInfo[locus], ind->getAdmixtureProps(), counts, 
			 YMinusEY, phi, DInvLink);

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

// ********** OUTPUT **********************************************************

void ScoreTests::Output(int iteration, const std::string * PLabels){
  PopLabels = PLabels;
  //Allelic association
  if( options->getTestForAllelicAssociation() ){
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      //case of simple locus
      if((* Lociptr)(j)->GetNumberOfLoci() == 1 )
	OutputTestsForAllelicAssociation(iteration, j, dim_[j], SumLocusLinkageAlleleScore[j], SumLocusLinkageAlleleScore2[j], 
					 SumLocusLinkageAlleleInfo[j]);
      //case of haplotype
      else
	OutputTestsForAllelicAssociation(iteration, j, (*Lociptr)(j)->GetNumberOfLoci(), SumScoreWithinHaplotype[ j ], 
					 SumScore2WithinHaplotype[ j ], SumInfoWithinHaplotype[ j ]);
      
      
    }//end j loop over comp loci
  }

  //haplotype association
  if( options->getTestForSNPsInHaplotype() ){
    OutputTestsForSNPsInHaplotype( iteration );
  }

  //ancestry association
  if( options->getTestForLinkageWithAncestry() ){
    
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

  //admixture association
  if( options->getTestForAdmixtureAssociation() ){
    OutputAdmixtureScoreTest( iteration );
  }

}

// next few function calculate score tests from the cumulative sums of
// the score, score squared, and information score and info can be
// scalars, or respectively a vector and a matrix should have just one
// method or class to do this
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

void ScoreTests::OutputTestsForSNPsInHaplotype( int iteration )
  // misleading name for this method - should be called OutputTestsForHaplotypeAssociation
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
      scale_matrix(ScoreVector, 1.0/( iteration - options->getBurnIn()), dim_[j], 1);
      
      CompleteMatrix = new double[dim_[j]*dim_[j]];
      copy(SumLocusLinkageAlleleInfo[j], SumLocusLinkageAlleleInfo[j]+ dim_[j]*dim_[j], CompleteMatrix);
      scale_matrix(CompleteMatrix, 1.0/( iteration - options->getBurnIn()), dim_[j], dim_[j]);
      
      ObservedMatrix = new double[dim_[j]*dim_[j]];
      for(unsigned d1 = 0; d1 < dim_[j]; ++d1)for(unsigned d2 = 0; d2 < dim_[j]; ++d2)
	ObservedMatrix[d1*dim_[j] + d2] = CompleteMatrix[d1*dim_[j]+d2] + ScoreVector[d1]*ScoreVector[d2] -
	  SumLocusLinkageAlleleScore2[j][d1*dim_[j]+d2]/( iteration - options->getBurnIn() );
      
      NumberOfMergedHaplotypes = dim_[j];
      for( int k = 0; k < NumberOfMergedHaplotypes; k++ ){
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
	*SNPsAssociationScoreStream  << double2R(ScoreVector[k]) << ",";
	*SNPsAssociationScoreStream  << double2R(CompleteMatrix[k*dim_[j]+k]) << ",";
	*SNPsAssociationScoreStream  << double2R(ObservedMatrix[k*dim_[j]+k]) << ",";
	*SNPsAssociationScoreStream  << double2R(100*ObservedMatrix[k*dim_[j]+k] / CompleteMatrix[k*dim_[j]+k]) << ",";
	*SNPsAssociationScoreStream  << double2R(ScoreVector[ k ] / sqrt( ObservedMatrix[k*dim_[j]+k] ));
	*SNPsAssociationScoreStream  << ",";
	// if not last allele at locus, output unquoted "NA" in chi-square column
	if( k != NumberOfMergedHaplotypes - 1 ){
	  *SNPsAssociationScoreStream  << "NA," << endl;
	}
      }
      // calculate summary chi-square statistic
      double chisq = GaussianConditionalQuadraticForm( NumberOfMergedHaplotypes - 1, ScoreVector, ObservedMatrix, dim_[j] );
      if(chisq == -1)	*SNPsAssociationScoreStream  << "NaN" << "," << endl;// -1 means ObservedMatrix is rank deficient
      else
	*SNPsAssociationScoreStream  << double2R(chisq) << "," << endl;

    }
  }
}

void ScoreTests::OutputTestsForAllelicAssociation( int iteration, int locus, unsigned dim, const double* score, const double* scoresq, 
						   const double* info)
{
  double Score, CompleteInfo, MissingInfo, ObservedInfo, PercentInfo, zscore;
  for(unsigned a = 0; a < dim; ++a){
    Score = score[a] / ( iteration - options->getBurnIn() );
    CompleteInfo = info[a] / ( iteration - options->getBurnIn() );
    MissingInfo = scoresq[a] / ( iteration - options->getBurnIn() ) - Score * Score;
    ObservedInfo = CompleteInfo - MissingInfo;
    
    if(CompleteInfo > 0.0) {
      PercentInfo = 100*ObservedInfo / CompleteInfo;
      zscore = Score / sqrt( ObservedInfo );
    } else {
      PercentInfo = 0;
      zscore = 0;
    }

    string locuslabel = (*Lociptr)(locus)->GetLabel(a);
    if(dim==1 || (*Lociptr)(locus)->GetNumberOfLoci()>1) genescorestream << locuslabel << ",";
    else genescorestream << locuslabel.substr(0, locuslabel.size()-1)<< "("<<a+1<<")\",";
    genescorestream << double2R(Score)        << ","
		    << double2R(CompleteInfo) << ","
		    << double2R(ObservedInfo) << ",";
    if(CompleteInfo > 0.0) {
      genescorestream << double2R(PercentInfo) << ",";
      genescorestream << double2R(zscore)    << "," << endl;
    }
    else{
      genescorestream << "NaN" << ","
		      << "NaN" << "," << endl;
    }
    
  }
}

void ScoreTests::OutputTestsForLocusLinkage( int iteration, ofstream* outputstream,
					     const double* Score, const double* VarScore,
					     const double* Score2, const double* Info )
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
      *outputstream << (*Lociptr)(j)->GetLabel(0) << ",";
      *outputstream << "\""<<PopLabels[k+k1] << "\","; //need offset to get second poplabel for 2pops
      
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
  int numPrintedIterations = (options->getTotalSamples() - options->getBurnIn()) / (options->getSampleEvery() * 10);
  /**
   * writes out the dimensions and labels of the 
   * R object previously written to genescorestream
   */
  int count;
  if(options->getTestForAllelicAssociation()){
    count = 0;
    for(unsigned int j = 0; j < Lociptr->GetNumberOfCompositeLoci(); j++ ){
      if((* Lociptr)(j)->GetNumberOfLoci() == 1 ) count += dim_[j];
      else count += (*Lociptr)(j)->GetNumberOfLoci();
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
