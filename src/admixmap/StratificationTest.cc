/** 
 *   ADMIXMAP
 *   StratificationTest.cc 
 *   Class to implement a test for residual population stratification
 *   Copyright (c) 2005, 2006, 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "StratificationTest.h"
#include "AdmixOptions.h"
#include "AdmixFilenames.h"
#include <numeric>
#include "gsl/gsl_eigen.h"
#include "gsl/gsl_linalg.h"
#include "gsl/gsl_blas.h"
#include "FreqArrays.h"
#include "bcppcl/LogWriter.h"

using namespace std;

StratificationTest::StratificationTest()
{
   NumberOfTestLoci = 0;
   count = 0;
   T = 0;
}
void StratificationTest::Initialize( AdmixOptions* const options, const Genome &Loci,  
				     const IndividualCollection* const IC, LogWriter &Log )
{
  Log.setDisplayMode(Quiet);
  if(options->getStratificationTest() ){

    //float DistanceFromLast = 0;
    //OLD CODE
//     for(unsigned int j = 0; j < Loci.GetNumberOfCompositeLoci(); j++ ){
//       DistanceFromLast += Loci.GetDistance(j);
//       if( DistanceFromLast > 10 ){ // set cutoff to 10 morgans so that only unlinked loci will be included
// 	// should have a more specific way to code unlinked loci than setting Distance=100
// 	// this algorithm selects the first locus on each chromosome
// 	// preferable to select the most informative locus (greatest heterozygosity)
// 	if( Loci(j)->GetNumberOfStates() == 2 ){ // test uses only diallelic loci
// 	  // should be able to use multi-allelic loci by grouping alleles into 2 bins
// 	  TestLoci.push_back(j);
// 	  NumberOfTestLoci++;
// 	  DistanceFromLast = 0;
// 	}
//       }
//     }
    //NEW CODE
    unsigned abslocus = 0;
    for(unsigned c = 0; c < Loci.GetNumberOfChromosomes(); ++c){
      int j = -1;
      double max = 0;
      double n1, n2;
      double ExpHet;//expected heterozygosity
      //select most informative locus (with greatest exp heterozygosity), from those with
      //less than 5% missing genotypes, on each chromosome
      for(unsigned locus = 0; locus < Loci.GetSizeOfChromosome(c); ++locus){
	if(Loci(abslocus)->GetNumberOfStates() == 2){// test uses only diallelic loci

 	  //count number of missing genotypes at locus
	  int count = IC->getNumberOfMissingGenotypes(locus);

	  if( ( (double)count / (double)(IC->getSize()) ) < 0.1){//exclude loci with >10% missing genotypes
	    n1 = (double) IC->GetSNPAlleleCounts(locus, 1);//# copies of allele1
	    n2 = (double) IC->GetSNPAlleleCounts(locus, 2);//# copies of allele2
	    ExpHet = 2.0 * (n1*n2) / ((n1+n2)*(n1+n2)); 
	    if( ExpHet > max){
	      max  = ExpHet;
	      j = abslocus;//gets locus number (on genome) of this locus
	    }
	  }
	}
	++abslocus;
      }
      //j is the most informative locus
      if(j > -1 ){
	TestLoci.push_back(j);
	NumberOfTestLoci++;
	//DistanceFromLast = 0;
      }
    }
    
    if( NumberOfTestLoci < 2 ){
      Log.setDisplayMode(On);
      Log << "\nToo few unlinked loci to run stratification test\n";
      options->setStratificationTest(false);
      if(outputstream.is_open())outputstream.close();
    }
    else{
      Log.setDisplayMode(Off);
      Log << NumberOfTestLoci << " loci used in stratification test.\n";
      for(int i = 0; i < NumberOfTestLoci; i++){
	Log << Loci(TestLoci[i])->GetLabel(0) << "\n";
      }
      ModelIndicator = options->isRandomMatingModel();
      
      outputstream << "T_obs" << "\t" << "T_rep\n";
    }
    
  }
  else{
    Log << "No test for residual population stratification.\n";
  }
}

void StratificationTest::calculate( const IndividualCollection* const individuals, const FreqArray& AlleleFreqs,
				    const vector<vector<int> > ChrmAndLocus, int Populations )
{
  // matrix of (observed minus expected copies allele 1) scores for each individual at each locus
  gsl_matrix *popX = gsl_matrix_calloc( individuals->getSize(), NumberOfTestLoci ); 
  // matrix of replicate scores
  gsl_matrix *popRepX = gsl_matrix_calloc( individuals->getSize(), NumberOfTestLoci );
  //bool flag = false;
  vector<unsigned short> genotype(2, 0);

  for( int j = 0; j < NumberOfTestLoci; j++ ){
    int jj = TestLoci[j];

    //const double* const freqs = AlleleFreqs[jj];  // array of length (NumberOfStates-1)*Populations

    for( int i = 0; i < individuals->getSize(); i++ ){
      const Individual* const ind = individuals->getIndividual(i);
      bool diploid = !(ind->isHaploidatLocus(jj));
      if(diploid){//skip Xloci in males
	if(ind->GenotypeIsMissing(jj)){// if genotype is missing, sample it
	  int ancestry[2];
	  ind->GetLocusAncestry( ChrmAndLocus[jj][0], ChrmAndLocus[jj][1], ancestry );
	  genotype = SimGenotypeConditionalOnAncestry(AlleleFreqs[jj] , ancestry );
	}
	else{//get sampled haplotype pair
	  // (the number of copies of allele1 is the same as in the observed genotype)
	  const int* genotypeArray = ind->getSampledHapPair(jj);
	  genotype[0] = genotypeArray[0]+1;
	  genotype[1] = genotypeArray[1]+1;
	}
	//if( genotype[0] != genotype[1] ){ // if heterozygous
	// genotype = SampleHeterozygotePhase( freqs, ancestry ); // sample phase conditional on ordered diploid ancestry 
	//} else 
	
	// ProbAllele1 = Prob( allele 1 ) conditional on individual admixture
	vector<double> ProbAllele1 = GenerateExpectedGenotype( ind, AlleleFreqs[jj], Populations );
	vector<unsigned short> repgenotype = SimGenotypeConditionalOnAdmixture( ProbAllele1 );
	// vector<unsigned short> repgenotype = SimGenotypeConditionalOnAncestry( freqs, ancestry );
	// Calculate score X = Obs0 + Obs1 - Expected0 - Expected1, where Obs0, Obs1 are coded 1 for allele 1, 0 for allele 2
	// Obs0 + Obs1 = 4 - genotype[0] - genotype[1]
	gsl_matrix_set(popX, i, j, (double)(4 - genotype[0] - genotype[1]) - ProbAllele1[0] - ProbAllele1[1]);
	gsl_matrix_set(popRepX, i, j ,(double)(4 - repgenotype[0] - repgenotype[1]) - ProbAllele1[0] - ProbAllele1[1]);
      }//end is diploid
    }//end indiv loop  
  }//end locus loop
  
  // covariance matrix for (observed minus expected copies allele 1) scores
  gsl_matrix *Cov = gsl_matrix_calloc(NumberOfTestLoci, NumberOfTestLoci);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, popX, popX, 0, Cov  );//Cov = popX' * popX

  // covariance matrix for replicate scores
  gsl_matrix *RepCov = gsl_matrix_calloc(NumberOfTestLoci, NumberOfTestLoci);
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, popRepX, popRepX, 0, RepCov  );//RepCov = popRepX' * popRepX

  //compute eigenvalues
  gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc (NumberOfTestLoci);//allocates workspace for eigenvalue calculations
  gsl_vector *EigenValues = gsl_vector_calloc(NumberOfTestLoci);
  gsl_eigen_symm (Cov, EigenValues, w);

  gsl_vector *RepEigenValues = gsl_vector_calloc(NumberOfTestLoci);
  gsl_eigen_symm(RepCov, RepEigenValues, w);

  gsl_eigen_symm_free(w); 

  //normalise eigenvalues
  const double Sum = accumulate(EigenValues->data, EigenValues->data + NumberOfTestLoci, 0.0);
  gsl_vector_scale (EigenValues, 1.0/Sum) ;
  const double RepSum = accumulate(RepEigenValues->data, RepEigenValues->data + NumberOfTestLoci, 0.0);
  gsl_vector_scale (RepEigenValues, 1.0/RepSum) ;
  //find max eigenvalues
  double maxEigenValue = gsl_vector_max(EigenValues);
  double maxRepEigenValue = gsl_vector_max(RepEigenValues);
  //output to file
  outputstream << maxEigenValue << "\t" << maxRepEigenValue << endl;
  //increment test statistic
  if( maxEigenValue < maxRepEigenValue ){
    T++;
  }
  count++;
  //clean up
  gsl_vector_free(EigenValues);
  gsl_vector_free(RepEigenValues);
  gsl_matrix_free(Cov);
  gsl_matrix_free(RepCov);
  gsl_matrix_free(popX);
  gsl_matrix_free(popRepX);
}

vector<double> StratificationTest::GenerateExpectedGenotype( const Individual* const ind, const double* freqs, const int Populations )
{
  vector<double> pA(2,0);
  for( int k = 0; k < Populations; k++ ){
    pA[0] += freqs[ k*2 ] * ind->getAdmixtureProps()[k];
    if( ModelIndicator )
      pA[1] += freqs[ k*2 ] * ind->getAdmixtureProps()[k + Populations];
    else
      pA[1] += freqs[ k*2 ] * ind->getAdmixtureProps()[k];       
  }
  return pA; // vector of length 2 specifying probs of allele 1 on each gamete
}

vector<unsigned short> StratificationTest::SimGenotypeConditionalOnAdmixture( const vector<double> ProbAllele1 )
{
  vector<unsigned short> repgenotype(2,0);
  if( ProbAllele1[0] > Rand::myrand() )
    repgenotype[0] = 1;
  else
    repgenotype[0] = 2;
  if( ProbAllele1[1] > Rand::myrand() )
    repgenotype[1] = 1;
  else
    repgenotype[1] = 2;
  return repgenotype;
}

vector<unsigned short> StratificationTest::SimGenotypeConditionalOnAncestry( const double* const freqs, const int ancestry[2] )
{
  vector<unsigned short> repgenotype(2,0);
  if( freqs[ ancestry[0]*2 ] > Rand::myrand() )
     repgenotype[0] = 1;
  else
    repgenotype[0] = 2;
  if( freqs[ ancestry[1]*2 ] > Rand::myrand() )
    repgenotype[1] = 1;
  else
    repgenotype[1] = 2;
  return repgenotype;
}

vector<unsigned short> StratificationTest::SampleHeterozygotePhase( const double* freqs, const int Ancestry[2] )
{
  vector<unsigned short> genotype(2, 0);
  double q1 = freqs[ Ancestry[0]*2 ] * ( 1 - freqs[ Ancestry[1]*2 ] );
  double q2 = freqs[ Ancestry[1]*2 ] * ( 1 - freqs[ Ancestry[0]*2 ] );
  if( Rand::myrand() > q1 / ( q1 + q2 ) ){
    genotype[0] = 2;
    genotype[1] = 1;
  }
  else{
    genotype[0] = 1;
    genotype[1] = 2;
  }
  return genotype;
}

void StratificationTest::OpenOutputFile( const std::string& ResultsDir, LogWriter &Log){
  Log.setDisplayMode(Quiet);
  std::string filename = ResultsDir + "/" + STRAT_TEST_FILE;
  outputstream.open(filename.c_str(), ios::out );
  if( !outputstream.is_open() ){
    throw string("ERROR: Couldn't open stratificationtestfile");
  }
  Log << "Writing results of test for residual population stratification to "
      << filename << "\n";
}

void StratificationTest::Output(LogWriter &Log){
  Log.setDisplayMode(Quiet);
  Log << "\nStratification test: posterior predictive check probability " << (float)T/count << "\n\n";
}
