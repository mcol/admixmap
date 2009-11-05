/** 
 *   ADMIXMAP
 *   AffectedsOnlyTest.cc 
 *   Class implements affecteds-only score test for linkage with locus ancestry
 *   Copyright (c) 2006, 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "AffectedsOnlyTest.h"
#include "Genome.h"
#include "bclib/DelimitedFileWriter.h"
#include "bclib/LogWriter.h"
#include <gsl/gsl_cdf.h>
#include <math.h>

using namespace std;

AffectedsOnlyTest::AffectedsOnlyTest(){
  SumAffectedsScore = 0;
  SumAffectedsScore2 = 0;
  SumAffectedsInfo = 0;
  SumAffectedsVarScore = 0;
  test=false;
  numPrintedIterations = 0;
  numUpdates = 0;
}

AffectedsOnlyTest::~AffectedsOnlyTest(){
  //delete arrays for affecteds-only score test
  delete[] SumAffectedsScore;
  delete[] SumAffectedsScore2;
  delete[] SumAffectedsInfo;
  delete[] SumAffectedsVarScore;

  if (test){

    delete[] AffectedsScore;
    delete[] AffectedsVarScore;
    delete[] AffectedsInfo;
    delete[] LikRatio1;
    delete[] LikRatio2;

    std::vector<std::vector<std::string> > labels(1);
    labels[0].push_back("Locus");
    labels[0].push_back("Population");
    labels[0].push_back("minusLogPValue");
    
    std::vector<int> dimensions(3,0);
    dimensions[0] = labels[0].size(); 
    dimensions[1] = L * K;
    dimensions[2] = (int)(numPrintedIterations);

    R.close(dimensions,labels);
  }
}

void AffectedsOnlyTest::Initialise(const char* filename, const int NumPopulations, const int NumLoci){
  test = true;
  //K is the number of populations for which to perform test. For 2-way admixture, we only want 2nd population.
  K = NumPopulations;
  firstpoplabel = 0;
  if(NumPopulations == 2){
    K = 1;
    firstpoplabel = 1;//skip first pop label if 2 pops
  }

  L = NumLoci;

  //open output file and start writing R object
  R.open(filename);
  
  const int dim = L*K;
  AffectedsScore = new double[dim];
  AffectedsVarScore = new double[dim];
  AffectedsInfo = new double[dim];
  LikRatio1 = new double[dim];
  LikRatio2 = new double[dim];
  fill(LikRatio1, LikRatio1+dim, 0.0);
  fill(LikRatio2, LikRatio2+dim, 0.0);

  SumAffectedsScore2 = new double[dim];
  SumAffectedsScore = new double[dim];
  SumAffectedsVarScore = new double[dim];
  SumAffectedsInfo = new double[dim];
  fill(SumAffectedsScore, SumAffectedsScore + dim, 0.0);
  fill(SumAffectedsScore2, SumAffectedsScore2 + dim, 0.0);
  fill(SumAffectedsInfo, SumAffectedsInfo + dim, 0.0);
  fill(SumAffectedsVarScore, SumAffectedsVarScore + dim, 0.0);
}

///resets arrays holding sums of scores and info over individuals to zero; invoked at start of each iteration after burnin.
void AffectedsOnlyTest::Reset(){
  if( test ){
    for(unsigned j = 0; j < L*K; ++j){
      AffectedsScore[j] = 0.0;
      AffectedsVarScore[j] = 0.0;
      AffectedsInfo[j] = 0.0;
    }
  }
}

void AffectedsOnlyTest::Output(const Vector_s& PopLabels, const Genome& Loci){
  OutputAffectedsOnlyTest(R, PopLabels, Loci, ",", false);
    ++numPrintedIterations;
}

void AffectedsOnlyTest::WriteFinalTable(const char* filename, const Vector_s& PopLabels, 
					const Genome& Loci, bclib::LogWriter& Log){
  bclib::DelimitedFileWriter finaltable(filename);
  Log << bclib::Quiet << "Affected-only tests for association written to " << filename << "\n";
  finaltable <<"Locus" << "Population" << "Score" << "CompleteInfo" << "ObservedInfo" 
	     << "PercentInfo" << "Missing1" << "Missing2" << "StdNormal" << "PValue"
	     << bclib::newline;

  OutputAffectedsOnlyTest(finaltable, PopLabels, Loci, "\t", true);
  finaltable.close();
}
void AffectedsOnlyTest::OutputAffectedsOnlyTest(bclib::DelimitedFileWriter& outfile, const Vector_s& PopLabels, 
						const Genome& Loci, const string& sep, bool final){
  for(unsigned int j = 0; j < L; j++ ){
    const string locuslabel = Loci(j)->GetLabel(0);
    for( unsigned k = 0; k < K; k++ ){//end at 1 for 2pops
      const std::string label = locuslabel + "\"" + sep + "\"" + PopLabels[k+firstpoplabel];//label is output in quotes
	OutputRaoBlackwellizedScoreTest(outfile, label, SumAffectedsScore[ j*K + k], SumAffectedsScore2[ j*K + k], 
					SumAffectedsVarScore[ j*K + k ],SumAffectedsInfo[ j*K + k ], final); 
    }
  }

}
/**
   Updates score, variance of score and info for a single individual at a single locus
*/
void AffectedsOnlyTest::Update(unsigned int locus, int k0, const double* const Theta, 
			       bool RandomMatingModel, bool diploid, 
			       const vector<vector<double> >& AProbs){
  // k0 should be passed as 1 if K==1 
  // values of ancestry risk ratio at which likelihood ratio is evaluated
  double r1 = 0.5;
  double r2 = 2.0;//hard-coding these for now, can make them vary later
  if( diploid ) { // diploid case
    double theta[2];//paternal and maternal admixture proportions
    double Pi[3];//probs of 0,1,2 copies of Pop k given admixture
    for( unsigned k = 0; k < K; k++ ){ //if K==1, only k=0 will be evaluated
      theta[0] = Theta[ k+k0 ];
      if( RandomMatingModel )
	if( K==1 ) {  // we want theta[1] to be set to Theta[3]
	  theta[1] = Theta[ 2 + k+k0 ];
	} else {
	  theta[1] = Theta[ K + k+k0 ];
	}
      else
	theta[1] = theta[0];
      
#define DEBUG_AOTEST 0
#if DEBUG_AOTEST
  fprintf( stderr, "AProbs: %zd\n", AProbs.size() );
  for ( size_t idx1 = 0 ; idx1 < AProbs.size() ; ++idx1 )
    {
    const vector<double> & apRow = AProbs[ idx1 ];
    fprintf( stderr, "%3zd[%zd] ", idx1, apRow.size() );
    for ( size_t idx2 = 0 ; idx2 < apRow.size() ; ++idx2 )
	fprintf( stderr, " %.9lf", apRow[idx2] );
    putc( '\n', stderr );
    }
  fprintf(stderr, "AOT-debug-2: t=%d,k=%d,k0=%d mu[0]=%.12lf mu[1]=%.12lf", locus, k, k0, theta[0], theta[1] );
  fprintf(stderr, " old cum-score=%.12lf", AffectedsScore[locus*K+k] );
#endif
      //accumulate score, score variance, and info
      AffectedsScore[locus *K + k]+= 0.5*( AProbs[1][k+k0] + 2.0*AProbs[2][k+k0] - theta[0] - theta[1] );
      AffectedsVarScore[locus * K + k]+= 0.25 *( AProbs[1][k+k0]*(1.0 - AProbs[1][k+k0]) + 4.0*AProbs[2][k+k0]*AProbs[0][k+k0]); 
      AffectedsInfo[locus * K +k]+= 0.25* ( theta[0]*( 1.0 - theta[0] ) + theta[1]*( 1.0 - theta[1] ) );
      
#if DEBUG_AOTEST
  fprintf(stderr, " cum-score=%.12lf cum-var=%.12lf cum-inf=%.12lf\n", AffectedsScore[locus*K+k],
    AffectedsVarScore[locus*K+k], AffectedsInfo[locus*K+k] );
#endif
      //probs of 0,1,2 copies of Pop k given admixture
      Pi[2] = theta[0] * theta[1];
      Pi[1] = theta[0] * (1.0 - theta[1]) + theta[1] * (1.0 - theta[0]);
      Pi[0] = (1.0 - theta[0]) * (1.0 - theta[1]);
      
      //compute contribution to likelihood ratio
      LikRatio1[locus *K + k] += (AProbs[0][k+k0] + sqrt(r1)*AProbs[1][k+k0] + r1 * AProbs[2][k+k0]) / 
	(Pi[0] + sqrt(r1)*Pi[1] + r1*Pi[2]);
      LikRatio2[locus *K + k] += (AProbs[0][k+k0] + sqrt(r2)*AProbs[1][k+k0] + r2 * AProbs[2][k+k0]) / 
	(Pi[0] + sqrt(r2)*Pi[1] + r2*Pi[2]);
    }
  } else { // haploid - effect of one extra copy from pop k0 is equivalent to two extra copies in diploid case
    double theta;//gamete admixture proportions
    // should call with maternal admixture proportions if random mating and X chr in male
    double Pi[2];//probs of 0,1 copies of Pop k given admixture
    for( unsigned k = 0; k < K; k++ ){
      theta = Theta[ k+k0 ]; 
      
      //accumulate score, score variance, and info
      AffectedsScore[locus *K + k] += AProbs[1][k+k0] - theta;
      AffectedsVarScore[locus * K + k] += AProbs[1][k+k0] * (1 - AProbs[1][k+k0]); 
      AffectedsInfo[locus * K +k]+= theta * (1.0 - theta);
      
      //probs of 0,1 copies of Pop k given admixture
      Pi[1] = theta;
      Pi[0] = 1.0 - theta;
      
      //compute contribution to likelihood ratio - check this formula
      LikRatio1[locus *K + k] += (AProbs[0][k+k0] + r1*AProbs[1][k+k0]) / (Pi[0] + r1*Pi[1]);
      LikRatio2[locus *K + k] += (AProbs[0][k+k0] + r2*AProbs[1][k+k0]) / (Pi[0] + r2*Pi[1]);
    }
  }
}

///accumulates E(score), E(score squared), variance of score and info over iterations
void AffectedsOnlyTest::Accumulate(){
  for(unsigned j = 0; j < L; ++j)
    for( unsigned k = 0; k < K; k++ ){
      SumAffectedsScore[j*K +k] += AffectedsScore[j*K + k];
      SumAffectedsVarScore[j*K +k] += AffectedsVarScore[j * K +k];
      SumAffectedsInfo[j*K +k] += AffectedsInfo[j * K +k];
      SumAffectedsScore2[j*K +k] +=  AffectedsScore[j*K +k] * AffectedsScore[j*K +k];
    }

  //increment update counter
  ++numUpdates;
}

///outputs ergodic averages of Likelihood Ratios as R object
void AffectedsOnlyTest::OutputLikRatios(const char* const filename, const Vector_s& PopLabels, const Genome& Loci){
  //open outut file
  bclib::RObjectWriter likratiostream(filename);

  double L1, L2;

  for(unsigned int j = 0; j < L; j++ ){
    const string locuslabel = "\"" + Loci(j)->GetLabel(0) + "\"";
    for( unsigned k = 0; k < K; k++ ){//end at 1 for 2pops
      const string poplabel = "\"" + PopLabels[k+firstpoplabel] + "\"";
      likratiostream << locuslabel
		     << poplabel; //need offset to get second poplabel for 2pops
      
      L1 = LikRatio1[ j*K + k] / ( numUpdates );
      L2 = LikRatio2[ j*K + k] / ( numUpdates );
      
      likratiostream << double2R(L1)
		     << double2R(L2)<< bclib::newline;
    }
  }
  
  std::vector<int> dim(2,0);
  dim[0] = 4;
  dim[1] = L * K;
  
  std::vector<std::vector<std::string> > labels(1);
  labels[0].push_back("Locus");
  labels[0].push_back("Population");
  labels[0].push_back("L1");
  labels[0].push_back("L2");

  likratiostream.close(dim, labels);
}
