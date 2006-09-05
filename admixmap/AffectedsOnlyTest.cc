/** 
 *   ADMIXMAP
 *   AffectedsOnlyTest.cc 
 *   Class implements affecteds-only score test for linkage with locus ancestry
 *   Copyright (c) 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */

#include "AffectedsOnlyTest.h"
#include "Genome.h"
#include <gsl/gsl_cdf.h>
#include <math.h>

using namespace std;

AffectedsOnlyTest::AffectedsOnlyTest(){
  SumAffectedsScore = 0;
  SumAffectedsScore2 = 0;
  SumAffectedsInfo = 0;
  SumAffectedsVarScore = 0;
}

AffectedsOnlyTest::~AffectedsOnlyTest(){
  //delete arrays for affecteds-only score test
  delete[] SumAffectedsScore;
  delete[] SumAffectedsScore2;
  delete[] SumAffectedsInfo;
  delete[] SumAffectedsVarScore;

}

void AffectedsOnlyTest::Initialise(const char* filename, const int NumPopulations, const int NumLoci, LogWriter &Log){
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
  OpenFile(Log, &outputfile, filename, "Affected-only tests for association", true);
  
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

void AffectedsOnlyTest::Output(int iterations, const Vector_s& PopLabels, const Genome& Loci, bool final, const char* filename){
  std::ofstream* outfile;
  if(final){
    outfile = new ofstream(filename);
    *outfile <<"Locus\tPopulation\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tMissing1\tMissing2\tStdNormal\tPValue\n";
  }
  else outfile = &outputfile;

  string sep = final? "\t" : ",";
  for(unsigned int j = 0; j < L; j++ ){
    const string locuslabel = Loci(j)->GetLabel(0);
    for( unsigned k = 0; k < K; k++ ){//end at 1 for 2pops
      const std::string label = locuslabel + "\"" + sep + "\"" + PopLabels[k+firstpoplabel];//label is output in quotes
	OutputRaoBlackwellizedScoreTest(iterations, outfile, label, SumAffectedsScore[ j*K + k], SumAffectedsScore2[ j*K + k], 
					SumAffectedsVarScore[ j*K + k ],SumAffectedsInfo[ j*K + k ], final); 
    }
  }
  if(final)outfile->close();
}

/**
 * writes out the dimensions and labels of the 
 * R-matrix already written to file
 */
void AffectedsOnlyTest::ROutput(const int numIterations){
  if (test){
    std::vector<int> dimensions(3,0);
    dimensions[0] = 10;
    dimensions[1] = L * K;
    dimensions[2] = (int)(numIterations);
    
    std::vector<std::string> labels(dimensions[0],"");
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
    R_output3DarrayDimensions(&outputfile,dimensions,labels);
  }
}

/**
   Updates score, variance of score and info for a single individual at a single locus
*/
void AffectedsOnlyTest::Update(unsigned int locus, int k0, const double* const Theta, 
			       bool RandomMatingModel, bool diploid, 
			       const vector<vector<double> > AProbs){
  // values of ancestry risk ratio at which likelihood ratio is evaluated
  double r1 = 0.5;
  double r2 = 2.0;//hard-coding these for now, can make them vary later
  if( diploid ) { // diploid case
    double theta[2];//paternal and maternal admixture proportions
    double Pi[3];//probs of 0,1,2 copies of Pop k given admixture
    for( unsigned k = 0; k < K; k++ ){
      theta[0] = Theta[ k+k0 ];
      if( RandomMatingModel )
	theta[1] = Theta[ K + k+k0 ];
      else
	theta[1] = theta[0];
      
      //accumulate score, score variance, and info
      AffectedsScore[locus *K + k]+= 0.5*( AProbs[1][k+k0] + 2.0*AProbs[2][k+k0] - theta[0] - theta[1] );
      AffectedsVarScore[locus * K + k]+= 0.25 *( AProbs[1][k+k0]*(1.0 - AProbs[1][k+k0]) + 4.0*AProbs[2][k+k0]*AProbs[0][k+k0]); 
      AffectedsInfo[locus * K +k]+= 0.25* ( theta[0]*( 1.0 - theta[0] ) + theta[1]*( 1.0 - theta[1] ) );
      
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
    double theta;//paternal and maternal admixture proportions
    double Pi[2];//probs of 0,1 copies of Pop k given admixture
    for( unsigned k = 0; k < K; k++ ){
      theta = Theta[ k+k0 ];
      
      //accumulate score, score variance, and info
      AffectedsScore[locus *K + k] += AProbs[1][k+k0] - theta;
      AffectedsVarScore[locus * K + k] += AProbs[0][k+k0] * AProbs[1][k+k0]; 
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

///accumulates E(score), E(score squared), variance of score and info over iterations, for locus j
void AffectedsOnlyTest::Accumulate(unsigned j){
  for( unsigned k = 0; k < K; k++ ){
    SumAffectedsScore[j*K +k] += AffectedsScore[j*K + k];
    SumAffectedsVarScore[j*K +k] += AffectedsVarScore[j * K +k];
    SumAffectedsInfo[j*K +k] += AffectedsInfo[j * K +k];
    SumAffectedsScore2[j*K +k] +=  AffectedsScore[j*K +k] * AffectedsScore[j*K +k];
  }
}

///outputs ergodic averages of Likelihood Ratios as R object
void AffectedsOnlyTest::OutputLikRatios(const char* const filename, int iterations, const Vector_s& PopLabels, const Genome& Loci){
  //open outut file
  std::ofstream likratiostream(filename);

  //start writing R object
  likratiostream<< "structure(.Data=c(" << endl;

  double L1, L2;

  for(unsigned int j = 0; j < L; j++ ){
    for( unsigned k = 0; k < K; k++ ){//end at 1 for 2pops
      likratiostream << "\"" << Loci(j)->GetLabel(0) << "\",";
      likratiostream << "\""<<PopLabels[k+firstpoplabel] << "\","; //need offset to get second poplabel for 2pops
      
      L1 = LikRatio1[ j*K + k] / ( iterations );
      L2 = LikRatio2[ j*K + k] / ( iterations );
      
      likratiostream << double2R(L1)<< ","
		   << double2R(L2)<< ","<<endl;
    }
  }
  
  std::vector<int> dim(2,0);
  dim[0] = 4;
  dim[1] = L * K;
  
  std::vector<std::string> labels(4,"");
  labels[0] = "Locus";
  labels[1] = "Population";
  labels[2] = "L1";
  labels[3] = "L2";

  R_output3DarrayDimensions(&likratiostream, dim, labels);
  
  likratiostream.close();
}
