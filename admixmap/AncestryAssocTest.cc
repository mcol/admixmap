/** 
 *   ADMIXMAP
 *   AncestryAssocTest.h 
 *   Class to implement score test for association of trait with ancestry
 *   Copyright (c) 2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "AncestryAssocTest.h"
#include "Genome.h"
#include "utils/linalg.h"
#include "utils/dist.h"
#include "utils/misc.h"

using namespace::std;
AncestryAssocTest::AncestryAssocTest(){
  SumAncestryScore = 0;
  SumAncestryInfo = 0;
  SumAncestryVarScore = 0;
  SumAncestryScore2 = 0;
  AncestryScore = 0;
  AncestryInfo = 0;
  AncestryVarScore = 0;
  AncestryInfoCorrection = 0;
  B = 0;
  PrevB = 0;
  Xcov = 0;
  K = 0;
  L = 0;
  numPrintedIterations = 0;
}
AncestryAssocTest::~AncestryAssocTest(){
  if(test){
    delete[] SumAncestryScore;
    delete[] SumAncestryInfo;
    delete[] SumAncestryVarScore;
    delete[] SumAncestryScore2;
    delete[] B;
    delete[] PrevB;
    delete[] Xcov;
    free_matrix(AncestryScore, L);
    free_matrix(AncestryInfo, L);  
    free_matrix(AncestryVarScore, L);
    free_matrix(AncestryInfoCorrection, L);
  }
}
void AncestryAssocTest::Initialise(const char* filename, const int NumPopulations, const int NumLoci, LogWriter &Log){
  test = true;
  OpenFile(Log, &outputfile, filename, "Tests for locus linkage", true);
  
  K = NumPopulations;
  L = NumLoci;
  int KK = NumPopulations;
  firstpoplabel = 0;
  if(K == 2){
    KK = 1;//only keeping scores for second population when 2 populations
    firstpoplabel = 1;//skip first pop label if 2 pops
  }
  SumAncestryScore = new double[L * KK];
  SumAncestryInfo = new double[L * KK];
  SumAncestryScore2 = new double[L * KK];
  SumAncestryVarScore = new double[L * KK];
  fill(SumAncestryScore, SumAncestryScore +L*KK, 0.0);
  fill(SumAncestryScore2, SumAncestryScore2 +L*KK, 0.0);
  fill(SumAncestryInfo, SumAncestryInfo + L*KK, 0.0);
  fill(SumAncestryVarScore, SumAncestryVarScore + L*KK, 0.0);

  AncestryScore = alloc2D_d(L, 2*K);
  AncestryInfo = alloc2D_d(L, 4*K*K);
  AncestryVarScore = alloc2D_d(L, K);
  AncestryInfoCorrection = alloc2D_d(L, K);
  B = new double[K * K];
  PrevB = new double[K * K];
  Xcov = new double[K];
}

void AncestryAssocTest::Reset(){
  if(test){
    for(unsigned j = 0; j < L; ++j){
      for(unsigned k = 0; k < 2*K; ++k)
	AncestryScore[j][k] = 0.0;
      for(unsigned k = 0; k < 4*K*K; ++k)
	AncestryInfo[j][k] = 0.0;
      for(unsigned k = 0; k < K; ++k){
	AncestryInfoCorrection[j][k] = 0.0;
	AncestryVarScore[j][k] = 0.0;
      }
    }
    for(unsigned k = 0; k < K*K; ++k){
      PrevB[k] = B[k];           //PrevB stores the sum for the previous iteration
      B[k] = 0.0;                //while B accumulates the sum for the current iteration 
    }
    for(unsigned k = 0; k < K; ++k){
      Xcov[k] = 0.0;
    }
  }
}
void AncestryAssocTest::Output(int iterations, const Vector_s& PopLabels, const Genome& Loci, bool final, const char* filename){
  std::ofstream* outfile;
  if(final){
    outfile = new ofstream(filename);
    *outfile <<"Locus\tPopulation\tScore\tCompleteInfo\tObservedInfo\tPercentInfo\tMissing1\tMissing2\tStdNormal\tPValue\n";
  }
  else {
    outfile = &outputfile;
    ++numPrintedIterations;
  }

  unsigned KK = K;
  if(K == 2)KK = 1;//only keeping scores for second population when 2 populations
  string sep = final? "\t" : ",";
  for(unsigned int j = 0; j < L; j++ ){
    const string locuslabel = Loci(j)->GetLabel(0);
    for( unsigned k = 0; k < KK; k++ ){//end at 1 for 2pops
      const std::string label = locuslabel + "\"" + sep + "\"" + PopLabels[k+firstpoplabel];//label is output in quotes
      OutputRaoBlackwellizedScoreTest(iterations, outfile, label, SumAncestryScore[ j*KK + k], SumAncestryScore2[ j*KK + k], 
 				      SumAncestryVarScore[ j*KK + k ],SumAncestryInfo[ j*KK + k ], final); 
    }
  }
  if(final)outfile->close();
}

void AncestryAssocTest::ROutput(const int ){
  if(test){
    int KK = K;
    if(KK ==2 )KK = 1;
    vector<int> dimensions(3,0);
    dimensions[0] = 10;
    dimensions[1] = L * KK;
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
    
    R_output3DarrayDimensions(&outputfile, dimensions, labels);
  }
}

void AncestryAssocTest::Update(int locus, const double* admixtureCovars, double phi, double YMinusEY, double DInvLink, 
			       const vector<vector<double> > AProbs) {
  //Updates score stats for test for association with locus ancestry
  //now use Rao-Blackwellized estimator by replacing realized ancestries with their expectations
  //Notes: 1/phi is dispersion parameter
  //       = lambda[0] for linear regression, = 1 for logistic
  //       YMinusEY = Y - E(Y) = Y - g^{-1}(\eta_i)
  //       VarX = Var(X)
  //       DInvLink = {d  g^{-1}(\eta)} / d\eta = derivative of inverse-link function
  // admixtureCovars are the centred admixture proportions (except first one) used in regression model
  //Xcov is a vector of covariates
  //Note that only the intercept and admixture proportions are used.
  // X is (A, cov)'  
  
  double *X = new double [2 * K], *Xcopy = new double[2*K], *XX = new double[4*K*K];
  //Xcopy is an exact copy of X; We need two copies as one will be destroyed
  double xBx[1];
  double* BX = new double[K];
  double* VarA = new double[K];
 
  X[ 2*K - 1] = 1;//intercept
  Xcov[K-1] = 1;
  //set covariates, admixture props for pops 2 to K 
  for( unsigned k = 0; k < K - 1; k++ ){
    X[ K + k] = admixtureCovars[k];//Theta[ k+1 ];
    BX[k] = Xcov[k] = admixtureCovars[k];//Theta[ k+1 ];
  }

  for( unsigned k = 0; k < K ; k++ ){
    Xcopy[k] = X[k] = AProbs[1][k] + 2.0 * AProbs[2][k];//Conditional expectation of ancestry
    VarA[k] = AProbs[1][k]*(1.0 - AProbs[1][k]) + 4.0*AProbs[2][k]*AProbs[0][k];//conditional variances
  }
  //KLUDGE: need to reset Xcopy each time since destroyed in computation of score
  Xcopy[2*K-1] = 1;
  for( unsigned k = 0; k < K-1; k++ )Xcopy[k + K] = admixtureCovars[k];//Theta[ k+1 ];

  try{
    // ** compute expectation of score **
    scale_matrix(Xcopy, YMinusEY*phi, 2*K, 1);      //Xcopy *= YMinusEY *phi
    add_matrix(AncestryScore[locus], Xcopy, 2*K, 1);//AncestryScore[locus] += Xcopy
    // ** compute uncorrected info **
    matrix_product(X, X, XX, 2*K, 1, 2*K);        //XX = X'X
    scale_matrix(XX, DInvLink*phi, 2*K, 2*K);     //XX = DInvLink * phi * X'X
    add_matrix(AncestryInfo[locus], XX, 2*K, 2*K);//AncestryInfo[locus] += XX
    // ** compute variance of score and correction term for info **    
    HH_solve(K, PrevB, Xcov, BX);          //BX = inv(PrevB) * Xcov
    matrix_product(Xcov, BX, xBx, 1, K, 1);//xBx = Xcov' * BX
  }
  catch(string s){
    delete[] X; delete[] Xcopy;
    delete[] XX;
    delete[] BX;
    delete[] VarA;
    std::string error_string = "Error in AncestryAssocTest::Update:\n";
    error_string.append(s);
    throw(error_string);//throw error message to top level
  }
  for(unsigned k = 0; k < K ; k++ ){
    AncestryInfoCorrection[locus][k] += VarA[k] * (DInvLink *phi - phi * phi * DInvLink * DInvLink * xBx[0]); 
    AncestryVarScore[locus][k] += VarA[k] * phi * phi * YMinusEY * YMinusEY;
  }
  delete[] X; delete[] Xcopy;
  delete[] XX;
  delete[] BX;
  delete[] VarA;
}

void AncestryAssocTest::Accumulate(unsigned j){
  double *score = new double[K], *info = new double[K*K];
  try{
    CentredGaussianConditional(K, AncestryScore[j], AncestryInfo[j], score, info, 2*K );
  }
  catch(string s){
    string error = "Error centring ancestry association scores:\n";
    error.append(s);
    throw(error);
  }
  
  //accumulate over iterations
  //for two populations, accumulate the scores for second population only
  unsigned KK = K, k1 = 0;
  if(K == 2){KK = 1; k1 = 1;}
     
  for( unsigned k = 0; k < KK ; k++ ){
    SumAncestryScore[j*KK +k] += score[k+k1];
    SumAncestryInfo[j*KK +k] += info[(k+k1)*K +k+k1] + AncestryInfoCorrection[j][k+k1];
    SumAncestryScore2[j*KK +k] += score[k+k1] * score[k+k1];
    SumAncestryVarScore[j*KK +k] += AncestryVarScore[j][k+k1];
  }
  delete[] score;
  delete[] info;

}

//TODO: fix to be called within update
void AncestryAssocTest::UpdateB(double DInvLink, double dispersion, const double* admixtureCovars){
  //increment B using new Admixture Props
  //Xcov is a vector of admixture props as covariates as in UpdateScoreForAncestry
    Xcov[K-1] = 1;//last entry is intercept
    for( unsigned k = 0; k < K - 1; k++ ){
      Xcov[k] = admixtureCovars[k]; //centred admixture covariates
    }
    double *temp = new double[K*K];
    matrix_product(Xcov, Xcov, temp, K, 1, K);
    scale_matrix(temp, DInvLink*dispersion, K, K);
    add_matrix(B, temp, K, K);
    delete[] temp;
}
