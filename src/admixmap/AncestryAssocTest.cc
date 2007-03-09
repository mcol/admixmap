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
  SumScore = 0;
  SumInfo = 0;
  SumVarScore = 0;
  SumScore2 = 0;
  Score = 0;
  Info = 0;
  VarScore = 0;
  InfoCorrection = 0;
  B = 0;
  PrevB = 0;
  useprevb = true;
  Xcov = 0;
  K = 0;
  L = 0;
  numPrintedIterations = 0;
}
AncestryAssocTest::AncestryAssocTest(bool use_prevb){
  AncestryAssocTest();
  useprevb = use_prevb;
}
AncestryAssocTest::~AncestryAssocTest(){
  if(test){
    delete[] SumScore;
    delete[] SumInfo;
    delete[] SumVarScore;
    delete[] SumScore2;
    delete[] B;
    if(useprevb)delete[] PrevB;
    delete[] Xcov;
    free_matrix(Score, L);
    free_matrix(Info, L);  
    free_matrix(VarScore, L);
    free_matrix(InfoCorrection, L);
  }
}
void AncestryAssocTest::Initialise(const char* filename, const int NumStrata, const int NumLoci, LogWriter &Log, bool use_prevb){
  test = true;
  OpenFile(Log, &outputfile, filename, "Tests for locus linkage", true);
  
  K = NumStrata;
  L = NumLoci;
  int KK = K;
  firstpoplabel = 0;
  if(K == 2){
    KK = 1;//only keeping scores for second population when 2 populations
    firstpoplabel = 1;//skip first pop label if 2 pops
  }
  SumScore = new double[L * KK];
  SumInfo = new double[L * KK];
  SumScore2 = new double[L * KK];
  SumVarScore = new double[L * KK];
  fill(SumScore, SumScore +L*KK, 0.0);
  fill(SumScore2, SumScore2 +L*KK, 0.0);
  fill(SumInfo, SumInfo + L*KK, 0.0);
  fill(SumVarScore, SumVarScore + L*KK, 0.0);

  Score = alloc2D_d(L, 2*K);
  Info = alloc2D_d(L, 4*K*K);
  VarScore = alloc2D_d(L, K);
  InfoCorrection = alloc2D_d(L, K);
  B = new double[K * K];
  useprevb = use_prevb;
  if(useprevb)PrevB = new double[K * K];
  else PrevB = B;
  Xcov = new double[K];
}

void AncestryAssocTest::Reset(){
  if(test){
    for(unsigned j = 0; j < L; ++j){
      for(unsigned k = 0; k < 2*K; ++k)
	Score[j][k] = 0.0;
      for(unsigned k = 0; k < 4*K*K; ++k)
	Info[j][k] = 0.0;
      for(unsigned k = 0; k < K; ++k){
	InfoCorrection[j][k] = 0.0;
	VarScore[j][k] = 0.0;
      }
    }
    for(unsigned k = 0; k < K*K; ++k){
      if(useprevb)PrevB[k] = B[k];           //PrevB stores the sum for the previous iteration
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
    *outfile <<"Locus\t";
    if(K>1)*outfile << "Population\t";
    *outfile << "Score\tCompleteInfo\tObservedInfo\tPercentInfo\tMissing1\tMissing2\tStdNormal\tPValue\n";
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
      std::string label = locuslabel;
      if(K>1)label += "\"" + sep + "\"" + PopLabels[k+firstpoplabel];//label is output in quotes
      OutputRaoBlackwellizedScoreTest(iterations, outfile, label, SumScore[ j*KK + k], SumScore2[ j*KK + k], 
 				      SumVarScore[ j*KK + k ],SumInfo[ j*KK + k ], final); 
    }
  }
  if(final)outfile->close();
}

void AncestryAssocTest::ROutput(const int ){
  if(test){
    int KK = K;
    if(KK ==2 )KK = 1;
    
    vector<string> labels;
    labels.push_back("Locus");
    labels.push_back("minusLog10PValue");
//     if(K>1)labels.push_back("Population");
//     labels.push_back( "Score");
//     labels.push_back("CompleteInfo");
//     labels.push_back("ObservedInfo");
//     labels.push_back("PercentInfo");
//     labels.push_back("Missing1");
//     labels.push_back("Missing2");
//     labels.push_back("StdNormal");
//     labels.push_back("PValue");

    vector<int> dimensions(3,0);
    dimensions[0] = labels.size();
    dimensions[1] = L * KK;
    dimensions[2] = (int)(numPrintedIterations);
    
    R_output3DarrayDimensions(&outputfile, dimensions, labels);
  }
}

void AncestryAssocTest::Update(int locus, const double* Covariates, double phi, double YMinusEY, double DInvLink, 
			       const vector<vector<double> > Probs) {
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
  //copy covariates into X
  for( unsigned k = 0; k < K - 1; k++ ){
    X[ K + k] = Covariates[k];
    BX[k] = Xcov[k] = Covariates[k];
  }

  for( unsigned k = 0; k < K ; k++ ){
    Xcopy[k] = X[k] = Probs[1][k] + 2.0 * Probs[2][k];//Conditional expectation
    VarA[k] = Probs[1][k]*(1.0 - Probs[1][k]) + 4.0*Probs[2][k]*Probs[0][k];//conditional variances
  }
  //KLUDGE: need to reset Xcopy each time since destroyed in computation of score
  Xcopy[2*K-1] = 1;
  for( unsigned k = 0; k < K-1; k++ )Xcopy[k + K] = Covariates[k];

  try{
    // ** compute expectation of score **
    scale_matrix(Xcopy, YMinusEY*phi, 2*K, 1);      //Xcopy *= YMinusEY *phi
    add_matrix(Score[locus], Xcopy, 2*K, 1);//Score[locus] += Xcopy
    // ** compute uncorrected info **
    matrix_product(X, X, XX, 2*K, 1, 2*K);        //XX = X'X
    scale_matrix(XX, DInvLink*phi, 2*K, 2*K);     //XX = DInvLink * phi * X'X
    add_matrix(Info[locus], XX, 2*K, 2*K);//Info[locus] += XX
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
    InfoCorrection[locus][k] += VarA[k] * (DInvLink *phi - phi * phi * DInvLink * DInvLink * xBx[0]); 
    VarScore[locus][k] += VarA[k] * phi * phi * YMinusEY * YMinusEY;
  }
  delete[] X; delete[] Xcopy;
  delete[] XX;
  delete[] BX;
  delete[] VarA;
}

void AncestryAssocTest::Accumulate(unsigned j){
  double *score = new double[K], *info = new double[K*K];
  try{
    CentredGaussianConditional(K, Score[j], Info[j], score, info, 2*K );
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
    SumScore[j*KK +k] += score[k+k1];
    SumInfo[j*KK +k] += info[(k+k1)*K +k+k1] + InfoCorrection[j][k+k1];
    SumScore2[j*KK +k] += score[k+k1] * score[k+k1];
    SumVarScore[j*KK +k] += VarScore[j][k+k1];
  }
  delete[] score;
  delete[] info;

}

//TODO: fix to be called within update
void AncestryAssocTest::UpdateB(double DInvLink, double dispersion, const double* Covariates){
  //increment B using covariates
  //Xcov is a vector of covariates as in UpdateScoreForAncestry
    Xcov[K-1] = 1;//last entry is intercept
    if(K==1){//case of no covariates, Xcov = 1
      B[0] += DInvLink*dispersion;
    }
    else{
      for( unsigned k = 0; k < K - 1; k++ ){
	Xcov[k] = Covariates[k]; //centred covariates
      }
      double *temp = new double[K*K];
      matrix_product(Xcov, Xcov, temp, K, 1, K);
      scale_matrix(temp, DInvLink*dispersion, K, K);
      add_matrix(B, temp, K, K);
    delete[] temp;
    }
}
