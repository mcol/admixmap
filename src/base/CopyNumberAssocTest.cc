//=============================================================================
//
// Copyright (C) 2006  David O'Donnell, Clive Hoggart and Paul McKeigue
// Portions Copyright (C) 2010  Marco Colombo
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License version 2 or later as published by
// the Free Software Foundation.
//
// This software is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
// details.
//
// You should have received a copy of the GNU General Public License along with
// this software; see the file COPYING.  If not, it can be found at
// http://www.gnu.org/copyleft/gpl.html or by writing to the Free Software
// Foundation, 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
//
//=============================================================================

//=============================================================================
/// \file CopyNumberAssocTest.cc
/// Implementation of the CopyNumberAssocTest class.
//=============================================================================

#include "CopyNumberAssocTest.h"
#include "Genome.h"
#include "bclib/linalg.h"
#include "bclib/dist.h"
#include "bclib/misc.h"
#include "bclib/DelimitedFileWriter.h"

using namespace::std;

CopyNumberAssocTest::CopyNumberAssocTest(){
  SumScore = 0;
  SumInfo = 0;
  SumVarScore = 0;
  SumScore2 = 0;
  Score = 0;
  Info = 0;
  VarScore = 0;
  InfoCorrection = 0;
  B = 0;
  Xcov = 0;
  NumStrata = 0;
  NumOutputStrata = 0;
  L = 0;
  numPrintedIterations = 0;
  numUpdates = 0;
  LU_B = 0;
  LU_P = 0;
}

CopyNumberAssocTest::~CopyNumberAssocTest(){
  if(test){
    delete[] SumScore;
    delete[] SumInfo;
    delete[] SumVarScore;
    delete[] SumScore2;
    delete[] B;
    delete[] LU_B;
    delete[] LU_P;
    delete[] Xcov;
    bclib::free_matrix(Score, L);
    bclib::free_matrix(Info, L);  
    bclib::free_matrix(VarScore, L);
    bclib::free_matrix(InfoCorrection, L);
  }
}
void CopyNumberAssocTest::Initialise(const char* filename, const int numStrata, const int NumLoci){
  test = true;
  R.open(filename);

  NumStrata = numStrata;
  L = NumLoci;
  NumOutputStrata = NumStrata;
  if(NumStrata == 2){
    NumOutputStrata = 1;//only keeping scores for second population when 2 populations
  }

  const int dim = L * NumOutputStrata;

  SumScore    = new double[dim];
  SumInfo     = new double[dim];
  SumScore2   = new double[dim];
  SumVarScore = new double[dim];
  fill(SumScore, SumScore + dim, 0.0);
  fill(SumScore2, SumScore2 + dim, 0.0);
  fill(SumInfo, SumInfo + dim, 0.0);
  fill(SumVarScore, SumVarScore + dim, 0.0);

  using bclib::alloc2D_d;
  Score = alloc2D_d(L, 2*NumStrata);
  Info = alloc2D_d(L, 4*NumStrata*NumStrata);
  VarScore = alloc2D_d(L, NumStrata);
  InfoCorrection = alloc2D_d(L, NumStrata);
  B = new double[NumStrata * NumStrata];
  LU_B = new double[NumStrata * NumStrata];
  LU_P = new size_t[NumStrata];
  Xcov = new double[NumStrata];
  fill(B, B + NumStrata * NumStrata, 0.0);
}

void CopyNumberAssocTest::Reset(){
  if(test){
    const unsigned KK = NumStrata * NumStrata;
    const unsigned K2 = 2* NumStrata;
    const unsigned KK4 = 4 * KK;

    for(unsigned j = 0; j < L; ++j){
      for(unsigned k = 0; k < K2; ++k)
	Score[j][k] = 0.0;
      for(unsigned k = 0; k < KK4; ++k)
	Info[j][k] = 0.0;
      for(unsigned k = 0; k < NumStrata; ++k){
	InfoCorrection[j][k] = 0.0;
	VarScore[j][k] = 0.0;
      }
    }

    // compute the LU factors for the B from the previous iteration
    copy(B, B + KK, LU_B);
    bclib::LU_decomp(NumStrata, LU_B, LU_P);

    // B accumulates the sum for the current iteration
    fill(B, B + KK, 0.0);
    fill(Xcov, Xcov + NumStrata, 0.0);
  }
}
void CopyNumberAssocTest::OutputCopyNumberAssocTest(unsigned j, unsigned k,  
                                         bclib::DelimitedFileWriter& outfile,
                                         const string& label, bool final) {
  OutputRaoBlackwellizedScoreTest(outfile, label, SumScore[ j*NumOutputStrata + k], SumScore2[ j*NumOutputStrata + k], 
				  SumVarScore[ j*NumOutputStrata + k ],SumInfo[ j*NumOutputStrata + k ], final); 
}

void CopyNumberAssocTest::Update(int locus, const double* Covariates,
                                 double phi, double YMinusEY, double DInvLink,
                                 bool diploid, const vector<vector<double> >& Probs) {
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
  
  double *X = new double [2 * NumStrata], *XX = new double[4*NumStrata*NumStrata];
  double xBx[1];
  double* BX = new double[NumStrata];
  double* VarA = new double[NumStrata];
 
  X[ 2*NumStrata - 1] = 1;//intercept
  Xcov[NumStrata-1] = 1;
  //copy covariates into X
  for( unsigned k = 0; k < NumStrata - 1; k++ ){
    X[ NumStrata + k] = Covariates[k];
    BX[k] = Xcov[k] = Covariates[k];
  }

  for( unsigned k = 0; k < NumStrata ; k++ ){
    if(diploid) {
    X[k] = Probs[1][k] + 2.0 * Probs[2][k]; // Conditional expectation
    VarA[k] = Probs[1][k]*(1.0 - Probs[1][k]) + 4.0*Probs[2][k]*Probs[0][k];//conditional variances
    } else {
      // haploid - effect of one extra copy from pop k is equivalent to two extra copies in diploid case
      X[k] = 2 * Probs[1][k]; // Conditional expectation
      VarA[k] = 4 * Probs[1][k] * (1.0 - Probs[1][k]);//conditional variances
    }
  }

  try{
    using namespace bclib;
    // ** compute uncorrected info **
    matrix_product(X, X, XX, 2*NumStrata, 1, 2*NumStrata);        //XX = X'X
    scale_matrix(XX, DInvLink*phi, 2*NumStrata, 2*NumStrata);     //XX = DInvLink * phi * X'X
    add_matrix(Info[locus], XX, 2*NumStrata, 2*NumStrata);//Info[locus] += XX

    // ** compute expectation of score **
    scale_matrix(X, YMinusEY*phi, 2*NumStrata, 1); // X *= YMinusEY * phi
    add_matrix(Score[locus], X, 2*NumStrata, 1);   // Score[locus] += X

    // ** compute variance of score and correction term for info **    
    LU_solve(NumStrata, LU_B, LU_P, Xcov, BX);     // BX = inv(PrevB) * Xcov
    matrix_product(Xcov, BX, xBx, 1, NumStrata, 1);//xBx = Xcov' * BX
  }
  catch(string s){
    delete[] X;
    delete[] XX;
    delete[] BX;
    delete[] VarA;
    std::string error_string = "Error in CopyNumberAssocTest::Update:\n";
    error_string.append(s);
    throw(error_string);//throw error message to top level
  }
  for(unsigned k = 0; k < NumStrata ; k++ ){
    InfoCorrection[locus][k] += VarA[k] * (DInvLink *phi - phi * phi * DInvLink * DInvLink * xBx[0]); 
    VarScore[locus][k] += VarA[k] * phi * phi * YMinusEY * YMinusEY;
  }
  delete[] X;
  delete[] XX;
  delete[] BX;
  delete[] VarA;
}

///accumulate score, info, scoresq and var score over iterations
void CopyNumberAssocTest::Accumulate(){

  double *CentredScore = new double[NumStrata];
  double *CentredInfo  = new double[NumStrata * NumStrata];

  //increment update counter
  for(unsigned j = 0; j < L; ++j)
    Accumulate(j, CentredScore, CentredInfo);
  ++numUpdates;

  delete[] CentredScore;
  delete[] CentredInfo;
}

void CopyNumberAssocTest::Accumulate(unsigned j, double *CentredScore, double *CentredInfo) {

  //centre score and info by conditioning on covariates
  try{
    bclib::CentredGaussianConditional(NumStrata, Score[j], Info[j], CentredScore, CentredInfo, 2*NumStrata );
  }
  catch(string s){
    string error = "Error centring ancestry association scores:\n";
    error.append(s);
    throw(error);
  }
  
  //accumulate over iterations
  //for two populations, accumulate the scores for second population only
  unsigned k1 = 0;
  if(NumStrata == 2)
    k1 = 1;
     
  for( unsigned k = 0; k < NumOutputStrata ; k++ ){
    SumScore[j*NumOutputStrata +k] += CentredScore[k+k1];
    SumInfo[j*NumOutputStrata +k] += CentredInfo[(k+k1)*NumStrata +k+k1] + InfoCorrection[j][k+k1];
    SumScore2[j*NumOutputStrata +k] += CentredScore[k+k1] * CentredScore[k+k1];
    SumVarScore[j*NumOutputStrata +k] += VarScore[j][k+k1];
  }
}

//TODO: fix to be called within update
void CopyNumberAssocTest::UpdateB(double DInvLink, double dispersion, const double* Covariates){
  //increment B using covariates
  //Xcov is a vector of covariates as in UpdateScoreForAncestry
    Xcov[NumStrata-1] = 1;//last entry is intercept
    if(NumStrata==1){//case of no covariates, Xcov = 1
      B[0] += DInvLink*dispersion;
    }
    else{
      for( unsigned k = 0; k < NumStrata - 1; k++ ){
	Xcov[k] = Covariates[k]; //centred covariates
      }
      double *temp = new double[NumStrata*NumStrata];
      bclib::matrix_product(Xcov, Xcov, temp, NumStrata, 1, NumStrata);
      bclib::scale_matrix(temp, DInvLink*dispersion, NumStrata, NumStrata);
      bclib::add_matrix(B, temp, NumStrata, NumStrata);
    delete[] temp;
    }
}
