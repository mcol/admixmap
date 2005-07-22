/** 
 *   ADMIXMAP
 *   Individual.cc 
 *   Class to represent an individual and update individual-level parameters
 *   Copyright (c) 2002, 2003, 2004, 2005 LSHTM
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
#include "Individual.h"
#include "StringConvertor.h"
#include <algorithm>

#define PR(x) cout << #x << " = " << x << endl;

double *Individual::AffectedsScore = 0;
double *Individual::AffectedsVarScore = 0;
double *Individual::AffectedsInfo = 0;
double **Individual::AncestryScore = 0;
double **Individual::AncestryInfo = 0;
double **Individual::AncestryVarScore = 0;
double **Individual::AncestryInfoCorrection = 0;
double *Individual::B;
double *Individual::PrevB;
double *Individual::Xcov;

unsigned int Individual::numChromosomes;
Genome *Individual::Loci;
int *Individual::sumxi;
double Individual::Sumrho0;
int Individual::Populations;

//temporary function for marg likelihood calculation
//given an array representing a matrix, returns a vector, containing one column
static Vector_d GetRow(double *x, int row, int ncols){
  Vector_d v(ncols);
  for(int i = 0; i < ncols; ++i)v(i) = x[row*ncols + i];
  return v;
}
// static Vector_d GetRow(int *x, int row, int ncols){
//   Vector_d v(ncols);
//   for(int i = 0; i < ncols; ++i)v(i) = x[row*ncols + i];
//   return v;
// }

Individual::Individual()
{//should initialise pointers here
}

Individual::Individual(int mynumber,AdmixOptions* options, InputData *Data, Genome& Loci,Chromosome **chrm)
{
  if( options->getRhoIndicator() ){
        TruncationPt = options->getTruncPt();
        if( options->isRandomMatingModel() )
	  if(!options->RhoFlatPrior() && !options->logRhoFlatPrior() )
                _rho.assign(2,options->getRhoalpha());
            else
                _rho.assign(2,1);
        else
	  if(!options->RhoFlatPrior() && !options->logRhoFlatPrior() )
                _rho.assign(1,options->getRhoalpha());
            else
                _rho.assign(2,1);
    }

    // Read sex value if present.
    sex = 1;
    if (options->getgenotypesSexColumn() == 1) {
	sex = Data->GetSexValue(mynumber);
    }

    int numCompositeLoci = Loci.GetNumberOfCompositeLoci();

    sumxi = new int[numCompositeLoci];

    LocusAncestry = new int*[ numChromosomes ]; // array of matrices in which each col stores 2 integers 

    Theta = 0;
    ThetaX = 0;
    ThetaProposal = 0;
    ThetaXProposal = 0;
 
    // SumLocusAncestry is sum of locus ancestry states over loci at which jump indicator xi is 1  
    SumLocusAncestry = new int[options->getPopulations()*2];

    if( options->isRandomMatingModel() ){//random mating model
      ThetaProposal = new double[ Populations * 2 ];
    }
    else{
      ThetaProposal = new double[ Populations];
    }

    //X chromosome objects
    SumLocusAncestry_X = 0;    
//    if(Loci.isX_data() ){
      if( sex == 1 ){
	ThetaXProposal = new double[ Populations];
	SumLocusAncestry_X = new int[Populations];
      }
      else{
	ThetaXProposal = new double[ Populations * 2 ];
	SumLocusAncestry_X = new int[Populations * 2 ];
      }
//    }

    // vector of possible haplotype pairs - expect 2 integers per locus 
                                                        // or 1 integer (haploid)
    PossibleHapPairs = new vector<hapPair >[numCompositeLoci];
 
    X_posn = 9999;
    string s1("\"X\""), s2("X");
    size_t AncestrySize = 0;
  // set size of locus ancestry array
    //gametes holds the number of gametes for each chromosome, either 1 or 2
    //X_posn is the number of the X chromosome. Note that this is not necessarily 22, even in humans
    for( unsigned int j = 0; j < numChromosomes; j++ ){
        if( chrm[j]->GetLabel(0) != s1  && chrm[j]->GetLabel(0) != s2){// if not X chromosome, set number of elements to 2, num loci
	  AncestrySize = 2 * chrm[j]->GetSize() ;
	  gametes.push_back(2);
        }
        else if( sex != 2 ){//male or missing
	  AncestrySize = chrm[j]->GetSize() ;
	  gametes.push_back(1);
	  X_posn = j;
        }
        else{//female
	  AncestrySize = 2 * chrm[j]->GetSize() ;
	  gametes.push_back(2);
	  X_posn = j;
        }
	LocusAncestry[j] = new int[ AncestrySize];
        for(unsigned i = 0; i < AncestrySize; ++i)LocusAncestry[j][i] = 0;

    }
    //retrieve genotypes
    Data->GetGenotype(mynumber, options->getgenotypesSexColumn(), Loci, &genotypes);

    // loop over composite loci to set possible haplotype pairs compatible with genotype 
    for(int j=0;j<numCompositeLoci;++j) {
       Loci(j)->setPossibleHaplotypePairs(genotypes[j], PossibleHapPairs[j]);
    }

}

void Individual::SetStaticMembers(Genome *pLoci, AdmixOptions *options){
  Loci = pLoci;
  numChromosomes = Loci->GetNumberOfChromosomes();
  Populations = options->getPopulations();
  int L = Loci->GetNumberOfCompositeLoci();

  AncestryScore = 0;
  AncestryInfo = 0;
  AncestryVarScore = 0;
  AncestryInfoCorrection = 0;
  AffectedsScore = 0;
  AffectedsVarScore = 0;
  AffectedsInfo = 0;
  B = 0;
  PrevB = 0;
  Xcov = 0;

  int K = Populations;

  if( options->getTestForAffectedsOnly() ){
    if(Populations == 2) K = 1;
    AffectedsScore = new double[L * K];
    AffectedsVarScore = new double[L * K];
    AffectedsInfo = new double[L * K];
  }
  if( options->getTestForLinkageWithAncestry() ){
    AncestryScore = alloc2D_d(L, 2*Populations);
    AncestryInfo = alloc2D_d(L, 4*Populations*Populations);
    AncestryVarScore = alloc2D_d(L, Populations);
    AncestryInfoCorrection = alloc2D_d(L, Populations);
    B = new double[Populations * Populations];
    PrevB = new double[Populations * Populations];
    Xcov = new double[Populations];
  }
}

Individual::~Individual()
{
  delete[] PossibleHapPairs;
  delete[] LocusAncestry;
  delete[] SumLocusAncestry;
  delete[] SumLocusAncestry_X;  
  delete[] Theta;
  delete[] ThetaX;
  delete[] ThetaProposal;
  delete[] ThetaXProposal;
}
void Individual::DeleteStaticMembers(){
  delete[] sumxi;
  delete[] AffectedsScore;
  delete[] AffectedsInfo;
  delete[] AffectedsVarScore;
  free_matrix(AncestryScore, Loci->GetNumberOfCompositeLoci());
  free_matrix(AncestryInfo, Loci->GetNumberOfCompositeLoci());  
  free_matrix(AncestryVarScore, Loci->GetNumberOfCompositeLoci());
  free_matrix(AncestryInfoCorrection, Loci->GetNumberOfCompositeLoci());
}

void Individual::ResetStaticSums(){
  for(unsigned int i = 0; i < Loci->GetNumberOfCompositeLoci(); ++i)sumxi[i] = 0;
  Sumrho0 = 0;
}

unsigned short **Individual::getGenotype(unsigned int locus){
  return genotypes[locus];
}

std::vector<hapPair > &Individual::getPossibleHapPairs(unsigned int locus){
  return PossibleHapPairs[locus];
}

//Indicates whether a genotype is missing at a locus
bool Individual::IsMissing(unsigned int locus)
{
  unsigned int count = 0;
  int NumberOfLoci = (*Loci)(locus)->GetNumberOfLoci();
  for(int i = 0; i < NumberOfLoci; i++){
    count += genotypes[locus][i][0];
  }
  return (count == 0);
}

double *Individual::getAdmixtureProps()
{
  return Theta;
}

//should call this initialise not set
void Individual::setAdmixtureProps(double *a, size_t size)
{
  delete[] Theta; //safe if Theta initialised to 0
  Theta = new double[size];
  for(unsigned i = 0; i < size; ++i)  Theta[i] = a[i];
}

// Matrix_d& Individual::getAdmixturePropsX()
// {
//   return ThetaX;
// }

void Individual::setAdmixturePropsX(double *a, size_t size)
{
  delete[] ThetaX;
  ThetaX = new double[size];
  for(unsigned i = 0; i < size; ++i)  ThetaX[i] = a[i];
}

int Individual::getSex()
{
   return sex;
}

int *Individual::getSumXi()
{
   return sumxi;
}

int Individual::getSumXi(int j){
  return sumxi[j];
}

double Individual::getSumrho0()
{
   return Sumrho0;
}

double Individual::getSumrho()
{
   double sumrho = 0;
   for( unsigned int g = 0; g < _rho.size(); g++ )
      sumrho += _rho[g];
   return sumrho;
}

vector<double> Individual::getRho()
{
   return _rho;
}

Vector_i Individual::GetLocusAncestry(int chrm, int locus){
  Vector_i Ancestry(2);
  Ancestry(0)  = LocusAncestry[chrm][locus];
  if((unsigned)chrm == X_posn)Ancestry(1) = Ancestry(0);
  else Ancestry(1) = LocusAncestry[chrm][Loci->GetSizesOfChromosomes()[chrm]  + locus];
  return Ancestry; 
}

int Individual::GetLocusAncestry(int chrm, int gamete, int locus){
  int g = (gametes[chrm] == 2) ? gamete : 0; //so that gamete = 1 works when gametes[chrm] = 1;
  return LocusAncestry[chrm][g * Loci->GetSizesOfChromosomes()[chrm]  + locus] ;
}

double Individual::getLogPosteriorProb()
{
   return LogPosterior;
}

// update the individual admixture values (mean of both gametes) used in the regression model
void Individual::UpdateAdmixtureForRegression( int i,int Populations, int NoCovariates,
                                               const double *poptheta, bool RandomMatingModel, Matrix_d *Covariates0)
{
  double avgtheta[Populations];
  if(RandomMatingModel )//average over gametes
    for(int k = 0; k < Populations; ++k) avgtheta[k] = (Theta[k] + Theta[k + Populations]) / 2.0;    
  
  else
    for(int k = 0; k < Populations; ++k) avgtheta[k] = Theta[k];
  for( int k = 0; k < Populations - 1; k++ )
    (*Covariates0)( i, NoCovariates - Populations + k + 1 )
      = avgtheta[ k + 1 ] - poptheta[ k + 1 ];
}

// Metropolis update for admixture proportions theta, taking log of acceptance probability ratio as argument
// uses log ratio because this is easier than ratio to calculate for linear regression model
// if no regression model, logpratio remains set to 0, so all proposals are accepted 
void Individual::Accept_Reject_Theta( double logpratio, bool xdata, int Populations, bool RandomMatingModel )
{
  bool test = true;
  // loop over populations: if element of Dirichlet parameter vector is 0, do not update corresponding element of 
  // admixture proportion vector
  for( int k = 0; k < Populations; k++ ){
    if( ThetaProposal[ k ] == 0.0 )
      test = false;
    else if( RandomMatingModel && ThetaProposal[ k + Populations ] == 0.0 )
      test = false;
  }

  size_t size_admix;
  if(RandomMatingModel) size_admix = Populations *2;
  else size_admix = Populations;
  // generic Metropolis rejection step
  if( logpratio < 0 ){
     if( test && log(myrand()) < logpratio ){
       setAdmixtureProps(ThetaProposal, size_admix);
        if( xdata )
	  setAdmixturePropsX(ThetaXProposal, size_admix);
     }
  }
  else{
    setAdmixtureProps(ThetaProposal, size_admix);
     if( xdata )
       setAdmixturePropsX(ThetaXProposal, size_admix);
  }
}

// these two functions return log of ratio of likelihoods of new and old values of population admixture
// in regression models.  individual admixture theta is standardized about the mean poptheta calculated during burn-in. 
 
// should have just one function to get the likelihood in the regression model, given a value of population admixture
// should be generic method for GLM, given Xbeta, Y and probability distribution 
// then should calculate ratio in Metropolis step 
//  
double Individual::AcceptanceProbForTheta_LogReg( int i, int TI, bool RandomMatingModel,int Populations,
					      int NoCovariates, Matrix_d &Covariates0, double **beta, double **ExpectedY,
						  Matrix_d *Outcome, const double *poptheta)
{
  double probratio, Xbeta = 0;
  // TI = Target Indicator, indicates which outcome var is used
  double avgtheta[Populations];
  // calculate mean of parental admixture proportions
  if( RandomMatingModel )
    for(int k = 0;k < Populations; ++k)avgtheta[k] = (ThetaProposal[k] + ThetaProposal[k + Populations])/ 2.0 - poptheta[k];
  else
    for(int k = 0;k < Populations; ++k)avgtheta[k] = ThetaProposal[k]  - poptheta[k];

  for( int jj = 0; jj < NoCovariates - Populations + 1; jj++ )
    Xbeta += Covariates0( i, jj ) * beta[TI][jj];
  for( int k = 0; k < Populations - 1; k++ ){
    //? Old code had 0 instead of TI in for index of beta.
    Xbeta += avgtheta[ k ] * beta[TI][NoCovariates - Populations + k + 1];}
  double newExpectedY = 1.0 / ( 1.0 + exp( -Xbeta ) );
  if( Outcome[ TI ]( i, 0 ) == 1 )
    probratio = newExpectedY / ExpectedY[ TI ][i];
  else
    probratio = ( 1 - newExpectedY ) / ( 1 - ExpectedY[ TI ][i] );

  return( log(probratio) );
} 

double Individual::AcceptanceProbForTheta_LinearReg( int i, int TI,  bool RandomMatingModel, int Populations,
						 int NoCovariates, Matrix_d &Covariates0, double **beta, double **ExpectedY,
						 Matrix_d *Outcome, const double *poptheta, double *lambda)
{
  double prob, Xbeta = 0;
   double avgtheta[Populations];
  if( RandomMatingModel )
    for(int k = 0;k < Populations; ++k)avgtheta[k] = (ThetaProposal[k] + ThetaProposal[k + Populations ])/ 2.0 - poptheta[k];
  else
    for(int k = 0;k < Populations; ++k)avgtheta[k] = ThetaProposal[k]  - poptheta[k];

  for( int jj = 0; jj < NoCovariates - Populations + 1; jj++ )
    Xbeta += Covariates0( i, jj ) * beta[ TI ][jj];
  for( int k = 0; k < Populations - 1; k++ ){
    Xbeta += avgtheta[ k ] * beta[ TI ][NoCovariates - Populations + k + 1];
  }

  prob = 0.5 * lambda[ TI ] * (( ExpectedY[ TI ][i] - Outcome[ TI ]( i, 0 ) ) * ( ExpectedY[ TI ][i] - Outcome[ TI ]( i, 0 ) )
			       - ( Xbeta - Outcome[ TI ]( i, 0 ) ) * ( Xbeta - Outcome[ TI ]( i, 0 ) ) );

  return( prob );
}

double Individual::AcceptanceProbForTheta_XChrm(std::vector<double> &sigma, int Populations )
{
   int gametes = 1;
   if( sex == 2 )
      gametes = 2;
   double p = 0, sum1 = 0.0, sum2 = 0.0;
   //Matrix_d ThetaOld = Theta, ThetaXOld = ThetaX;
   for( int g = 0; g < gametes; g++ ){
     sum1 = sum2 = 0.0;
     for( int k = 0; k < Populations; k++ ) {
       sum1 += ThetaProposal[g*Populations + k];
       sum2 += Theta[g*Populations + k];
     }
       p += gsl_sf_lngamma( sigma[g]*sum1 )
         - gsl_sf_lngamma( sigma[g]*sum2 );
      for( int k = 0; k < Populations; k++ ){
         p += gsl_sf_lngamma( sigma[g]*Theta[g*Populations + k] ) - gsl_sf_lngamma( sigma[g]*ThetaProposal[g*Populations +k] );
         p += (sigma[g]*ThetaProposal[g*Populations + k]-1.0)*log(ThetaXProposal[g*Populations +k]) - 
	   (sigma[g]*Theta[g*Populations + k]-1.0) *log(ThetaX[g*Populations + k]);
      }
   }
   return p;
}

void Individual::SampleParameters( int i, double *SumLogTheta, AlleleFreqs *A, int iteration , Matrix_d *Outcome,
				  int NumOutcomes,  int* OutcomeType, double **ExpectedY, double *lambda, int NoCovariates,
				   Matrix_d &Covariates0, double **beta, const double *poptheta,
				   AdmixOptions* options, Chromosome **chrm, 
				   vector<vector<double> > &alpha, double rhoalpha, 
				   double rhobeta,vector<double> sigma, double DInvLink, double dispersion)
//Outcome = Outcome variable(s)
//ExpectedY = expected outcome variable
//lambda = precision in linear regression model (if there is one)
//NoCovariates = # covariates
//alpha = pop admixture Dirichlet parameters
//rhoalpha, rhobeta = shape and scale parameters in prior for rho
//sigma = 
//DInvLink = Derivative Inverse Link function in regression model, used in ancestry score test
//dispersion = dispersion parameter in regression model (if there is one) = lambda for linear reg, 1 for logistic
{
  // ** reset SumLocusAncestry and ThetaProposal **
  // should clean this up, takes 4 lines with std::vectors
  for(int j = 0; j < Populations *2; ++j)SumLocusAncestry[j] = 0;
//  if(Loci->isX_data() ){
    int J = Populations;
    if(sex != 1) J *=2;
    for(int j = 0; j < J ;++j)SumLocusAncestry_X[j] = 0;
//  }

  size_t size_theta;
   if( options->isRandomMatingModel() )
     size_theta = Populations*2; // double the size for 2 gametes in RMM
   else//assortative mating
     size_theta = Populations;

   for(unsigned k = 0; k < size_theta; ++k){
     ThetaProposal[k] = 0.0;
   }
   if(ThetaXProposal){
     size_theta = Populations;if(sex != 1) size_theta *= 2;
     for(unsigned k = 0; k < size_theta; ++k) ThetaXProposal[k] = 0.0;
   }

  //SumN is the number of arrivals between each pair of adjacent loci  
  unsigned int SumN[] = {0,0};
  unsigned int SumN_X[] = {0,0};
    
  bool isdiploid;
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    //Update Forward/Backward probs in HMM
    isdiploid = UpdateForBackProbs(j, chrm[j], options);

    //update score tests for linkage with ancestry for *previous* iteration
    if(iteration > options->getBurnIn()){
      //Update affecteds only scores
      if( options->getAnalysisTypeIndicator() == 0 ) 
	UpdateScoreForLinkageAffectedsOnly(j, options->isRandomMatingModel(), 
					   chrm );
      else if( options->getTestForAffectedsOnly() && Outcome[0](i,0) == 1 ){
	UpdateScoreForLinkageAffectedsOnly(j, options->isRandomMatingModel(), 
					   chrm );
      }
      //update ancestry score tests
      if( options->getTestForLinkageWithAncestry() ){   
	UpdateScoreForAncestry(j,dispersion, Outcome[0](i,0) - ExpectedY[0][i], DInvLink,chrm);
      }    
    }

    //Sample locus ancestry
    // sampling locus ancestry requires calculation of forward probability vectors alpha in HMM 
    chrm[j]->SampleLocusAncestry(LocusAncestry[j], Theta, isdiploid);
 
   //loop over loci on current chromosome and update allele counts
    for( unsigned int jj = 0; jj < chrm[j]->GetSize(); jj++ ){
      int locus =  chrm[j]->GetLocus(jj);
      if( !(IsMissing(j)) ){
	  int h[2];//to store sampled hap pair
	  //might be a shortcut for haploid data since there is only one compatible hap pair, no need to sample
	  (*Loci)(locus)->SampleHapPair(h, PossibleHapPairs[locus], GetLocusAncestry(j,jj));
	  A->UpdateAlleleCounts(locus,h, GetLocusAncestry(j,jj), isdiploid);
      }
     }   

    //sample number of arrivals and update sumxi, sumrho0 and SumLocusAncestry
    bool isX = (j == X_posn);
    chrm[j]->SampleJumpIndicators(LocusAncestry[j], gametes[j], sumxi, &Sumrho0, SumLocusAncestry, SumLocusAncestry_X, isX,
				  SumN, SumN_X, options->getRhoIndicator());
  }//end chromosome loop

  double L = Loci->GetLengthOfGenome(), L_X=0.0;
  if( Loci->isX_data() ) L_X = Loci->GetLengthOfXchrm();

  // sample sum of intensities parameter rho
  if( options->getRhoIndicator() ){
     SampleRho( options->getXOnlyAnalysis(), options->isRandomMatingModel(), Loci->isX_data(), rhoalpha, rhobeta, L, L_X, 
	       SumN, SumN_X);
  }

  if( options->getAnalysisTypeIndicator() > 1 ){
    // sample missing values of outcome variable
    double u;
    for( int k = 0; k < NumOutcomes; k++ ){
      if( Outcome[k].IsMissingValue( i, 0 ) ){
	if( !OutcomeType[k] ) // linear regression
	  Outcome[k]( i, 0 ) = gennor( ExpectedY[k][i], 1 / sqrt( lambda[k] ) );
	else{// logistic regression
	  u = myrand();
	  if( u * ExpectedY[k][i] < 1 )
	    Outcome[k]( i, 0 ) = 1;
	  else
	    Outcome[k]( i, 0 ) = 0;
	}
      }
    }
  }

  //sample admixture proportions, Theta
  SampleTheta(i, SumLogTheta,Outcome, NumOutcomes, OutcomeType, ExpectedY, lambda, NoCovariates,
	      Covariates0, beta, poptheta, options, alpha, sigma);

  //increment B using new Admixture Props
  //Xcov is a vector of admixture props as covariates as in UpdateScoreForAncestry
  if(iteration >= options->getBurnIn() && options->getTestForLinkageWithAncestry()){
    Xcov[Populations-1] = 1;//last entry is intercept
    for( int k = 0; k < Populations - 1; k++ ){
      Xcov[k] = Theta[ k]; 
    }
    double *temp = new double[Populations*Populations];
    matrix_product(Xcov, Xcov, temp, Populations, 1, Populations);
    scale_matrix(temp, DInvLink*dispersion, Populations, Populations);
    add_matrix(B, temp, Populations, Populations);
    delete[] temp;
  }

  //calculate log posterior if necessary 
  if( options->getMLIndicator() && i == 0 && iteration > options->getBurnIn() )
    CalculateLogPosterior(options,Loci->isX_data(), alpha, rhoalpha, rhobeta, L, L_X, SumN, SumN_X);
  

}

// samples individual admixture proportions
void Individual::SampleTheta( int i, double *SumLogTheta, Matrix_d *Outcome,
                                 int NumOutcomes,  int* OutcomeType, double **ExpectedY, double *lambda, int NoCovariates,
                                  Matrix_d &Covariates0, double **beta, const double *poptheta,
                                  AdmixOptions* options, vector<vector<double> > &alpha, vector<double> sigma)
{
  // propose new value for individual admixture proportions
  // should be modified to allow a population mixture component model
  ProposeTheta(options, sigma, alpha);       
  int K = Populations;
  
  double logpratio = 0;
  //calculate Metropolis acceptance probability ratio for proposal theta    
  //should have one function to do this

  //linear regression case
  if( options->getAnalysisTypeIndicator() == 2 && !options->getTestForAdmixtureAssociation() ){
    logpratio = AcceptanceProbForTheta_LinearReg( i, 0, options->isRandomMatingModel(), K,
                                                  NoCovariates, Covariates0, beta, ExpectedY, Outcome, poptheta,lambda);
  }
  //logistic regression case
  else if( (options->getAnalysisTypeIndicator() == 3 || options->getAnalysisTypeIndicator() == 4) && !options->getTestForAdmixtureAssociation() ){
    logpratio = AcceptanceProbForTheta_LogReg( i, 0, options->isRandomMatingModel(), K,
				       NoCovariates, Covariates0, beta, ExpectedY, Outcome, poptheta);
    }
  //case of both linear and logistic regressions
  else if( options->getAnalysisTypeIndicator() == 5 ){
    for( int k = 0; k < NumOutcomes; k++ ){
      if( OutcomeType[ k ] )
	logpratio += AcceptanceProbForTheta_LogReg( i, k, options->isRandomMatingModel(), K,
                                                    NoCovariates, Covariates0, beta, ExpectedY, Outcome, poptheta);
      else
	logpratio += AcceptanceProbForTheta_LinearReg( i, k, options->isRandomMatingModel(), K,
                                                       NoCovariates, Covariates0, beta, ExpectedY, Outcome, poptheta,lambda);
      }
  }
  //case of X only data
  if( Loci->isX_data() && !options->getXOnlyAnalysis() )
    logpratio += AcceptanceProbForTheta_XChrm( sigma, K);

 //Accept or reject proposed value - if no regression model, proposal will be accepted because logpratio = 0    
  Accept_Reject_Theta(logpratio, Loci->isX_data(), K, options->isRandomMatingModel() );

 // update the value of admixture proportions used in the regression model  
  if( options->getAnalysisTypeIndicator() > 1 )
    UpdateAdmixtureForRegression(i, K, NoCovariates, poptheta, options->isRandomMatingModel(),&(Covariates0));

  for( int k = 0; k < K; k++ ){
    SumLogTheta[ k ] += log( Theta[ k ] );
      if(options->isRandomMatingModel() && !options->getXOnlyAnalysis() )
	SumLogTheta[ k ] += log( Theta[ K + k ] );
    }
}

// Samples individual admixture proportions conditional on sampled values of ancestry at loci where 
// jump indicator xi is 1, population admixture distribution parameters alpha, and likelihood from regression 
// model (if there is one) 
// Proposes new value for individual admixture proportions 
// as conjugate Dirichlet posterior conditional on prior parameter vector alpha and 
// multinomial likelihood given by sampled values of ancestry at loci where jump indicator xi is 1
// proposes new values for both gametes if random mating model 
//
// should have an alternative function to sample population mixture component membership and individual admixture proportions
// conditional on genotype, not sampled locus ancestry
void Individual::ProposeTheta(AdmixOptions *options, vector<double> sigma, vector<vector<double> > &alpha){
  size_t K = Populations;
  double temp[K];//used to hold dirichlet parameters of theta posterior
  // if no regression model, sample admixture proportions theta as a conjugate Dirichlet posterior   
  if( options->getXOnlyAnalysis() ){
    for(size_t k = 0; k < K; ++k)
       temp[k] = alpha[0][k] + SumLocusAncestry_X[k];
    gendirichlet(K, temp, ThetaProposal );
  }
  else if( options->isRandomMatingModel() ){//random mating model
    for( unsigned int g = 0; g < 2; g++ ){
      if( options->getAnalysisTypeIndicator() > -1 ){
         for(size_t k = 0; k < K; ++k){
            temp[k] = alpha[0][k] + SumLocusAncestry[k + K*g];
            if( g == 0 )
               temp[k] += SumLocusAncestry_X[k];
//            cout << SumLocusAncestry_X[k] << " ";
         }
//         cout << endl;
      }
      else{//single individual
         for(size_t k = 0; k < K; ++k){
            temp[k] = alpha[g][k] + SumLocusAncestry[k + K*g];
            if( g == 0 )
               temp[k] += SumLocusAncestry_X[k];
         }
      }
      gendirichlet(K, temp, ThetaProposal+g*K );
    }
    if( Loci->isX_data() ){
      for( unsigned int g = 0; g < gametes[X_posn]; g++ ){
	for(size_t k = 0; k < K; ++k)
	  temp[k] = SumLocusAncestry_X[g*K + k] + ThetaProposal[g*K + k]*sigma[g];
	gendirichlet(K, temp, ThetaXProposal + g*K );
	}
    }
  }
  else{//assortative mating model
    for(size_t k = 0; k < K; ++k)temp[k] = alpha[0][k] + SumLocusAncestry[k] + SumLocusAncestry[k + K];
    gendirichlet(K, temp, ThetaProposal );
  }
}

//Updates forward and backward probabilities in HMM for chromosome j 
bool Individual::UpdateForBackProbs(unsigned int j, Chromosome *chrm, AdmixOptions *options){
  bool isdiploid;
  bool calcbackprobs = (options->getTestForAffectedsOnly() || options->getTestForLinkageWithAncestry());
  //Update Forward/Backward probs in HMM
  if( j != X_posn ){
    chrm->UpdateParameters(this, Theta, options, _rho, false, true, calcbackprobs);
    isdiploid = true;
  }
  else if( options->getXOnlyAnalysis() ){
    chrm->UpdateParameters( this, Theta, options, _rho, false, false, calcbackprobs);
    isdiploid = false;
  }
  else if( sex == 1 ){
    chrm->UpdateParameters( this, ThetaX, options, _rho, false, false, calcbackprobs);
    isdiploid = false;
  }
  else{
    chrm->UpdateParameters( this, ThetaX, options, _rho, false, true, calcbackprobs);
    isdiploid = true;
  }
  return isdiploid;
}

void Individual::SampleRho(bool XOnly, bool RandomMatingModel, bool X_data, double rhoalpha, double rhobeta, double L, double L_X, 
			   unsigned int SumN[], unsigned int SumN_X[]){
  // Samples sum of intensities parameter as conjugate gamma with Poisson likelihood
  // SumN is the number of arrivals between each pair of adjacent loci
  if(XOnly ){
    do{
      _rho[0] = gengam( rhobeta + L_X, rhoalpha + (double)SumN_X[0] );
    }while( _rho[0] > TruncationPt || _rho[0] < 1.0 );
  }
  else if(RandomMatingModel ){
    for( unsigned int g = 0; g < 2; g++ ){
      do{
	_rho[g] = gengam( rhobeta + L, rhoalpha + (double)SumN[g] );
      }while( _rho[g] > TruncationPt || _rho[g] < 1.0 );
    }
    if(X_data  ){
      for( unsigned int g = 0; g < gametes[X_posn]; g++ ){
	do{
	  _rho_X[g] = gengam( rhobeta + L_X, rhoalpha + (double)SumN_X[g] );
	}while( _rho[g] > TruncationPt || _rho[g] < 1.0 );
      }
    }
  }
  else{
    _rho[0] = gengam( rhobeta + 2*L, rhoalpha + (double)(SumN[0] + SumN[1]) );
  }
}

void Individual::ResetScores(AdmixOptions *options){
  int KK = Populations;
  if(Populations == 2)KK = 1;

  if( options->getTestForAffectedsOnly() )
    for(unsigned j = 0; j < Loci->GetNumberOfCompositeLoci()*KK; ++j){
      AffectedsScore[j] = 0.0;
      AffectedsVarScore[j] = 0.0;
      AffectedsInfo[j] = 0.0;
    }
  if( options->getTestForLinkageWithAncestry() ){
    for(unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); ++j){
      for(int k = 0; k < 2*Populations; ++k)
	AncestryScore[j][k] = 0.0;
      for(int k = 0; k < 4*Populations*Populations; ++k)
	AncestryInfo[j][k] = 0.0;
      for(int k = 0; k < Populations; ++k){
	AncestryInfoCorrection[j][k] = 0.0;
	AncestryVarScore[j][k] = 0.0;
      }
    }
    for(int k = 0; k < Populations*Populations; ++k){
      PrevB[k] = B[k];           //PrevB stores the sum for the previous iteration
      B[k] = 0.0;                //while B accumulates the sum for the current iteration 
    }
    for(int k = 0; k < Populations; ++k){
      Xcov[k] = 0.0;
    }
  }
}

void Individual::UpdateScoreForLinkageAffectedsOnly(int j, bool RandomMatingModel, Chromosome **chrm){
  // Different from the notation in McKeigue et  al. (2000). McKeigue
  // uses P0, P1, P2, which relate to individual admixture as follows;
  // P0 = ( 1 - theta0 ) * ( 1 - theta1 )
  // P1 = ( 1 - theta0 ) * theta1 + theta0 * ( 1 - theta1 )
  // P2 = theta0 * theta1

  // Test in McKeigue et al. (2000) was on r defined on the positive
  // real line.  This test is on log r, which should have better
  // asymptotic properties.

  //we don't bother computing scores for the first populstion when there are two
  int KK = Populations,k1 = 0;
  if(Populations ==2) {KK = 1;k1 = 1;}
  
  double theta[2];//paternal and maternal admixture proportions
  double AProbs[Populations][3];

  for( int k = 0; k < KK; k++ ){
    theta[0] = Theta[ k+k1 ];
    if( RandomMatingModel )
      theta[1] = Theta[ Populations + k+k1 ];
    else
      theta[1] = theta[0];

    int locus;
    for( unsigned int jj = 0; jj < chrm[j]->GetSize(); jj++ ){
      locus = chrm[j]->GetLocus(jj); 
      //retrieve AncestryProbs from HMM
      chrm[j]->getAncestryProbs( jj, AProbs );
      //accumulate score, score variance, and info
      AffectedsScore[locus *KK + k]+= 0.5*( AProbs[k+k1][1] + 2.0*AProbs[k+k1][2] - theta[0] - theta[1] );
      AffectedsVarScore[locus * KK + k]+= 0.25 *( AProbs[k+k1][1]*(1.0 - AProbs[k+k1][1]) + 4.0*AProbs[k+k1][2]*AProbs[k+k1][0]); 
      AffectedsInfo[locus * KK +k]+= 0.25* ( theta[0]*( 1.0 - theta[0] ) + theta[1]*( 1.0 - theta[1] ) );
      ++locus;
    }
  }
}

void Individual::UpdateScoreForAncestry(int j,double phi, double YMinusEY, double DInvLink, Chromosome **chrm)
{
  //Updates score stats for test for association with locus ancestry
  //now use Rao-Blackwellized estimator by replacing realized ancestries with their expectations
  //Notes: 1/phi is dispersion parameter
  //       = lambda[0] for linear regression, = 1 for logistic
  //       YMinusEY = Y - E(Y) = Y - g^{-1}(\eta_i)
  //       VarX = Var(X)
  //       DInvLink = {d  g^{-1}(\eta)} / d\eta = derivative of inverse-link function
  //Xcov is a vector of covariates
  //Note that only the intercept and admixture proportions are used.
  // X is (A, cov)'  
  
  double AProbs[Populations][3];
  double X[2 * Populations], Xcopy[2*Populations], XX[4*Populations*Populations];
  double xBx[1], BX[Populations];
  double VarA[Populations];
 
  X[ 2*Populations - 1] = 1;//intercept
  Xcov[Populations-1] = 1;
  //set covariates 
  for( int k = 0; k < Populations - 1; k++ ){
    X[ k + Populations] = Theta[ k ];
    BX[k] = Xcov[k] = Theta[ k ];
  }

  int locus; 
  for( unsigned int jj = 0; jj < chrm[j]->GetSize(); jj++ ){
    locus = chrm[j]->GetLocus(jj);      
    chrm[j]->getAncestryProbs( jj, AProbs );//conditional locus ancestry probs      
    
    for( int k = 0; k < Populations ; k++ ){
      Xcopy[k] = X[k] = AProbs[k][1] + 2.0 * AProbs[k][2];//Conditional expectation of ancestry
      Xcopy[k + Populations] = Theta[ k ];
      VarA[k] = AProbs[k][1]*(1.0 - AProbs[k][1]) + 4.0*AProbs[k][2]*AProbs[k][0];//conditional variances
      }
    Xcopy[2*Populations-1] = 1;
    // ** compute expectation of score **
    scale_matrix(Xcopy, YMinusEY*phi, 2*Populations, 1);      //Xcopy *= YMinusEY *phi
    add_matrix(AncestryScore[locus], Xcopy, 2*Populations, 1);//AncestryScore[locus] += Xcopy
 
    // ** compute uncorrected info **
    matrix_product(X, X, XX, 2*Populations, 1, 2*Populations);        //XX = X'X
    scale_matrix(XX, DInvLink*phi, 2*Populations, 2*Populations);     //XX = DInvLink * phi * X'X
    add_matrix(AncestryInfo[locus], XX, 2*Populations, 2*Populations);//AncestryInfo[locus] += XX

    // ** compute variance of score and correction term for info **    
    HH_solve(Populations, PrevB, Xcov, BX);          //BX = inv(PrevB) * Xcov
    matrix_product(Xcov, BX, xBx, 1, Populations, 1);//xBx = Xcov' * BX
    for( int k = 0; k < Populations ; k++ ){
      AncestryInfoCorrection[locus][k] += VarA[k] * (DInvLink *phi - phi * phi * DInvLink * DInvLink * xBx[0]); 
      AncestryVarScore[locus][k] += VarA[k] * phi * phi * YMinusEY * YMinusEY;
    }
    ++locus;
  }//end locus loop
}

void Individual::SumScoresForLinkageAffectedsOnly(int j, double *SumAffectedsScore, 
				      double *SumAffectedsVarScore, double *SumAffectedsScore2, double *SumAffectedsInfo){
  int KK = Populations;
  if(KK == 2) KK = 1;

  for( int k = 0; k < KK; k++ ){
    SumAffectedsScore[j*KK +k] += AffectedsScore[j*KK + k];
    SumAffectedsVarScore[j*KK +k] += AffectedsVarScore[j * KK +k];
    SumAffectedsInfo[j*KK +k] += AffectedsInfo[j * KK +k];
    SumAffectedsScore2[j*KK +k] +=  AffectedsScore[j*KK +k] * AffectedsScore[j*KK +k];
  }
}

void Individual::SumScoresForAncestry(int j, double *SumAncestryScore, double *SumAncestryInfo, double *SumAncestryScore2,
				      double *SumAncestryVarScore){

  double *score = new double[Populations], *info = new double[Populations*Populations];
  CentredGaussianConditional(Populations, AncestryScore[j], AncestryInfo[j], score, info, 2*Populations );
  
  //accumulate over iterations
  //for two populations, we only accumulate the scores for the second population
  int KK = Populations, k1 = 0;
  if(Populations == 2){KK = 1; k1 = 1;}
     
  for( int k = 0; k < KK ; k++ ){
    SumAncestryScore[j*KK +k] += score[k+k1];
    SumAncestryInfo[j*KK +k] += info[(k+k1)*Populations +k+k1] + AncestryInfoCorrection[j][k+k1];
    SumAncestryScore2[j*KK +k] += score[k+k1] * score[k+k1];
    SumAncestryVarScore[j*KK +k] += AncestryVarScore[j][k+k1];
  }
  delete[] score;
  delete[] info;
}

// unnecessary duplication of code - should use same method as for > 1 population
void Individual::OnePopulationUpdate( int i, Matrix_d *Outcome, int NumOutcomes, int* OutcomeType, double **ExpectedY, 
				      double *lambda, int AnalysisTypeIndicator )
{
  for( int k = 0; k < NumOutcomes; k++ ){
    if( AnalysisTypeIndicator > 1 ){
      if( Outcome[k].IsMissingValue( i, 0 ) ){
	if( !OutcomeType[k] )
	  Outcome[k]( i, 0 ) = gennor( ExpectedY[k][i], 1 / sqrt( lambda[k] ) );
	else{
	  if( myrand() * ExpectedY[k][i] < 1 )
	    Outcome[k]( i, 0 ) = 1;
	  else
	    Outcome[k]( i, 0 ) = 0;
	}
      }
    }
  }
  // sampled alleles should be stored in Individual objects, then summed over individuals to get counts
  // is this extra update necessary? isn't this method called anyway, irrespective of whether there is only one population?        
 //  for( int j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
//     A->UpdateAlleleCounts(j,getPossibleHaplotypes(j), ancestry );
//   }
}

void Individual::InitializeChib(double *theta, double *thetaX, vector<double> rho, vector<double> rhoX, 
				AdmixOptions *options, AlleleFreqs *A, Chromosome **chrm, double rhoalpha, double rhobeta, 
				vector<vector<double> > &alpha, chib *MargLikelihood, LogWriter *Log)
//Computes LogPrior and LogLikelihood used for Chib Algorithm
{
  int K = Populations;
  size_t theta_size = Populations;
  if(options->isRandomMatingModel()) theta_size *=2;

   double LogPrior=0, LogLikelihoodAtEst;
   Log->write("Calculating posterior at individual admixture\n");
   Log->write( theta, theta_size);Log->write("\nand sumintensities\n");Log->write( rho[0]);Log->write(rho[1]);Log->write("\n");
   // ** case of xonly data **
   if( options->getXOnlyAnalysis() ){
     LogLikelihoodAtEst = getLogLikelihoodXOnly( options, chrm, theta, rho);
     if( options->RhoFlatPrior() ){ //flat prior on sumintensities
         LogPrior = -log( options->getTruncPt() - 1.0 );
      }
     else if( options->logRhoFlatPrior() ){//flat prior on log sumintensities
         LogPrior = -log( rho[0]*(log( options->getTruncPt() ) ) );
      }
     else{//gamma prior on sumintensities
         LogPrior = getGammaLogDensity( rhoalpha, rhobeta, rho[0] );
         LogPrior -= log( gsl_cdf_gamma_Q(rhobeta, rhoalpha, 1.0) );
      }
     //alpha contribution
      LogPrior += getDirichletLogDensity( alpha[0], GetRow(theta, 0, K) );
   }
   // ** case of some x data **
   else if( Loci->isX_data() ){
     LogLikelihoodAtEst = getLogLikelihoodAtEst( options, chrm, theta, rho, thetaX, rhoX);
      if( options->RhoFlatPrior() ){//flat prior on sumintensities
         LogPrior = -4.0*log( options->getTruncPt() - 1.0 );
      }
      else if( options->logRhoFlatPrior() ){//flat prior on log sumintensities
         LogPrior = -log( rho[0]*(log( options->getTruncPt() ) ) )
            -log( rho[1]*(log( options->getTruncPt() ) ) )
            -log( rhoX[0]*(log( options->getTruncPt() ) ) )
            -log( rhoX[1]*(log( options->getTruncPt() ) ) );
      }
      else{//gamma prior on sumintensities
         LogPrior = getGammaLogDensity( rhoalpha, rhobeta, rho[0] )
            + getGammaLogDensity( rhoalpha, rhobeta, rhoX[0] )
            + getGammaLogDensity( rhoalpha, rhobeta, rho[1] )
            + getGammaLogDensity( rhoalpha, rhobeta, rhoX[1] );
         LogPrior /= gsl_cdf_gamma_Q(rhobeta, rhoalpha, 1.0);
      }
     //alpha contribution
      LogPrior += getDirichletLogDensity( alpha[0],GetRow(theta, 0, K) )
         + getDirichletLogDensity( alpha[0], GetRow(thetaX, 0, K) )
         + getDirichletLogDensity( alpha[1], GetRow(theta, 1, K) )
         + getDirichletLogDensity( alpha[1], GetRow(thetaX, 1, K) );
      LogLikelihoodAtEst = getLogLikelihoodAtEst( options, chrm, theta, rho, thetaX, rhoX );
   }
   else{
      if( Populations > 1 ){
	LogLikelihoodAtEst = getLogLikelihoodAtEst( options, chrm, theta, rho, thetaX, rhoX );
	if( options->isAdmixed(0) ){
            if( options->RhoFlatPrior() ){
               LogPrior = -log( options->getTruncPt() - 1.0 );
            }
            else if( options->logRhoFlatPrior() ){
               LogPrior = -log( rho[0]*(log( options->getTruncPt() ) ) );
            }
            else{
               LogPrior = getGammaLogDensity( rhoalpha, rhobeta, rho[0] );
               LogPrior -= log( gsl_cdf_gamma_Q(rhobeta, rhoalpha, 1.0) );
            }
            LogPrior += getDirichletLogDensity( alpha[0], GetRow(theta, 0, K) );
         }
	if( options->isAdmixed(1) ){

            if( options->RhoFlatPrior() ){
               LogPrior -= log( options->getTruncPt() - 1.0 );
            }
            else if( options->logRhoFlatPrior() ){
               LogPrior -= log( rho[1]*(log( options->getTruncPt() ) ) );
            }
            else{
               LogPrior += getGammaLogDensity( rhoalpha, rhobeta, rho[1] );
               LogPrior -= log( gsl_cdf_gamma_Q(rhobeta, rhoalpha, 1.0) );
            }
            LogPrior += getDirichletLogDensity( alpha[1], GetRow(theta, 1, K) );
         }
      }
      else{
	LogLikelihoodAtEst = getLogLikelihoodOnePop();
      }
   }
   if( A->IsRandom() ){
      for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
         for( int k = 0; k < Populations; k++ ){
	   //CompositeLocus *locus = (CompositeLocus*)(*Loci)(j);
	   LogPrior += getDirichletLogDensity( A->GetPriorAlleleFreqs(j, k), A->getAlleleFreqsMAP(j,k) );
         }
      }
   }
   MargLikelihood->setLogPrior( LogPrior );
   MargLikelihood->setLogLikelihood( LogLikelihoodAtEst );   
}

 // this function computes marginal likelihood by the Chib algorithm.  Can be replaced with 
  // more efficient algorithm based on HMM likelihood
// Chib method for marginal likelihood should be rewritten to use the HMM likelihood, without sampling locus ancestry or arrivals
void Individual::ChibLikelihood(int iteration, double *LogLikelihood, double *SumLogLikelihood, double *MaxLogLikelihood,
				AdmixOptions *options, Chromosome **chrm, vector<vector<double> > &alpha, 
				double rhoalpha, double rhobeta, double *thetahat,
				double *thetahatX, vector<double> &rhohat,
				vector<double> &rhohatX, LogWriter *Log, chib *MargLikelihood, AlleleFreqs* A){
            
  //if( iteration <= options->getBurnIn() ){
  int K = Populations;
  size_t theta_size = Populations;
  if(options->isRandomMatingModel()) theta_size *=2;
  
  *LogLikelihood = getLogLikelihood(options, chrm);//gets loglikelihood at current parameter values

    if( K > 1 ){
      if( !options->RhoFlatPrior() && !options->logRhoFlatPrior() ){
	if( options->isAdmixed(0) ){
	  *LogLikelihood+=getGammaLogDensity( rhoalpha, rhobeta, _rho[0] );}
	if(  options->isAdmixed(1) )
	  *LogLikelihood+=getGammaLogDensity( rhoalpha, rhobeta, _rho[1] );
      }
      else if( options->logRhoFlatPrior() ){
	if( options->isAdmixed(0) )
	  *LogLikelihood -= log( _rho[0] );
	if(  options->isAdmixed(1) )
	  *LogLikelihood -= log( _rho[1] );
      }
      *LogLikelihood+=
	getDirichletLogDensity(alpha[0],
			       GetRow(Theta, 0, K))
	+getDirichletLogDensity(alpha[1],
				GetRow(Theta, 1, K));
      if( *LogLikelihood > *MaxLogLikelihood ){
	Log->write("Admixture (gamete 1):");
	Log->write(Theta, K);Log->write("\n");
	Log->write("Admixture (gamete 2):");
	Log->write(Theta+K, K);Log->write("\n");
	Log->write("sumintensities: ");
	Log->write( _rho[0]);Log->write( _rho[1]);
	Log->write("\nLogLikelihood:");Log->write( *LogLikelihood);
	Log->write("\niteration: ");Log->write(iteration);Log->write("\n\n");

	*MaxLogLikelihood = *LogLikelihood;
	if( iteration <= options->getBurnIn() ){
	  for(unsigned k = 0; k < theta_size; ++k)thetahat[k] = Theta[k];
	  rhohat = _rho;
	  A->setAlleleFreqsMAP();
	  for( unsigned  j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
	    CompositeLocus *locus = (CompositeLocus*)(A->getLocus(j));
	    //locus->setAlleleFreqsMAP();
	    if( locus->GetNumberOfLoci() > 2 )
	      locus->setHaplotypeProbsMAP();
	  }
	  if( Loci->isX_data() ){
	  for(unsigned k = 0; k < theta_size; ++k)thetahatX[k] = ThetaX[k];
	    rhohatX = _rho_X;
	  }
	}
      }
    }
    else{//populations <=0
      if( *LogLikelihood > *MaxLogLikelihood ){

	*MaxLogLikelihood = *LogLikelihood;
	if( iteration <= options->getBurnIn() ){
	  A->setAlleleFreqsMAP();
	  for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){

	    //locus->setAlleleFreqsMAP();
	    if( (*Loci)(j)->GetNumberOfLoci() > 2 )
	      (*Loci)(j)->setHaplotypeProbsMAP();
	  }
	}
      }
    }
    if( options->getAnalysisTypeIndicator() == -1 ){
      if( iteration == options->getBurnIn() ){
	InitializeChib(thetahat, thetahatX, rhohat, rhohatX,
		       options, A, chrm, rhoalpha, rhobeta, 
		       alpha, MargLikelihood, Log);
      }
      if( iteration > options->getBurnIn() ){
	double LogPosterior = 0;
	if( Populations > 1 )
	  LogPosterior = getLogPosteriorProb();
	if( A->IsRandom() ){
	  for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
	    for( int k = 0; k < Populations; k++ ){
	      //CompositeLocus *locus = (CompositeLocus*)(*Loci)(j);
	      LogPosterior += getDirichletLogDensity( A->GetPriorAlleleFreqs(j,k) + A->GetAlleleCounts(j,k), 
						      A->getAlleleFreqsMAP(j, k) );
	    }
	  }
	}
	MargLikelihood->addLogPosteriorObs( LogPosterior );
	*SumLogLikelihood += *LogLikelihood;
      }
    }

}

// //TODO: need to fix this
 double 
Individual::getLogLikelihoodXOnly( AdmixOptions* options, Chromosome **chrm, double *admixture, vector<double> rho)
{
   double LogLikelihood = 0.0;
   _rhoHat = rho;
   for(int i = 0; i < AdmixtureHat.GetNumberOfRows(); ++i)
     for(int j = 0; j < AdmixtureHat.GetNumberOfCols(); ++j)
       AdmixtureHat(i,j) = admixture[j*Populations + i];

   chrm[0]->UpdateParameters( this, admixture, options, _rhoHat,  true, false, false);

   LogLikelihood += chrm[0]->getLogLikelihood();

   return LogLikelihood;
}

//computes log likelihood at parameter estimates
//This one is called by InitializeChib, called in turn by ChibLikelihood
//only used to compute marginal likelihood   
double Individual::getLogLikelihoodAtEst( AdmixOptions* options, Chromosome **chrm, double *admixture, vector<double> rho, double *admixture_X, vector<double> rho_X)
{
   double LogLikelihood = 0.0;
   _rhoHat = rho;
   if(options->isRandomMatingModel()){
     AdmixtureHat.SetNumberOfElements(Populations, 2);
     XAdmixtureHat.SetNumberOfElements(Populations, 2);
   }
   else{
     AdmixtureHat.SetNumberOfElements(Populations, 1);
     XAdmixtureHat.SetNumberOfElements(Populations, 1);
   }
    for(int i = 0; i < AdmixtureHat.GetNumberOfRows(); ++i)
     for(int j = 0; j < AdmixtureHat.GetNumberOfCols(); ++j)
       AdmixtureHat(i,j) = admixture[j*Populations + i];

    //update forward probs in HMM using current parameter values   
   for( unsigned int j = 0; j < numChromosomes; j++ ){      
      if( j != X_posn ){
	chrm[j]->UpdateParameters( this, admixture, options, _rhoHat, true, true, false);
      }
      else{//X chromosome
         _rhoHat_X = rho_X;
	 for(int i = 0; i < XAdmixtureHat.GetNumberOfRows(); ++i)
	   for(int jj = 0; jj < XAdmixtureHat.GetNumberOfCols(); ++jj)
	     XAdmixtureHat(i,jj) = admixture_X[jj*Populations + i];
         if( sex == 1 ){//male
	   chrm[j]->UpdateParameters( this, admixture_X, options, _rhoHat_X, true, false, false);
	 }
         else{//female
	   chrm[j]->UpdateParameters( this, admixture_X, options, _rhoHat_X, true, true, false);
	 }
      }
      //accumulate loglikelihood from HMM
      LogLikelihood += chrm[j]->getLogLikelihood();
   }

   return LogLikelihood;
}

double Individual::getLogLikelihoodOnePop()
{
   double Likelihood = 0.0;
   double *Prob;
   Prob = new double[1];//one pop so 1x1 array
   for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
     if(!IsMissing(j)){
       (*Loci)(j)->GetGenotypeProbs(Prob,getPossibleHapPairs(j), true);
       Likelihood += log( Prob[0] );
     }
   }
   return Likelihood;
}

//updates forward probs in HMM and retrieves loglikelihood
double Individual::getLogLikelihood( AdmixOptions* options, Chromosome **chrm)
{
   double LogLikelihood = 0.0;

   for( unsigned int j = 0; j < numChromosomes; j++ ){      
      if( j != X_posn ){

	chrm[j]->UpdateParameters( this, Theta, options, _rho, false, true, false);
      }
      else if( options->getXOnlyAnalysis() ){

	chrm[j]->UpdateParameters( this, Theta, options, _rho, false, false, false);
      }
      else{//X chromosome
	if( sex == 1 ){//male
	  chrm[j]->UpdateParameters( this, Theta, options, _rho_X, false, false, false);
	}
	else{//female
	  chrm[j]->UpdateParameters( this, Theta, options, _rho_X, false, true, false);
	}
      }
      LogLikelihood += chrm[j]->getLogLikelihood();
     }

   return LogLikelihood;
}

void Individual::CalculateLogPosterior(AdmixOptions *options, bool isX_data, vector<vector<double> > &alpha, 
				       double rhoalpha, double rhobeta, double L, 
				       double L_X, unsigned int SumN[],unsigned int SumN_X[]){

  LogPosterior = 0.0; 
  double IntConst1;
 {
   vector<double> alphaparams1(Populations), alphaparams0(Populations);
    if( options->getXOnlyAnalysis() ){
      LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN_X[0],
					  rhobeta + L_X, _rhoHat[0] );
      if(!options->RhoFlatPrior() && !options->logRhoFlatPrior() )
	IntConst1 = IntegratingConst(rhoalpha+(double)SumN_X[0], rhobeta+L_X, 1.0, options->getTruncPt() );
      else
	IntConst1 = gsl_cdf_gamma_Q(rhobeta+L_X, rhoalpha+(double)SumN_X[0], 1.0);
      LogPosterior -= log(IntConst1);
      //      alphaparams0 = alpha[0] + GetRow(SumLocusAncestry_X, 0 , Populations);
      transform(alpha[0].begin(), alpha[0].end(), SumLocusAncestry, alphaparams0.begin(), std::plus<double>());
      LogPosterior += getDirichletLogDensity(alphaparams0,AdmixtureHat.GetColumn(0));
    }
    else if( isX_data ){
      for( unsigned int g = 0; g < 2; g++ ){
	LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN[g],
					    rhobeta + L, _rhoHat[g] );
	if(!options->RhoFlatPrior() && !options->logRhoFlatPrior() )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN[g], rhobeta+L, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[g], 1.0);
	LogPosterior -= log(IntConst1);
	//alphaparams0 = alpha[g] + GetRow(SumLocusAncestry, g, Populations);
	transform(alpha[g].begin(), alpha[g].end(), SumLocusAncestry +g*Populations, alphaparams0.begin(), std::plus<double>());
	LogPosterior += getDirichletLogDensity(alphaparams0,AdmixtureHat.GetColumn(g));
      }
      for( unsigned int g = 0; g < gametes[X_posn]; g++ ){
	LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN_X[g],
					    rhobeta + L_X, _rhoHat_X[g] );
	if( options->RhoFlatPrior() || options->logRhoFlatPrior() )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN_X[g], rhobeta+L_X, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L_X, rhoalpha+(double)SumN_X[g], 1.0);
	LogPosterior -= log(IntConst1);
	//alphaparams0 = alpha[g] + GetRow(SumLocusAncestry_X, g, Populations);
	transform(alpha[g].begin(), alpha[g].end(), SumLocusAncestry_X +g*Populations, alphaparams0.begin(), std::plus<double>());
	LogPosterior += getDirichletLogDensity(alphaparams0,XAdmixtureHat.GetColumn(g));
      }
    }
    else if( options->isSymmetric() ){
      vector<double> x(2,0.0);
      x[0] += getGammaLogDensity( rhoalpha + (double)SumN[0],
				  rhobeta + L, _rhoHat[0] );
      x[1] += getGammaLogDensity( rhoalpha + (double)SumN[1],
				  rhobeta + L, _rhoHat[0] );
      x[0] += getGammaLogDensity( rhoalpha + (double)SumN[1],
				  rhobeta + L, _rhoHat[1] );
      x[1] += getGammaLogDensity( rhoalpha + (double)SumN[0],
				  rhobeta + L, _rhoHat[1] );
      double IntConst1, IntConst2;
      if( options->RhoFlatPrior() || options->logRhoFlatPrior() ){
	IntConst1 = IntegratingConst(rhoalpha+(double)SumN[0], rhobeta+L, 1.0, options->getTruncPt() );
	IntConst2 = IntegratingConst(rhoalpha+(double)SumN[1], rhobeta+L, 1.0, options->getTruncPt() );
      }
      else{
	IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[0], 1.0);
	IntConst2 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[1], 1.0);
      }
      x[0] -= log(IntConst1);
      x[1] -= log(IntConst2);
      //alphaparams0 = alpha[0] + GetRow(SumLocusAncestry, 0, Populations);
      transform(alpha[0].begin(), alpha[0].end(), SumLocusAncestry, alphaparams0.begin(), std::plus<double>());
      x[0] += getDirichletLogDensity(alphaparams0,AdmixtureHat.GetColumn(0));
      x[1] += getDirichletLogDensity(alphaparams0,AdmixtureHat.GetColumn(1));
      //alphaparams1 = alpha[1] + GetRow(SumLocusAncestry, 1, Populations);
      transform(alpha[1].begin(), alpha[1].end(), SumLocusAncestry+Populations, alphaparams1.begin(), std::plus<double>());
      x[0] += getDirichletLogDensity(alphaparams1,AdmixtureHat.GetColumn(1));
      x[1] += getDirichletLogDensity(alphaparams1,AdmixtureHat.GetColumn(0));
      if( isnan(x[0]) || isinf(x[0]) )
	LogPosterior = x[1] - log(2.0);
      else if( isnan(x[1]) || isinf(x[1])  )
	LogPosterior = x[0] - log(2.0);
      else if( x[0] < x[1] )
	LogPosterior = x[1] + log( 1 + exp( x[0] - x[1] ) ) - log(2.0);
      else
	LogPosterior = x[0] + log( 1 + exp( x[1] - x[0] ) ) - log(2.0);
      //if( isnan(LogPosterior) || isinf(LogPosterior) ){
      //PR(alphaparams0);
      //PR(alphaparams1);
      //PR(AdmixtureHat.GetColumn(0));
      //   PR(AdmixtureHat.GetColumn(1));
      //   PR(x[0]);
      //   PR(x[1]);
      //   exit(0);
      //}
    }
    else{
      if(  options->isAdmixed(0) ){
	LogPosterior = getGammaLogDensity( rhoalpha + (double)SumN[0],
					   rhobeta + L, _rhoHat[0] );
	if( options->RhoFlatPrior() || options->logRhoFlatPrior() )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN[0], rhobeta+L, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[0], 1.0);
           LogPosterior -= log( IntConst1 );
           //alphaparams0 = alpha[0] + GetRow(SumLocusAncestry, 0, Populations);
	   transform(alpha[0].begin(), alpha[0].end(), SumLocusAncestry, alphaparams0.begin(), std::plus<double>());
           LogPosterior+=getDirichletLogDensity(alphaparams0,AdmixtureHat.GetColumn(0));
      }
      if(  options->isAdmixed(1) ){
	LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN[1],
					    rhobeta + L, _rhoHat[1] );
	if( options->RhoFlatPrior() || options->logRhoFlatPrior() )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN[1], rhobeta+L, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[1], 1.0);
	LogPosterior -= log( IntConst1 );
	//alphaparams1 = alpha[1] + GetRow(SumLocusAncestry, 1, Populations);
	transform(alpha[1].begin(), alpha[1].end(), SumLocusAncestry+Populations, alphaparams1.begin(), std::plus<double>());
	LogPosterior+=getDirichletLogDensity(alphaparams1,AdmixtureHat.GetColumn(1));
      }
    }
  }
}

double Individual::IntegratingConst( double alpha, double beta, double a, double b )
{
   double I = gsl_cdf_gamma_P( b*beta, alpha, 1 ) - gsl_cdf_gamma_P( a*beta, alpha, 1);
   return I;
}

