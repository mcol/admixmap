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
#include <sstream>

#define PR(x) cout << #x << " = " << x << endl;

// ******* Static Member Declarations
double *Individual::LikRatio1;
double *Individual::LikRatio2;
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
int Individual::Populations;

//******** Constructors **********
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
	_rho.assign(1,1);
  }
  sumlogrho.assign(_rho.size(), 0.0);
  _rhoHat.assign(_rho.size(), 0.0);
  _rhoHat_X.assign(_rho.size(), 0.0);
  
  // Read sex value if present.
  sex = male;
  if (options->getgenotypesSexColumn() == 1) {
    sex = Data->GetSexValue(mynumber);
  }
  
  int numCompositeLoci = Loci.GetNumberOfCompositeLoci();
  
  LocusAncestry = new int*[ numChromosomes ]; // array of matrices in which each col stores 2 integers 
  
  Theta = 0;
  ThetaX = 0;
  ThetaProposal = 0;
  ThetaXProposal = 0;
  ThetaHat = 0;
  ThetaXHat = 0;
  SumSoftmaxTheta = 0;
  
  // SumLocusAncestry is sum of locus ancestry states over loci at which jump indicator xi is 1  
  SumLocusAncestry = new int[options->getPopulations()*2];
  
  
  if( options->isRandomMatingModel() ){//random mating model
    ThetaProposal = new double[ Populations * 2 ];
    Theta = new double[ Populations * 2 ];
    SumSoftmaxTheta = new double[ Populations * 2 ];
    fill(SumSoftmaxTheta, SumSoftmaxTheta + Populations*2, 0.0);
    ThetaHat = new double[ Populations * 2 ];;
  }
  else{
    ThetaProposal = new double[ Populations];
    Theta = new double[ Populations ];
    SumSoftmaxTheta = new double[ Populations ];
    fill(SumSoftmaxTheta, SumSoftmaxTheta + Populations, 0.0);
    ThetaHat = new double[ Populations ];
  }
  
  //X chromosome objects
  SumLocusAncestry_X = 0;    
  //    if(Loci.isX_data() ){
  if( sex == male ){
    ThetaXProposal = new double[ Populations];
    ThetaX = new double[Populations];
    ThetaXHat = new double[Populations];
    SumLocusAncestry_X = new int[Populations];
  }
  else{
    ThetaXProposal = new double[ Populations * 2 ];
    ThetaX = new double[ Populations * 2 ];
    ThetaXHat = new double[ Populations * 2 ];
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
    else if( sex != female ){//male or missing
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
  
  // ** set up StepSizeTuner object for random walk updates of admixture **
  NumberOfUpdates = 0;
  w = 1;
  step0 = 1.0; // initial sd of proposal distribution
  //need to choose sensible value for this initial RW sd
  step = step0;
  ThetaTuner.SetParameters( step0, 0.00, 10.0, 0.44);  
   
}

//********** Destructor **********
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
  delete[] SumSoftmaxTheta;
}

//********** Allocation and deletion of static objects for score tests
void Individual::SetStaticMembers(Genome *pLoci, AdmixOptions *options){
  Loci = pLoci;
  numChromosomes = Loci->GetNumberOfChromosomes();
  Populations = options->getPopulations();
  int L = Loci->GetNumberOfCompositeLoci();

  LikRatio1 = 0;
  LikRatio2 = 0;
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
    LikRatio1 = new double[L*K];
    LikRatio2 = new double[L*K];
    fill(LikRatio1, LikRatio1+L*K, 0.0);
    fill(LikRatio2, LikRatio2+L*K, 0.0);
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

void Individual::DeleteStaticMembers(){
  delete[] AffectedsScore;
  delete[] AffectedsInfo;
  delete[] AffectedsVarScore;
  delete[] LikRatio1;
  delete[] LikRatio2;
  free_matrix(AncestryScore, Loci->GetNumberOfCompositeLoci());
  free_matrix(AncestryInfo, Loci->GetNumberOfCompositeLoci());  
  free_matrix(AncestryVarScore, Loci->GetNumberOfCompositeLoci());
  free_matrix(AncestryInfoCorrection, Loci->GetNumberOfCompositeLoci());
}
//********** Initialisation of Admixture porportions *********

//should call this initialise not set
void Individual::setAdmixtureProps(double *a, size_t size)
{
  for(unsigned i = 0; i < size; ++i)  Theta[i] = a[i];
}

void Individual::setAdmixturePropsX(double *a, size_t size)
{
  for(unsigned i = 0; i < size; ++i)  ThetaX[i] = a[i];
}

//******************** Accessors ***********************************************************
unsigned short **Individual::getGenotype(unsigned int locus){
  return genotypes[locus];
}

std::vector<hapPair > &Individual::getPossibleHapPairs(unsigned int locus){
  return PossibleHapPairs[locus];
}

double *Individual::getAdmixtureProps()
{
  return Theta;
}

Sex Individual::getSex()
{
   return sex;
}


double Individual::getSumrho()
//returns sum of sumintensities over gametes
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

void Individual::GetLocusAncestry(int chrm, int locus, int Ancestry[2]){
  Ancestry[0]  = LocusAncestry[chrm][locus];
  if((unsigned)chrm == X_posn)Ancestry[1] = Ancestry[0];
  else Ancestry[1] = LocusAncestry[chrm][Loci->GetSizesOfChromosomes()[chrm]  + locus];
}

//returns value of LocusAncestry at a locus for a particular gamete
int Individual::GetLocusAncestry(int chrm, int gamete, int locus){
  int g = (gametes[chrm] == 2) ? gamete : 0; //so that gamete = 1 works when gametes[chrm] = 1;
  return LocusAncestry[chrm][g * Loci->GetSizesOfChromosomes()[chrm]  + locus] ;
}

int *Individual::getSumLocusAncestry(){
  return SumLocusAncestry;
}

//********** Missing Genotype Indicator *********************

//Indicates whether genotype is missing at all simple loci within a composite locus
bool Individual::IsMissing(unsigned int locus)
{
  unsigned int count = 0;
  int NumberOfLoci = (*Loci)(locus)->GetNumberOfLoci();
  for(int i = 0; i < NumberOfLoci; i++){
    count += genotypes[locus][i][0];
  }
  return (count == 0);
}

//****************** Log-Likelihoods **********************

// //TODO: need to fix this
double Individual::getLogLikelihoodXOnly( AdmixOptions* options, Chromosome **chrm, double *admixture, vector<double> rho)
{
   double LogLikelihood = 0.0;
   _rhoHat = rho;
   size_t size = Populations;
   if(options->isRandomMatingModel())size *= 2;
   copy(admixture, admixture + size, ThetaHat); 

   chrm[0]->UpdateParameters( this, admixture, options, _rhoHat,  true, false, false);

   LogLikelihood += chrm[0]->getLogLikelihood();

   return LogLikelihood;
}

double Individual::getLogLikelihoodOnePop(bool chibindicator)
//Single population log-likelihood
{
   double Likelihood = 0.0;
   double *Prob;
   Prob = new double[1];//one pop so 1x1 array
   for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
     if(!IsMissing(j)){
       (*Loci)(j)->GetGenotypeProbs(Prob,getPossibleHapPairs(j), chibindicator);
       Likelihood += log( Prob[0] );
     }
   }
   return Likelihood;
}

double Individual::getLogLikelihood( AdmixOptions* options, Chromosome **chrm){
  //use current parameter values
  return getLogLikelihood(options, chrm, Theta, ThetaX,_rho, _rho_X, false);
}

double Individual::getLogLikelihoodAtPosteriorMeans(AdmixOptions* options, Chromosome **chrm){
  //TODO: X chromosome objects
  //obtain ergodic averages of (softmax)admixture props and (log)sumintensities and transform back
  //to original scales
  unsigned size = Populations; if(options->isRandomMatingModel())size *=2;
  for(unsigned i = 0; i < size; ++i)SumSoftmaxTheta[i] /= (options->getTotalSamples() - options->getBurnIn());
  for(unsigned i = 0; i < _rho.size(); ++i)sumlogrho[i] = exp(sumlogrho[i]/(options->getTotalSamples() - options->getBurnIn()));

  //apply softmax transformation to obtain thetabar
  unsigned G = 1;
  if( options->isRandomMatingModel() )G = 2;//random mating model
  double ThetaBar[G*Populations];
  for( unsigned int g = 0; g < G; g++ ){
    bool b[Populations];
    for(int k = 0; k < Populations; ++k)if(Theta[g*Populations + k] > 0.0){
      b[k] = true; //to skip elements set to zero
    }
    else b[k] = false;

    softmax(Populations, ThetaBar+g*Populations, SumSoftmaxTheta+g*Populations, b);
  }

  return getLogLikelihood(options, chrm, ThetaBar, ThetaBar, sumlogrho, sumlogrho, false);
}
double Individual::getLogLikelihood( AdmixOptions* options, Chromosome **chrm, double *theta, double *thetaX,
				     vector<double > rho, vector<double> rho_X, bool chibindicator = false)
//updates forward probs in HMM and retrieves loglikelihood at supplied theta and rho
//(optional)chibindicator = true for computing LogL at MLEs; instructs CompositeLocus to use HapPairProbsMAP
//instead of HapPairProbs, when allelefreqs are not fixed, in calculating GenotypeProbs.
{

   double LogLikelihood = 0.0;
   if(Populations == 1)LogLikelihood = getLogLikelihoodOnePop(chibindicator);//in case this version called when only one population
   else{
     for( unsigned int j = 0; j < numChromosomes; j++ ){      
       UpdateForBackProbs(j, chrm[j], options, theta, thetaX, rho, rho_X, chibindicator, false);
       LogLikelihood += chrm[j]->getLogLikelihood();
     }
   }

   return LogLikelihood;
}

//************** Updating (Public) ***************************************************************************************

// unnecessary duplication of code - ? should embed within method for > 1 population
void Individual::OnePopulationUpdate( int i, DataMatrix *Outcome, int NumOutcomes, DataType* OutcomeType, double **ExpectedY,
				      double *lambda, Chromosome **chrm, AlleleFreqs *A )
{
  // sample missing values of outcome variable
  for( int k = 0; k < NumOutcomes; k++ ){
       if( Outcome->isMissing( i, k ) ){
	if( OutcomeType[k] == Continuous )
	  Outcome->set( i, k, gennor( ExpectedY[k][i], 1 / sqrt( lambda[k] ) ));
	else{
	  if( myrand() * ExpectedY[k][i] < 1 )
	    Outcome->set( i, k, 1);
	  else
	    Outcome->set( i, k, 0);
	}
      }
  }
  // update allele counts
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    //loop over loci on current chromosome and update allele counts
    for( unsigned int jj = 0; jj < chrm[j]->GetSize(); jj++ ){
      int locus =  chrm[j]->GetLocus(jj);
      if( !(IsMissing(j)) ){
	int anc[2] = {0, 0}; //ancestry states for single population
	// GetLocusAncestry(j,jj,anc);
	int h[2]; //to store sampled hap pair
	(*Loci)(locus)->SampleHapPair(h, PossibleHapPairs[locus], anc);
	A->UpdateAlleleCounts(locus, h, anc, true); // should fix this to work with haploid data: last argument should be isdiploid
      }
    }
  }   

}

void Individual::SampleParameters( int i, double *SumLogTheta, double *LogLikelihood, AlleleFreqs *A, int iteration , DataMatrix *Outcome,
				  int NumOutcomes, DataType* OutcomeType, double **ExpectedY, double *lambda, int NoCovariates,
				   DataMatrix *Covariates, double **beta, const double *poptheta,
				   AdmixOptions* options, Chromosome **chrm, 
				   vector<vector<double> > &alpha, double rhoalpha, 
				   double rhobeta,vector<double> sigma, double DInvLink, double dispersion)
/*arguments:
  i = this individuals's number-1
  SumLogTheta = array in IndividualCollection holding sums of log admixture props
  LogLikelihood = pointer to sum of individual LogLikelihoods (dould do this as static variable)
  AlleleFreqs = pointer to AlleleFreqs, needed to update allele counts
  iteration = current iteration
  Outcome = Outcome variable(s)
  NumOutcomes = number of outcomes
  OutcomeType = array of types of outcome (binary/continuous)   
  ExpectedY = expected outcome variable
  lambda = precision in linear regression model (if there is one)
  NoCovariates = # covariates, including admixture
  Covariates = Covariates array, including admixture
  beta = regression parameters
  poptheta = ergodic average of population admixture, used to centre the values of individual admixture in the regression model
  options = pointer to program options
  chrm = array of Chromosome pointers
  alpha = pop admixture Dirichlet parameters
  rhoalpha, rhobeta = shape and scale parameters in prior for rho
  sigma = ?? used in X-chromosome update
  DInvLink = Derivative Inverse Link function in regression model, used in ancestry score test
  dispersion = dispersion parameter in regression model (if there is one) = lambda for linear reg, 1 for logistic
*/
{
  // ** reset SumLocusAncestry and ThetaProposal **
  for(int j = 0; j < Populations *2; ++j)SumLocusAncestry[j] = 0;
//  if(Loci->isX_data() ){
    int J = Populations;
    if(sex != male) J *=2;
    for(int j = 0; j < J ;++j)SumLocusAncestry_X[j] = 0;
//  }

  size_t size_theta;
   if( options->isRandomMatingModel() )
     size_theta = Populations*2; // double the size for 2 gametes in RMM
   else//assortative mating
     size_theta = Populations;

   for(unsigned k = 0; k < size_theta; ++k){
     ThetaProposal[k] = 0.0;//may be unnecessary
   }
   if(ThetaXProposal){
     size_theta = Populations;if(sex != male) size_theta *= 2;
     for(unsigned k = 0; k < size_theta; ++k) ThetaXProposal[k] = 0.0;
   }

   if(!(iteration %2))//update theta with random walk on odd-numbered iterations
     SampleTheta(i, iteration, SumLogTheta,Outcome, chrm, NumOutcomes, OutcomeType, ExpectedY, lambda, NoCovariates,
 		 Covariates, beta, poptheta, options, alpha, sigma, DInvLink, dispersion, true);

  //SumN is the number of arrivals between each pair of adjacent loci
   SumN[0] = SumN[1] = 0;
   SumN_X[0] = SumN_X[1] = 0;  
 
  bool isdiploid;
  bool calcbackprobs = (options->getTestForAffectedsOnly() || options->getTestForLinkageWithAncestry());
  for( unsigned int j = 0; j < numChromosomes; j++ ){

    //Update Forward/Backward probs in HMM
    isdiploid = UpdateForBackProbs(j, chrm[j], options, Theta, ThetaX, _rho, _rho_X, false, calcbackprobs);

    //update score tests for linkage with ancestry for *previous* iteration
    if(iteration > options->getBurnIn()){
      //Update affecteds only scores
      if( options->getTestForAffectedsOnly())
	if(options->getNumberOfOutcomes() == 0 || Outcome->get(i, 0) == 1)
	UpdateScoreForLinkageAffectedsOnly(j, options->isRandomMatingModel(),chrm );

      //update ancestry score tests
      if( options->getTestForLinkageWithAncestry() ){   
	UpdateScoreForAncestry(j,dispersion, Outcome->get(i, 0) - ExpectedY[0][i], DInvLink,chrm);
      }    
    }

    //Sample locus ancestry
    // sampling locus ancestry requires calculation of forward probability vectors alpha in HMM 
    chrm[j]->SampleLocusAncestry(LocusAncestry[j], Theta, isdiploid);
 
   //loop over loci on current chromosome and update allele counts
    for( unsigned int jj = 0; jj < chrm[j]->GetSize(); jj++ ){
      int locus =  chrm[j]->GetLocus(jj);
      if( !(IsMissing(j)) ){
	  int anc[2];//to store ancestry states
	  GetLocusAncestry(j,jj,anc);
	  int h[2];//to store sampled hap pair
	  //might be a shortcut for haploid data since there is only one compatible hap pair, no need to sample
	  (*Loci)(locus)->SampleHapPair(h, PossibleHapPairs[locus], anc);
	  A->UpdateAlleleCounts(locus,h, anc, isdiploid);
      }
     }   

    //sample number of arrivals and SumLocusAncestry
    bool isX = (j == X_posn);
    chrm[j]->SampleJumpIndicators(LocusAncestry[j], gametes[j], SumLocusAncestry, SumLocusAncestry_X, isX,
				  SumN, SumN_X, options->getRhoIndicator());

    //accumulate LogLikelihood from HMM
    *LogLikelihood += chrm[j]->getLogLikelihood();
  }//end chromosome loop

  // sample sum of intensities parameter rho
  if( options->getRhoIndicator() ){
     SampleRho( options->getXOnlyAnalysis(), options->isRandomMatingModel(), Loci->isX_data(), rhoalpha, rhobeta, 
	       SumN, SumN_X);
  }
  if(iteration > options->getBurnIn()){
    for(unsigned i = 0; i < _rho.size(); ++i)sumlogrho[i] += log(_rho[i]);
  }

  if( options->getNumberOfOutcomes() > 0 ){
    // sample missing values of outcome variable
    double u;
    for( int k = 0; k < NumOutcomes; k++ ){
      if( Outcome->isMissing( i, k ) ){
	if( OutcomeType[k] == Continuous ) // linear regression
	  Outcome->set( i, k, gennor( ExpectedY[k][i], 1 / sqrt( lambda[k] ) ));
	else{// logistic regression
	  u = myrand();
	  if( u * ExpectedY[k][i] < 1 )
	    Outcome->set( i, k, 1);
	  else
	    Outcome->set( i, k, 0);
	}
      }
    }
  }

}
// ****** End Public Interface *******

//************** Updating (Private) ***************************************************************************************

void Individual::SampleTheta( int i, int iteration, double *SumLogTheta, DataMatrix *Outcome, Chromosome ** C,
                                 int NumOutcomes, DataType* OutcomeType, double **ExpectedY, double *lambda, int NoCovariates,
                                  DataMatrix *Covariates, double **beta, const double *poptheta,
			      AdmixOptions* options, vector<vector<double> > &alpha, vector<double> sigma,
			      double DInvLink, double dispersion, bool RW)
// samples individual admixture proportions
{
  double logpratio = 0;

  // propose new value for individual admixture proportions
  // should be modified to allow a population mixture component model
  if(RW) {
    NumberOfUpdates++;

    logpratio += ProposeThetaWithRandomWalk(options, C, alpha);
  }
  else 
  ProposeTheta(options, sigma, alpha);       

  int K = Populations;

  //calculate Metropolis acceptance probability ratio for proposal theta    
  if(!options->getTestForAdmixtureAssociation()){
    RegressionType RegType;
    for( int k = 0; k < NumOutcomes; k++ ){
      if(OutcomeType[k] == Binary)RegType = Logistic; else RegType = Linear;
      logpratio +=  LogAcceptanceRatioForRegressionModel( i, RegType, k, options->isRandomMatingModel(), K,
							 NoCovariates, Covariates, beta, ExpectedY, Outcome, poptheta,lambda);
    }
  }
  //case of X only data
  if( Loci->isX_data() && !options->getXOnlyAnalysis() )
    logpratio += AcceptanceProbForTheta_XChrm( sigma, K);

 //Accept or reject proposed value - if no regression model, proposal will be accepted because logpratio = 0
  Accept_Reject_Theta(logpratio, Loci->isX_data(), K, options->isRandomMatingModel(), RW );

 // update the value of admixture proportions used in the regression model  
  if( options->getNumberOfOutcomes() > 0 )
    UpdateAdmixtureForRegression(i, K, NoCovariates, poptheta, options->isRandomMatingModel(), Covariates);

  if(iteration > options->getBurnIn()){
    unsigned G = 1;
    if( options->isRandomMatingModel() )G = 2;//random mating model
 
    for( unsigned int g = 0; g < G; g++ ){
      bool b[Populations];
      for(int k = 0; k < Populations; ++k)if(Theta[g*Populations + k] > 0.0){
	b[k] = true; //to skip elements set to zero
      }
      else b[k] = false;

      double a[Populations];
      inv_softmax(Populations, Theta+g*Populations, a, b);

      transform(a, a+Populations, SumSoftmaxTheta+g*Populations, SumSoftmaxTheta+g*Populations, std::plus<double>());
      //transform(ThetaX, ThetaX+size_admix, ThetaXHat, ThetaXHat, std::plus<double>());
    }
  }

  for( int k = 0; k < K; k++ ){
    SumLogTheta[ k ] += log( Theta[ k ] );
      if(options->isRandomMatingModel() && !options->getXOnlyAnalysis() )
	SumLogTheta[ k ] += log( Theta[ K + k ] );
    }

//   //increment B using new Admixture Props
//   //Xcov is a vector of admixture props as covariates as in UpdateScoreForAncestry
   if(iteration >= options->getBurnIn() && options->getTestForLinkageWithAncestry()){
     UpdateB(DInvLink, dispersion);
   }

}

double Individual::ProposeThetaWithRandomWalk(AdmixOptions *options, Chromosome **C, vector<vector<double> > &alpha){
  //TODO: X-chromosome case
    double LogLikelihoodRatio = 0.0;
    double LogPriorRatio = 0.0;

    //cout<<"Stepsize = "<<step<<", AccRate = "<<ThetaTuner.getExpectedAcceptanceRate()<<endl;

    //generate proposals
    unsigned G = 1;
    if( options->isRandomMatingModel() )G = 2;//random mating model
    for( unsigned int g = 0; g < G; g++ ){
      //perform softmax transformation
      bool b[Populations];
      int last = -1;//index of last nonzero element of theta 
      for(int k = 0; k < Populations; ++k)if(Theta[g*Populations + k] > 0.0){
	b[k] = true; //to skip elements set to zero
	++last;//last should be >0 unless unadmixed on this gamete
      }
      else b[k] = false;

      double a[Populations];
      inv_softmax(Populations, Theta+g*Populations, a, b);

      //cout<<"a(current) = ";
      a[Populations-1] = 0.0;
	//random walk step
      for(int k = 0; k < last; ++k)if(Theta[g*Populations + k]>0.0){
	//cout<<a[k]<<" ";
	a[k] = gennor(a[k], step0);
	a[last] -= a[k];
      }

      //reverse transformation
      softmax(Populations, ThetaProposal+g*Populations, a, b);
      
      //cout<<endl<<"a(prop) = "<<a[0]<<" "<<a[1]<<endl;

      //compute prior ratio
      LogPriorRatio += getDirichletLogDensity_Softmax(alpha[0], ThetaProposal+g*Populations) - 
	getDirichletLogDensity_Softmax(alpha[0], Theta+g*Populations);
    }
    //get log likelihood at current parameter values
    LogLikelihoodRatio -= getLogLikelihood(options, C);
 
    //cout<<"Current theta: "<<endl;
    //for(int k = 0; k < Populations; ++k)cout<<Theta[k]<<" ";
    //cout<<endl;
    //cout<<"Proposal theta: "<<endl;
    //for(int k = 0; k < Populations; ++k)cout<<ThetaProposal[k]<<" ";
    //cout<<endl;
    //get log likelihood at proposal theta and current rho
    LogLikelihoodRatio += getLogLikelihood(options, C, ThetaProposal, ThetaXProposal,_rho, _rho_X, false);
    //cout<<"LogLikelihood ratio = "<<LogLikelihoodRatio<<endl
    //<<"Log Acceptance Prob = "<<LogLikelihoodRatio + LogPriorRatio<<endl;

    return LogLikelihoodRatio + LogPriorRatio;
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
      if( options->getIndAdmixHierIndicator() ){
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

// returns log of ratio of likelihoods of new and old values of population admixture
// in regression models.  individual admixture theta is standardized about the mean poptheta calculated during burn-in. 
// returns log of ratio of likelihoods of new and old values of population admixture
// in regression models.  individual admixture theta is standardized about the mean poptheta calculated during burn-in. 
double Individual::LogAcceptanceRatioForRegressionModel( int i, RegressionType RegType, int TI,  bool RandomMatingModel, int Populations,
							 int NoCovariates, DataMatrix *Covariates, double **beta, double **ExpectedY,
							 DataMatrix *Outcome, const double *poptheta, double *lambda)
{
  double logprobratio = 0.0, Xbeta = 0.0;
   double avgtheta[Populations];avgtheta[0] = 0.0;
  if( RandomMatingModel )
    for(int k = 1;k < Populations; ++k)avgtheta[k] = (ThetaProposal[k] + ThetaProposal[k + Populations ])/ 2.0 - poptheta[k];
  else
    for(int k = 1;k < Populations; ++k)avgtheta[k] = ThetaProposal[k]  - poptheta[k];

  for( int jj = 0; jj < NoCovariates - Populations + 1; jj++ )
    Xbeta += Covariates->get( i, jj ) * beta[ TI ][jj];
  for( int k = 1; k < Populations; k++ ){
    Xbeta += avgtheta[ k ] * beta[ TI ][NoCovariates - Populations + k ];
  }
  if(RegType == Linear){
    logprobratio = 0.5 * lambda[ TI ] * (( ExpectedY[ TI ][i] - Outcome->get( i, TI ) ) * ( ExpectedY[ TI ][i] - Outcome->get( i, TI ) )
					 - ( Xbeta - Outcome->get( i, TI ) ) * ( Xbeta - Outcome->get( i, TI) ) );
  }
  else if(RegType == Logistic){
    double newExpectedY = 1.0 / ( 1.0 + exp( -Xbeta ) );
    if( Outcome->get( i, TI ) == 1 )
      logprobratio = newExpectedY / ExpectedY[ TI ][i];
    else
      logprobratio = ( 1 - newExpectedY ) / ( 1 - ExpectedY[ TI ][i] );
    logprobratio = log(logprobratio);//what is actually calculated above is the prob ratio. We take the log here rather than compute 4 logs above
  }
  return( logprobratio );
}
double Individual::AcceptanceProbForTheta_XChrm(std::vector<double> &sigma, int Populations )
{
   int gametes = 1;
   if( sex == female )
      gametes = 2;
   double p = 0, sum1 = 0.0, sum2 = 0.0;

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

// update the individual admixture values (mean of both gametes) used in the regression model
void Individual::UpdateAdmixtureForRegression( int i,int Populations, int NoCovariates,
                                               const double *poptheta, bool RandomMatingModel, DataMatrix *Covariates)
{
  double avgtheta[Populations];
  if(RandomMatingModel )//average over gametes
    for(int k = 0; k < Populations; ++k) avgtheta[k] = (Theta[k] + Theta[k + Populations]) / 2.0;    
  
  else
    for(int k = 0; k < Populations; ++k) avgtheta[k] = Theta[k];
  for( int k = 1; k < Populations ; k++ )
    Covariates->set( i, NoCovariates - Populations + k, avgtheta[ k ] - poptheta[ k ] );
}

// Metropolis update for admixture proportions theta, taking log of acceptance probability ratio as argument
// uses log ratio because this is easier than ratio to calculate for linear regression model
// if no regression model, logpratio remains set to 0, so all proposals are accepted 
void Individual::Accept_Reject_Theta( double logpratio, bool xdata, int Populations, bool RandomMatingModel, bool RW )
{
  bool test = true;
  double AccProb = exp(logpratio);
  // loop over populations: if any element of Dirichlet parameter vector is too small, do not update admixture proportions
  for( int k = 0; k < Populations; k++ ){
    if( Theta[ k ] > 0.0 && ThetaProposal[ k ] < 0.0001 ){
      test = false;
      AccProb = 0.0;
    }
    else if( RandomMatingModel && Theta[k + Populations] > 0.0 && ThetaProposal[ k + Populations ] < 0.0001 ){
      test = false;
      AccProb = 0.0;
    }
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
  else{//logpratio >= 0 => always accept
    if(test){
    AccProb = 1.0;
    setAdmixtureProps(ThetaProposal, size_admix);
     if( xdata )
       setAdmixturePropsX(ThetaXProposal, size_admix);
    }
  }

  if(RW){
    //update sampler object every w updates
    if( !( NumberOfUpdates % w ) ){
      step = ThetaTuner.UpdateStepSize( AccProb );
    }
  }
}

//Updates forward and backward probabilities in HMM for chromosome j 
bool Individual::UpdateForBackProbs(unsigned int j, Chromosome *chrm, AdmixOptions *options, 
				    double* theta, double *thetaX,
				    vector<double> rho, vector<double> rhoX, bool chibindicator, bool calcbackprobs){
  bool isdiploid;

  //Update Forward/Backward probs in HMM
  if( j != X_posn ){
    chrm->UpdateParameters(this, theta, options, rho, chibindicator, true, calcbackprobs);
    isdiploid = true;
  }
  else if( options->getXOnlyAnalysis() ){
    chrm->UpdateParameters( this, theta, options, rho, chibindicator, false, calcbackprobs);
    isdiploid = false;
  }
  else if( sex == male ){
    chrm->UpdateParameters( this, thetaX, options, rhoX, chibindicator, false, calcbackprobs);
    isdiploid = false;
  }
  else{
    chrm->UpdateParameters( this, thetaX, options, rhoX, chibindicator, true, calcbackprobs);
    isdiploid = true;
  }
  return isdiploid;
}

void Individual::SampleRho(bool XOnly, bool RandomMatingModel, bool X_data, double rhoalpha, double rhobeta, 
			   unsigned int SumN[], unsigned int SumN_X[]){
  double L = Loci->GetLengthOfGenome(), L_X=0.0;
  if( Loci->isX_data() ) L_X = Loci->GetLengthOfXchrm();

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

  if( options->getTestForAffectedsOnly() ){
    for(unsigned j = 0; j < Loci->GetNumberOfCompositeLoci()*KK; ++j){
      AffectedsScore[j] = 0.0;
      AffectedsVarScore[j] = 0.0;
      AffectedsInfo[j] = 0.0;
    }
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

//********************** Score Tests ***************************

void Individual::UpdateScoreForLinkageAffectedsOnly(int j, bool RandomMatingModel, Chromosome **chrm){
  // Different from the notation in McKeigue et  al. (2000). McKeigue
  // uses P0, P1, P2, which relate to individual admixture as follows;
  // P0 = ( 1 - theta0 ) * ( 1 - theta1 )
  // P1 = ( 1 - theta0 ) * theta1 + theta0 * ( 1 - theta1 )
  // P2 = theta0 * theta1

  // Test in McKeigue et al. (2000) was on r defined on the positive
  // real line.  This test is on log r, which should have better
  // asymptotic properties.

  double r1 = 0.5;
  double r2 = 2.0;//hard-coding these for now, can make them vary later

  //we don't bother computing scores for the first population when there are two
  int KK = Populations,k1 = 0;
  if(Populations ==2) {KK = 1;k1 = 1;}
  
  double theta[2];//paternal and maternal admixture proportions
  double AProbs[Populations][3];

  double Pi[3];//probs of 0,1,2 copies of Pop1 given admixture
  int offset = 0;
  if(!RandomMatingModel)offset = Populations;


  int locus;
  for( unsigned int jj = 0; jj < chrm[j]->GetSize(); jj++ ){
    locus = chrm[j]->GetLocus(jj); 
    //retrieve AncestryProbs from HMM
    chrm[j]->getAncestryProbs( jj, AProbs );

    for( int k = 0; k < KK; k++ ){
      theta[0] = Theta[ k+k1 ];
      if( RandomMatingModel )
	theta[1] = Theta[ Populations + k+k1 ];
      else
	theta[1] = theta[0];
      
      //accumulate score, score variance, and info
      AffectedsScore[locus *KK + k]+= 0.5*( AProbs[k+k1][1] + 2.0*AProbs[k+k1][2] - theta[0] - theta[1] );
      AffectedsVarScore[locus * KK + k]+= 0.25 *( AProbs[k+k1][1]*(1.0 - AProbs[k+k1][1]) + 4.0*AProbs[k+k1][2]*AProbs[k+k1][0]); 
      AffectedsInfo[locus * KK +k]+= 0.25* ( theta[0]*( 1.0 - theta[0] ) + theta[1]*( 1.0 - theta[1] ) );

      //probs of 0,1,2 copies of Pop1 given admixture
      Pi[2] = theta[0] * theta[1];
      Pi[1] = theta[0] * (1.0 - theta[1]);
      Pi[0] = (1.0 - theta[0]) * (1.0 - theta[1]);

      //compute contribution to likelihood ratio
      LikRatio1[locus *KK + k] += (AProbs[k+k1][0] + sqrt(r1)*AProbs[k+k1][1] + r1 * AProbs[k+k1][2]) / 
	(Pi[0] + sqrt(r1)*Pi[1] + r1*Pi[2]);
      LikRatio2[locus *KK + k] += (AProbs[k+k1][0] + sqrt(r2)*AProbs[k+k1][1] + r2 * AProbs[k+k1][2]) / 
	(Pi[0] + sqrt(r2)*Pi[1] + r2*Pi[2]);
    }
    
    ++locus;
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
  //Xcopy is an exact copy of X; We need two copies as one will be destroyed
  double xBx[1], BX[Populations];
  double VarA[Populations];
 
  X[ 2*Populations - 1] = 1;//intercept
  Xcov[Populations-1] = 1;
  //set covariates, admixture props for pops 2 to K 
  for( int k = 0; k < Populations - 1; k++ ){
    X[ Populations + k] = Theta[ k+1 ];
    BX[k] = Xcov[k] = Theta[ k+1 ];
  }

  int locus; 
  for( unsigned int jj = 0; jj < chrm[j]->GetSize(); jj++ ){
    locus = chrm[j]->GetLocus(jj);      
    chrm[j]->getAncestryProbs( jj, AProbs );//conditional locus ancestry probs      
    
    for( int k = 0; k < Populations ; k++ ){
      Xcopy[k] = X[k] = AProbs[k][1] + 2.0 * AProbs[k][2];//Conditional expectation of ancestry
      VarA[k] = AProbs[k][1]*(1.0 - AProbs[k][1]) + 4.0*AProbs[k][2]*AProbs[k][0];//conditional variances
      }
    //KLUDGE: need to reset Xcopy each time since destroyed in computation of score
    Xcopy[2*Populations-1] = 1;
    for( int k = 0; k < Populations-1; k++ )Xcopy[k + Populations] = Theta[ k+1 ];

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

void Individual::UpdateB(double DInvLink, double dispersion){
  //increment B using new Admixture Props
  //Xcov is a vector of admixture props as covariates as in UpdateScoreForAncestry
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

//******************** Chib Algorithm ***************************************

void Individual::InitializeChib(double *theta, double *thetaX, vector<double> rho, vector<double> rhoX, 
				AdmixOptions *options, AlleleFreqs *A, Chromosome **chrm, double rhoalpha, double rhobeta, 
				vector<vector<double> > &alpha, chib *MargLikelihood, LogWriter *Log)
//Computes LogPrior and LogLikelihood used for Chib Algorithm
  //called once after burnin to set estimates
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
      LogPrior += getDirichletLogDensity( alpha[0], theta );
   }
   // ** case of some x data **
   else if( Loci->isX_data() ){
     _rhoHat = rho;
     _rhoHat_X = rhoX;
     copy(theta, theta + theta_size, ThetaHat); //ThetaHat = theta
     copy(theta, theta + theta_size, ThetaXHat); //ThetaXHat = theta

     LogLikelihoodAtEst = getLogLikelihood( options, chrm, theta, thetaX, rho, rhoX, true);
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
      LogPrior += getDirichletLogDensity( alpha[0], theta )
               + getDirichletLogDensity( alpha[0], thetaX )
               + getDirichletLogDensity( alpha[1], theta+K )
               + getDirichletLogDensity( alpha[1], thetaX + K );

      _rhoHat = rho;
      _rhoHat_X = rhoX;
      copy(theta, theta + theta_size, ThetaHat); 
      copy(theta, theta + theta_size, ThetaXHat);
      LogLikelihoodAtEst = getLogLikelihood( options, chrm, theta, thetaX, rho, rhoX, true);
   }
   else{
      if( Populations > 1 ){
 	_rhoHat = rho;
	copy(theta, theta + theta_size, ThetaHat); 
	LogLikelihoodAtEst = getLogLikelihood( options, chrm, theta, thetaX, rho, rhoX, true);
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
            LogPrior += getDirichletLogDensity( alpha[0], theta );
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
            LogPrior += getDirichletLogDensity( alpha[1], theta + K );
         }
      }
      else{
	LogLikelihoodAtEst = getLogLikelihoodOnePop(true);
      }
   }
   if( A->IsRandom() ){
      for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
         for( int k = 0; k < Populations; k++ ){
	   LogPrior += getDirichletLogDensity( A->GetPriorAlleleFreqs(j, k), A->getAlleleFreqsMAP(j,k) );
         }
      }
   }
   MargLikelihood->setLogPrior( LogPrior );
   MargLikelihood->setLogLikelihood( LogLikelihoodAtEst );   
}

 // this function computes marginal likelihood by the Chib algorithm.  
void Individual::ChibLikelihood(int iteration, double *LogLikelihood, double *SumLogLikelihood, double *MaxLogLikelihood,
				AdmixOptions *options, Chromosome **chrm, vector<vector<double> > &alpha, double globalrho,
				double rhoalpha, double rhobeta, double *thetahat,
				double *thetahatX, vector<double> &rhohat,
				vector<double> &rhohatX, LogWriter *Log, chib *MargLikelihood, AlleleFreqs* A){
vector<double> rho(2);
if(options->getRhoIndicator())
  rho = _rho;
else // global rho
  rho[0] = rho[1] = globalrho;
            
  //if( iteration <= options->getBurnIn() ){
  int K = Populations;
  size_t theta_size = Populations;
  if(options->isRandomMatingModel()) theta_size *=2;
  
  *LogLikelihood = getLogLikelihood(options, chrm);//gets loglikelihood at current parameter values

    if( K > 1 ){
      if( !options->RhoFlatPrior() && !options->logRhoFlatPrior() ){
	if( options->isAdmixed(0) ){
	  *LogLikelihood+=getGammaLogDensity( rhoalpha, rhobeta, rho[0] );}
	if(  options->isAdmixed(1) )
	  *LogLikelihood+=getGammaLogDensity( rhoalpha, rhobeta, rho[1] );
      }
      else if( options->logRhoFlatPrior() ){
	if( options->isAdmixed(0) )
	  *LogLikelihood -= log( rho[0] );
	if(  options->isAdmixed(1) )
	  *LogLikelihood -= log( rho[1] );
      }
      *LogLikelihood+= getDirichletLogDensity(alpha[0], Theta) + getDirichletLogDensity(alpha[1], Theta + K);

      if( *LogLikelihood > *MaxLogLikelihood ){
	Log->write("Admixture (gamete 1):");
	Log->write(Theta, K);Log->write("\n");
	Log->write("Admixture (gamete 2):");
	Log->write(Theta+K, K);Log->write("\n");
	Log->write("sumintensities: ");
	Log->write( rho[0]);Log->write( rho[1]);
	Log->write("\nLogLikelihood:");Log->write( *LogLikelihood);
	Log->write("\niteration: ");Log->write(iteration);Log->write("\n\n");

	*MaxLogLikelihood = *LogLikelihood;

	//during burnin
	if( iteration <= options->getBurnIn() ){
	  //set allelefreqsMAP to current values of allelefreqs
	  A->setAlleleFreqsMAP();
	  //set HapPairProbsMAP to current values of HapPairProbs
	  for( unsigned  j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
	    //if( locus->GetNumberOfLoci() > 2 )
	      (*Loci)(j)->setHaplotypeProbsMAP();
	  }
	  //
	  for(unsigned k = 0; k < theta_size; ++k)thetahat[k] = Theta[k];
	  rhohat = rho;

	  if( Loci->isX_data() ){
	  for(unsigned k = 0; k < theta_size; ++k)thetahatX[k] = ThetaX[k];
	    rhohatX = _rho_X;
	  }
	}//end during burnin block
      }
    }//end if K>1

    else{//single population
      if( *LogLikelihood > *MaxLogLikelihood ){

	*MaxLogLikelihood = *LogLikelihood;
	if( iteration <= options->getBurnIn() ){
	  A->setAlleleFreqsMAP();
	  for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
	    if( (*Loci)(j)->GetNumberOfLoci() > 2 )
	      (*Loci)(j)->setHaplotypeProbsMAP();
	  }
	}
      }
    }
    //if( options->getAnalysisTypeIndicator() == -1 ){
      if( iteration == options->getBurnIn() ){
	InitializeChib(thetahat, thetahatX, rhohat, rhohatX, options, A, chrm, rhoalpha, rhobeta, alpha, MargLikelihood, Log);
      }
      if( iteration > options->getBurnIn() ){
	double Log_Posterior = 0;
	if( Populations > 1 )
	  Log_Posterior = LogPosterior;//LogPOsterior as calculate by CalculateLogPosterior
	if( A->IsRandom() ){
	  for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
	    for( int k = 0; k < Populations; k++ ){
	      Log_Posterior += getDirichletLogDensity( A->GetPriorAlleleFreqs(j,k) + A->GetAlleleCounts(j,k), 
						      A->getAlleleFreqsMAP(j, k) );
	    }
	  }
	}
	MargLikelihood->addLogPosteriorObs( Log_Posterior );
	*SumLogLikelihood += *LogLikelihood;
     }
      //}

}

//called on first infividual after burnin
void Individual::CalculateLogPosterior(AdmixOptions *options, vector<vector<double> > &alpha, 
				       double rhoalpha, double rhobeta){
  double L = Loci->GetLengthOfGenome(), L_X=0.0;
  if( Loci->isX_data() ) L_X = Loci->GetLengthOfXchrm();
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
      transform(alpha[0].begin(), alpha[0].end(), SumLocusAncestry_X, alphaparams0.begin(), std::plus<double>());
      LogPosterior += getDirichletLogDensity(alphaparams0, ThetaHat);
    }
    else if( Loci->isX_data() ){
      for( unsigned int g = 0; g < 2; g++ ){
	LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN[g], rhobeta + L, _rhoHat[g] );
	if(!options->RhoFlatPrior() && !options->logRhoFlatPrior() )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN[g], rhobeta+L, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[g], 1.0);
	LogPosterior -= log(IntConst1);
	transform(alpha[g].begin(), alpha[g].end(), SumLocusAncestry +g*Populations, alphaparams0.begin(), std::plus<double>());
	LogPosterior += getDirichletLogDensity(alphaparams0, ThetaHat + g*Populations);
      }
      for( unsigned int g = 0; g < gametes[X_posn]; g++ ){
	LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN_X[g], rhobeta + L_X, _rhoHat_X[g] );
	if( options->RhoFlatPrior() || options->logRhoFlatPrior() )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN_X[g], rhobeta+L_X, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L_X, rhoalpha+(double)SumN_X[g], 1.0);
	LogPosterior -= log(IntConst1);
	transform(alpha[g].begin(), alpha[g].end(), SumLocusAncestry_X +g*Populations, alphaparams0.begin(), std::plus<double>());
	LogPosterior += getDirichletLogDensity(alphaparams0, ThetaXHat+g*Populations);
      }
    }
    else if( options->isSymmetric() ){//both gametes have same admixture
      vector<double> x(2,0.0);
      x[0] += getGammaLogDensity( rhoalpha + (double)SumN[0], rhobeta + L, _rhoHat[0] );
      x[1] += getGammaLogDensity( rhoalpha + (double)SumN[1], rhobeta + L, _rhoHat[0] );
      x[0] += getGammaLogDensity( rhoalpha + (double)SumN[1], rhobeta + L, _rhoHat[1] );
      x[1] += getGammaLogDensity( rhoalpha + (double)SumN[0], rhobeta + L, _rhoHat[1] );
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
      transform(alpha[0].begin(), alpha[0].end(), SumLocusAncestry, alphaparams0.begin(), std::plus<double>());
      x[0] += getDirichletLogDensity(alphaparams0, ThetaHat);
      x[1] += getDirichletLogDensity(alphaparams0, ThetaHat+Populations);
      transform(alpha[1].begin(), alpha[1].end(), SumLocusAncestry+Populations, alphaparams1.begin(), std::plus<double>());
      x[0] += getDirichletLogDensity(alphaparams1, ThetaHat+Populations);
      x[1] += getDirichletLogDensity(alphaparams1, ThetaHat);
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
    else{//different admixtures for each gamete
      if(  options->isAdmixed(0) ){//admixed first gamete
	LogPosterior = getGammaLogDensity( rhoalpha + (double)SumN[0],
					   rhobeta + L, _rhoHat[0] );
	if( options->RhoFlatPrior() || options->logRhoFlatPrior() )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN[0], rhobeta+L, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[0], 1.0);
           LogPosterior -= log( IntConst1 );
	   transform(alpha[0].begin(), alpha[0].end(), SumLocusAncestry, alphaparams0.begin(), std::plus<double>());
	   LogPosterior += getDirichletLogDensity(alphaparams0, ThetaHat);
      }
      if(  options->isAdmixed(1) ){//admixed second gamete
	LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN[1],
					    rhobeta + L, _rhoHat[1] );
	if( options->RhoFlatPrior() || options->logRhoFlatPrior() )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN[1], rhobeta+L, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[1], 1.0);
	LogPosterior -= log( IntConst1 );
	transform(alpha[1].begin(), alpha[1].end(), SumLocusAncestry+Populations, alphaparams1.begin(), std::plus<double>());
	LogPosterior += getDirichletLogDensity(alphaparams1, ThetaHat+Populations);
      }
    }
  }
}

double Individual::IntegratingConst( double alpha, double beta, double a, double b )
{
   double I = gsl_cdf_gamma_P( b*beta, alpha, 1 ) - gsl_cdf_gamma_P( a*beta, alpha, 1);
   return I;
}

//******************* Likelihood Ratios (Affecteds-only score test) ************************

static string double2R( double x )
{
  if( isnan(x) )
    return "NaN";
  else{
    stringstream ret;
    ret << x;
    return( ret.str() );
  }
}

void Individual::OutputLikRatios(const char *filename, int iterations, std::string *PopLabels){
  //open outut file
  std::ofstream outputstream(filename);

  //start writing R object
  outputstream<< "structure(.Data=c(" << endl;

  //output ergodic averages of LikRatios
  int KK = Populations, k1 = 0;
  if(KK == 2 ){
    KK = 1;k1 = 1;
  }
  double L1, L2;

  for(unsigned int j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
    for( int k = 0; k < KK; k++ ){//end at 1 for 2pops
      outputstream << (*Loci)(j)->GetLabel(0) << ",";
      outputstream << "\""<<PopLabels[k+k1] << "\","; //need offset to get second poplabel for 2pops
      
      L1 = LikRatio1[ j*KK + k] / ( iterations );
      L2 = LikRatio2[ j*KK + k] / ( iterations );
      
      outputstream << double2R(L1)<< ","
		   << double2R(L2)<< ","<<endl;
    }
  }
  
  vector<int> dim(2,0);
  dim[0] = 4;
  dim[1] = Loci->GetNumberOfCompositeLoci() * KK;
  
  vector<string> labels(4,"");
  labels[0] = "Locus";
  labels[1] = "Population";
  labels[2] = "L1";
  labels[3] = "L2";
  
  outputstream << ")," << endl;
  outputstream << ".Dim = c(";
  for(unsigned int i=0;i<dim.size();i++){
    outputstream << dim[i];
    if(i != dim.size() - 1){
      outputstream << ",";
    }
  }
  outputstream << ")," << endl;
  outputstream << ".Dimnames=list(c(";
  for(unsigned int i=0;i<labels.size();i++){
    outputstream << "\"" << labels[i] << "\"";
    if(i != labels.size() - 1){
      outputstream << ",";
    }
  }
  outputstream << "), character(0)))" << endl;

  outputstream.close();

}
