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
const Genome *Individual::Loci;
int Individual::Populations;

//******** Constructors **********
Individual::Individual()
{//should initialise pointers here
}

Individual::Individual(int number, const AdmixOptions* const options, const InputData* const Data, const Genome& Loci, 
		       const Chromosome* const * chrm)
{
  myNumber = number;

  if( !options->isGlobalRho() ){
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
  
  // Read sex value if present.
  sex = male;
  if (options->getgenotypesSexColumn() == 1) {
    sex = Data->GetSexValue(myNumber);
  }
  
  int numCompositeLoci = Loci.GetNumberOfCompositeLoci();
  
  LocusAncestry = new int*[ numChromosomes ]; // array of matrices in which each col stores 2 integers 
  
  Theta = 0;
  ThetaX = 0;
  ThetaProposal = 0;
  ThetaXProposal = 0;
  SumSoftmaxTheta = 0;
  
  // SumLocusAncestry is sum of locus ancestry states over loci at which jump indicator xi is 1  
  SumLocusAncestry = new int[options->getPopulations()*2];
  
  
  if( options->isRandomMatingModel() ){//random mating model
    ThetaProposal = new double[ Populations * 2 ];
    Theta = new double[ Populations * 2 ];
    SumSoftmaxTheta = new double[ Populations * 2 ];
    fill(SumSoftmaxTheta, SumSoftmaxTheta + Populations*2, 0.0);
  }
  else{
    ThetaProposal = new double[ Populations];
    Theta = new double[ Populations ];
    SumSoftmaxTheta = new double[ Populations ];
    fill(SumSoftmaxTheta, SumSoftmaxTheta + Populations, 0.0);
  }
  
  //X chromosome objects
  SumLocusAncestry_X = 0;    
  //    if(Loci.isX_data() ){
  if( sex == male ){
    ThetaXProposal = new double[ Populations];
    ThetaX = new double[Populations];
    //ThetaXHat = new double[Populations];
    SumLocusAncestry_X = new int[Populations];
  }
  else{
    ThetaXProposal = new double[ Populations * 2 ];
    ThetaX = new double[ Populations * 2 ];
    //ThetaXHat = new double[ Populations * 2 ];
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
  Data->GetGenotype(myNumber, options->getgenotypesSexColumn(), Loci, &genotypes);
  
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

  logLikelihood.value = 0.0;
  logLikelihood.ready = false;
  logLikelihood.HMMisOK = false;   
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
void Individual::SetStaticMembers(const Genome* const pLoci, const AdmixOptions* const options){
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
void Individual::setAdmixtureProps(const double* const a, size_t size)
{
  for(unsigned i = 0; i < size; ++i)  Theta[i] = a[i];
}

void Individual::setAdmixturePropsX(const double* const a, size_t size)
{
  for(unsigned i = 0; i < size; ++i)  ThetaX[i] = a[i];
}
void Individual::HMMIsBad(bool loglikisbad){
  logLikelihood.HMMisOK = false;
  if(loglikisbad)logLikelihood.ready = false;
}
//******************** Accessors ***********************************************************
const unsigned short* const* Individual::getGenotype(unsigned int locus)const{
  return genotypes[locus];
}

const std::vector<hapPair > &Individual::getPossibleHapPairs(unsigned int locus)const{
  return PossibleHapPairs[locus];
}

const double* Individual::getAdmixtureProps()const
{
  return Theta;
}

const Sex Individual::getSex()const 
{
   return sex;
}


double Individual::getSumrho()const
//returns sum of sumintensities over gametes
{
   double sumrho = 0;
   for( unsigned int g = 0; g < _rho.size(); g++ )
      sumrho += _rho[g];
   return sumrho;
}

const vector<double> Individual::getRho()const
{
   return _rho;
}

void Individual::GetLocusAncestry(int chrm, int locus, int Ancestry[2])const{
  Ancestry[0]  = LocusAncestry[chrm][locus];
  if((unsigned)chrm == X_posn)Ancestry[1] = Ancestry[0];
  else Ancestry[1] = LocusAncestry[chrm][Loci->GetSizesOfChromosomes()[chrm]  + locus];
}

//returns value of LocusAncestry at a locus for a particular gamete
int Individual::GetLocusAncestry(int chrm, int gamete, int locus)const{
  int g = (gametes[chrm] == 2) ? gamete : 0; //so that gamete = 1 works when gametes[chrm] = 1;
  return LocusAncestry[chrm][g * Loci->GetSizesOfChromosomes()[chrm]  + locus] ;
}

const int *Individual::getSumLocusAncestry()const{
  return SumLocusAncestry;
}

//********** Missing Genotype Indicator *********************

//Indicates whether genotype is missing at all simple loci within a composite locus
bool Individual::IsMissing(unsigned int locus)const
{
  unsigned int count = 0;
  int NumberOfLoci = (*Loci)(locus)->GetNumberOfLoci();
  for(int i = 0; i < NumberOfLoci; i++){
    count += genotypes[locus][i][0];
  }
  return (count == 0);
}

//****************** Log-Likelihoods **********************
double Individual::getLogLikelihood( const AdmixOptions* const options, Chromosome **chrm){
  //use current parameter values
  //only call with annealindicator = true once per individual per iteration to accumulate unannealed loglikelihood
  double LogLikelihood = 0.0;
  double Probs[Populations*Populations];
  int locus = 0;
  if(!logLikelihood.ready){
      for( unsigned int j = 0; j < numChromosomes; j++ ){
	if(Populations == 1){
	  for(unsigned jj = 0; jj < chrm[j]->GetNumberOfCompositeLoci(); ++jj){
	    (*Loci)(locus)->GetGenotypeProbs(Probs, getPossibleHapPairs(locus), false);
	    ++locus;
	  }
	  LogLikelihood += log( Probs[0] );
	}
	else{
	  //need to run forward recursions in HMM if current forward probs do not correspond to current parameter values
	  if(!logLikelihood.HMMisOK){
	    chrm[j]->SetGenotypeProbs(this, false);
	    UpdateHMMForwardProbs(j, chrm[j], options, Theta, ThetaX, _rho, _rho_X);
	  }
	  LogLikelihood += chrm[j]->getLogLikelihood();
	}
      }
    logLikelihood.value = LogLikelihood;
  }
  logLikelihood.ready = true;
  logLikelihood.HMMisOK = true;//because forward probs now correspond to current parameter values 
  //and UpdateHMMForwardProbs sets this to false
  return logLikelihood.value;
}

double Individual::getLogLikelihoodAtPosteriorMeans(const AdmixOptions* const options, Chromosome **chrm){
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

  double LogLikelihood = 0.0;
  double Probs[Populations*Populations];
  int locus = 0;
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    if(Populations == 1){
      for(unsigned jj = 0; jj < chrm[j]->GetNumberOfCompositeLoci(); ++jj){
	(*Loci)(j)->GetGenotypeProbs(Probs, getPossibleHapPairs(j),false);
	++locus;
      }
      LogLikelihood += log( Probs[0] );
    }
    else{
      chrm[j]->SetGenotypeProbs(this, false);//will set genotype probs to posterior means as happairprobs are already so
      UpdateHMMForwardProbs(j, chrm[j], options, ThetaBar, ThetaBar, sumlogrho, sumlogrho);
      LogLikelihood += chrm[j]->getLogLikelihood();
    }
  }
  return LogLikelihood;
}
double Individual::getLogLikelihood( const AdmixOptions* const options, Chromosome **chrm, const double* const theta, 
				     const double* const thetaX,
				     const vector<double > rho, const vector<double> rho_X, bool chibindicator = false)
//updates forward probs in HMM and retrieves loglikelihood at supplied theta and rho
//(optional)chibindicator = true for computing LogL at MLEs; instructs CompositeLocus to use HapPairProbsMAP
//instead of HapPairProbs, when allelefreqs are not fixed, in calculating GenotypeProbs.
{
   double LogLikelihood = 0.0;
   double Probs[Populations*Populations];
   int locus = 0;
   for( unsigned int j = 0; j < numChromosomes; j++ ){
     if(Populations == 1){ // likelihood is calculated as product of genotype probs at all loci
       for(unsigned jj = 0; jj < chrm[j]->GetNumberOfCompositeLoci(); ++jj){
	 (*Loci)(locus)->GetGenotypeProbs(Probs, getPossibleHapPairs(locus), chibindicator);
	 ++locus;
       }
       LogLikelihood += log( Probs[0] );
     }
     else{ // likelihood calculated from HMM
       chrm[j]->SetGenotypeProbs(this, chibindicator); 
       UpdateHMMForwardProbs(j, chrm[j], options, theta, thetaX, rho, rho_X);
       LogLikelihood += chrm[j]->getLogLikelihood();
     }
   }

//    double LogLikelihood = 0.0;
//    if(Populations == 1)LogLikelihood = getLogLikelihoodOnePop(chibindicator);//in case this version called when only one population
//    else{
//      for( unsigned int j = 0; j < numChromosomes; j++ ){
//        chrm[j]->SetGenotypeProbs(this, chibindicator);      
//        UpdateHMMForwardProbs(j, chrm[j], options, theta, thetaX, rho, rho_X, false);
//        LogLikelihood += chrm[j]->getLogLikelihood();
//      }
//    }
    return LogLikelihood;
}
double Individual::getLogLikelihoodOnePop(){ //convenient for a single population as no arguments required
  double logLikelihood = 0.0;
  double *Prob;
  Prob = new double[1];//one pop so 1x1 array
  for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
    if(!IsMissing(j)){
      (*Loci)(j)->GetGenotypeProbs(Prob,getPossibleHapPairs(j), false);
      logLikelihood += log( Prob[0] );
    }
  }
  return logLikelihood;
}
//************** Updating (Public) ***************************************************************************************

// unnecessary duplication of code - ? should embed within method for > 1 population
void Individual::OnePopulationUpdate( int i, DataMatrix *Outcome, int NumOutcomes, const DataType* const OutcomeType, 
				      const double* const* ExpectedY, const double* const lambda, const Chromosome* const*chrm, 
				      AlleleFreqs *A )
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

void Individual::SampleParameters( int i, double *SumLogTheta, AlleleFreqs *A, int iteration , DataMatrix *Outcome,
			 int NumOutcomes, const DataType* const OutcomeType, const double* const * ExpectedY, 
			 const double* const lambda, int NoCovariates,
			 DataMatrix *Covariates, const double* const* beta, const double *poptheta, const AdmixOptions* const options,
			 Chromosome **chrm, const vector<vector<double> > &alpha,  
			 double rhoalpha, double rhobeta, const vector<double> sigma, 
				   double DInvLink, double dispersion, bool anneal=false)
/*arguments:
  i = this individuals's number-1
  SumLogTheta = array in IndividualCollection holding sums of log admixture props
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
		Covariates, beta, poptheta, options, alpha, sigma, DInvLink, dispersion, true, anneal);
  
  //SumN is the number of arrivals between each pair of adjacent loci
  SumN[0] = SumN[1] = 0;
  SumN_X[0] = SumN_X[1] = 0;  
  
  bool calcbackprobs = (options->getTestForAffectedsOnly() || options->getTestForLinkageWithAncestry());
  for( unsigned int j = 0; j < numChromosomes; j++ ){
    
    //Update Forward/Backward probs in HMM
    if( !logLikelihood.HMMisOK ){
      chrm[j]->SetGenotypeProbs(this, false);
      UpdateHMMForwardProbs(j, chrm[j], options, Theta, ThetaX, _rho, _rho_X);
    }
    if(calcbackprobs)chrm[j]->UpdateHMMBackwardProbs(Theta);//TODO: pass correct theta for haploid case
    
    //update score tests for linkage with ancestry for *previous* iteration
    if(iteration > options->getBurnIn()){
      //Update affecteds only scores
      if( options->getTestForAffectedsOnly()){
	//determine which regression is logistic, in case of 2 outcomes
	unsigned col = 0;
	if(options->getNumberOfOutcomes() >1 && OutcomeType[0]!=Binary)col = 1;
	if(options->getNumberOfOutcomes() == 0 || Outcome->get(i, col) == 1)
	  UpdateScoreForLinkageAffectedsOnly(j, options->isRandomMatingModel(),chrm );
      }

      //update ancestry score tests
      if( options->getTestForLinkageWithAncestry() ){   
	UpdateScoreForAncestry(j, dispersion, Outcome->get(i, 0) - ExpectedY[0][i], DInvLink,chrm);
      }    
    }

    //Sample locus ancestry
    // sampling locus ancestry requires calculation of forward probability vectors alpha in HMM 
    chrm[j]->SampleLocusAncestry(LocusAncestry[j], Theta);
 
   //loop over loci on current chromosome and update allele counts
    for( unsigned int jj = 0; jj < chrm[j]->GetSize(); jj++ ){
      int locus =  chrm[j]->GetLocus(jj);
      if( !(IsMissing(j)) ){
	  int anc[2];//to store ancestry states
	  GetLocusAncestry(j,jj,anc);
	  int h[2];//to store sampled hap pair
	  //might be a shortcut for haploid data since there is only one compatible hap pair, no need to sample
	  (*Loci)(locus)->SampleHapPair(h, PossibleHapPairs[locus], anc);
	  A->UpdateAlleleCounts(locus,h, anc, chrm[j]->isDiploid());
      }
     }   

    //sample number of arrivals and SumLocusAncestry
    bool isX = (j == X_posn);
    chrm[j]->SampleJumpIndicators(LocusAncestry[j], gametes[j], SumLocusAncestry, SumLocusAncestry_X, isX,
				  SumN, SumN_X, options->isGlobalRho());

  }//end chromosome loop

  // sample sum of intensities parameter rho
  if( !options->isGlobalRho() ){
     SampleRho( options->getXOnlyAnalysis(), options->isRandomMatingModel(), Loci->isX_data(), rhoalpha, rhobeta, 
	       SumN, SumN_X);
  }
  if(!anneal && iteration > options->getBurnIn()){
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

void Individual::SampleTheta( int i, int iteration, double *SumLogTheta, const DataMatrix* const Outcome, Chromosome ** C,
			      int NumOutcomes, const DataType* const OutcomeType, const double* const* ExpectedY, 
			      const double* const lambda, int NoCovariates,
			      DataMatrix *Covariates, const double* const* beta, const double* const poptheta,
			      const AdmixOptions* const options, const vector<vector<double> > &alpha, const vector<double> sigma,
			      double DInvLink, double dispersion, bool RW, bool anneal=false)
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

  if(!anneal && iteration > options->getBurnIn()){
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

double Individual::ProposeThetaWithRandomWalk(const AdmixOptions* const options, Chromosome **C, const vector<vector<double> > &alpha){
  //TODO: X-chromosome case
  double LogLikelihoodRatio = 0.0;
  double LogPriorRatio = 0.0;
  
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
    
    a[Populations-1] = 0.0;
    //random walk step
    for(int k = 0; k < last; ++k)if(Theta[g*Populations + k]>0.0){
      a[k] = gennor(a[k], step0);
      a[last] -= a[k];
    }
    
    //reverse transformation
    softmax(Populations, ThetaProposal+g*Populations, a, b);
    
    //compute prior ratio
    LogPriorRatio += getDirichletLogDensity_Softmax(alpha[0], ThetaProposal+g*Populations) - 
      getDirichletLogDensity_Softmax(alpha[0], Theta+g*Populations);
  }
  //get log likelihood at current parameter values
  LogLikelihoodRatio -= getLogLikelihood(options, C);
  
  //get log likelihood at proposal theta and current rho
  LogLikelihoodRatio += getLogLikelihood(options, C, ThetaProposal, ThetaXProposal,_rho, _rho_X, false);

  return LogLikelihoodRatio + LogPriorRatio;// = log Posterior ratio
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
void Individual::ProposeTheta(const AdmixOptions* const options, const vector<double> sigma, const vector<vector<double> > &alpha){
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
double Individual::LogAcceptanceRatioForRegressionModel( int i, RegressionType RegType, int TI,  bool RandomMatingModel, 
							 int Populations, int NoCovariates, 
							 const DataMatrix* const Covariates, const double* const* beta, 
							 const double* const* ExpectedY, const DataMatrix* const Outcome, 
							 const double* const poptheta, const double* const lambda)
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
double Individual::AcceptanceProbForTheta_XChrm(const std::vector<double> &sigma, int Populations )
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
                                               const double* const poptheta, bool RandomMatingModel, 
					       DataMatrix *Covariates)
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
	if(RW)logLikelihood.HMMisOK = true;
	//with a random-walk update, the current values in HMM now correspond to the proposal,
	//which have now become the current values
	//next time the loglikelihood at current parameter values is required, it can be pulled directly from HMM
	else {
	  //with a conjugate update the values in HMM currently correspond to old parameter values so the HMMs will need
	  //to be updated again to get loglikelihood
	  logLikelihood.HMMisOK = false;
	}
	  logLikelihood.ready = false;
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

//Updates forward probabilities in HMM for chromosome j
//also sets Diploid flag in Chromosome (last arg of UpdateParameters)
void Individual::UpdateHMMForwardProbs(unsigned int j, Chromosome* const chrm, const AdmixOptions* const options, 
				       const double* const theta, const double* const thetaX, 
				       const vector<double> rho, const vector<double> rhoX){

  if( j != X_posn ){// NOT an X chromosome
    chrm->UpdateHMMForwardProbs(theta, options, rho, true);
  }
  else if( options->getXOnlyAnalysis() ){//X only data
    chrm->UpdateHMMForwardProbs(theta, options, rho, false);
  }
  else if( sex == male ){// X chromosome in male individual
    chrm->UpdateHMMForwardProbs(thetaX, options, rhoX, false);
  }
  else{// X chromosome in female individual or individual whose sex is unknown
    chrm->UpdateHMMForwardProbs(thetaX, options, rhoX, true);
  }
  logLikelihood.HMMisOK = false;//because forward probs in HMM have been changed
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
  else{//assortative mating
    _rho[0] = gengam( rhobeta + 2*L, rhoalpha + (double)(SumN[0] + SumN[1]) );
  }
  //now that rho has changed, current stored value of loglikelihood is no longer valid and 
  //HMMs will need to be updated before getting loglikelihood
  logLikelihood.HMMisOK = false;
  logLikelihood.ready = false;
}

void Individual::ResetScores(const AdmixOptions* const options){
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

void Individual::UpdateScoreForLinkageAffectedsOnly(int j, bool RandomMatingModel, const Chromosome* const* chrm){
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

void Individual::UpdateScoreForAncestry(int j,double phi, double YMinusEY, double DInvLink, const Chromosome* const*chrm)
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


 // this function computes marginal likelihood by the Chib algorithm.  
void Individual::Chib(int iteration, double *SumLogLikelihood, double *MaxLogLikelihood,
		      const AdmixOptions* const options, Chromosome **chrm, const vector<vector<double> > &alpha, 
		      double globalrho, double rhoalpha, double rhobeta, double *thetahat,
		      double *thetahatX, vector<double> &rhohat,
		      vector<double> &rhohatX, LogWriter* Log, chib *MargLikelihood, AlleleFreqs* A){
  vector<double> rho(2);
  if(!options->isGlobalRho())
    rho = _rho;
  else // global rho
    rho[0] = rho[1] = globalrho;

  int K = Populations;
  size_t theta_size = Populations;
  if(options->isRandomMatingModel()) theta_size *=2;
  double logLikelihood = 0.0;

  // *** Every iteration ***
  
  logLikelihood = getLogLikelihood(options, chrm);//gets loglikelihood at current parameter values
  
  // *** during BurnIn ***
  if( iteration <= options->getBurnIn() ){
    if( Populations > 1 ){  
      if( logLikelihood > *MaxLogLikelihood ){
	
	  Log->write("Admixture (gamete 1):");
	  Log->write(Theta, K);Log->write("\n");
	  Log->write("Admixture (gamete 2):");
	  Log->write(Theta+K, K);Log->write("\n");
	  Log->write("sumintensities: ");
	  Log->write( rho[0]);Log->write( rho[1]);
	  Log->write("\nLogLikelihood:");Log->write( logLikelihood);
	  Log->write("\niteration: ");Log->write(iteration);Log->write("\n\n");
	  
	  //set parameter estimates at max loglikelihood
	  for(unsigned k = 0; k < theta_size; ++k)thetahat[k] = Theta[k];
	  rhohat = rho;
	  
	  if( Loci->isX_data() ){
	    for(unsigned k = 0; k < theta_size; ++k)thetahatX[k] = ThetaX[k];
	    rhohatX = _rho_X;
	  }
      }//end if Loglikelihood > Max
    }//end if K>1
  
  
    if( logLikelihood > *MaxLogLikelihood ){
      *MaxLogLikelihood = logLikelihood;
      
	//set allelefreqsMAP to current values of allelefreqs
	A->setAlleleFreqsMAP();
	//set HapPairProbsMAP to current values of HapPairProbs
	for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
	  //if( (*Loci)(j)->GetNumberOfLoci() > 2 )
	  (*Loci)(j)->setHaplotypeProbsMAP();
	}
    }
  }

  // *** At end of BurnIn ***
  if( iteration == options->getBurnIn() ){
    
    //loglikelihood at estimates
    MargLikelihood->setLogLikelihood(getLogLikelihood( options, chrm, thetahat, thetahatX, rhohat, rhohatX, true));
  }
  
  // *** After BurnIn ***
  if( iteration > options->getBurnIn() ){
    //logprior at estimates
    MargLikelihood->addLogPrior(LogPrior(thetahat, thetahatX, rhohat, rhohatX, options, A, rhoalpha, rhobeta, alpha) );
    double LogPosterior = 0;
    if( Populations > 1 )
      LogPosterior = CalculateLogPosterior(options, thetahat, thetahatX, rhohat, rhohatX, alpha, rhoalpha, rhobeta);
    if( A->IsRandom() ){
	  for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
	    for( int k = 0; k < Populations; k++ ){
	      vector<double> args = A->GetPriorAlleleFreqs(j,k);
	      vector<int> counts = A->GetAlleleCounts(j,k);
	      transform(counts.begin(), counts.end(), args.begin(), args.begin(), plus<double>());//PriorAlleleFreqs + AlleleCounts
	      LogPosterior += getDirichletLogDensity( args, A->getAlleleFreqsMAP(j, k) );
	    }
	  }
    }
    MargLikelihood->addLogPosteriorObs( LogPosterior );
    *SumLogLikelihood += logLikelihood;
  }
}

double Individual::LogPrior(const double* const theta, const double* const thetaX, const vector<double> rho, const vector<double> rhoX, 
			    const AdmixOptions* const options, const AlleleFreqs* const A, double rhoalpha, double rhobeta, 
			    const vector<vector<double> > &alpha)const
//Computes LogPrior at supplied parameter values
{
  int K = Populations;
  size_t theta_size = Populations;
  if(options->isRandomMatingModel()) theta_size *=2;

   double LogPrior=0.0;

   // ** case of xonly data **
   if( options->getXOnlyAnalysis() ){

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
     //prior on theta
      LogPrior += getDirichletLogDensity( alpha[0], theta );
   }
   // ** case of some x data **
   else if( Loci->isX_data() ){
     //prior on rho
      if( options->RhoFlatPrior() ){//flat prior
	LogPrior = -4.0*log( options->getTruncPt() - 1.0 );//sum of 4 flat priors
      }
      else if( options->logRhoFlatPrior() ){//flat prior on log rho
	LogPrior = -log( rho[0]*(log( options->getTruncPt() ) ) )//gamete1
	  -log( rho[1]*(log( options->getTruncPt() ) ) )         //gamete2
	  -log( rhoX[0]*(log( options->getTruncPt() ) ) )       //X chromosome
            -log( rhoX[1]*(log( options->getTruncPt() ) ) );
      }
      else{//gamma prior
	LogPrior = getGammaLogDensity( rhoalpha, rhobeta, rho[0] )//gamete1
	  + getGammaLogDensity( rhoalpha, rhobeta, rhoX[0] )      //X chromosome
	  + getGammaLogDensity( rhoalpha, rhobeta, rho[1] )       //gamete2
	  + getGammaLogDensity( rhoalpha, rhobeta, rhoX[1] );    //X chr
         LogPrior /= gsl_cdf_gamma_Q(rhobeta, rhoalpha, 1.0);
      }
     //prior on theta
      LogPrior += getDirichletLogDensity( alpha[0], theta )//first gamete
	+ getDirichletLogDensity( alpha[0], thetaX )// X chromosome
	+ getDirichletLogDensity( alpha[1], theta+K )//second gamete
	+ getDirichletLogDensity( alpha[1], thetaX + K );//Xchromosome
   }
   else{
      if( Populations > 1 ){

	if( options->isAdmixed(0) ){//gamete 1
	  //prior on rho
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
	    //prior on theta
            LogPrior += getDirichletLogDensity( alpha[0], theta );
         }
	if( options->isAdmixed(1) ){//gamete 2

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
	    //prior on theta
            LogPrior += getDirichletLogDensity( alpha[1], theta + K );
         }
      }
   }
   //prior on allele freqs, at AlleleFreqsMAP
   if( A->IsRandom() ){
      for( unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); j++ ){
         for( int k = 0; k < Populations; k++ ){
	   LogPrior += getDirichletLogDensity( A->GetPriorAlleleFreqs(j, k), A->getAlleleFreqsMAP(j,k) );
         }
      }
   }
   return LogPrior;
}



//called on first individual after burnin
double Individual::CalculateLogPosterior(const AdmixOptions* const options, const double* const theta, const double* const thetaX, 
					 const vector<double> rho, const vector<double> rhoX,
					 const vector<vector<double> > &alpha, double rhoalpha, double rhobeta)const{
  double LogPosterior = 0.0;
  double L = Loci->GetLengthOfGenome(), L_X = 0.0;
  if( Loci->isX_data() ) L_X = Loci->GetLengthOfXchrm();

  double IntConst1;
 {
   vector<double> alphaparams1(Populations), alphaparams0(Populations);
    if( options->getXOnlyAnalysis() ){
      LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN_X[0], rhobeta + L_X, rho[0] );
      if(!options->RhoFlatPrior() && !options->logRhoFlatPrior() )
	IntConst1 = IntegratingConst(rhoalpha+(double)SumN_X[0], rhobeta+L_X, 1.0, options->getTruncPt() );
      else
	IntConst1 = gsl_cdf_gamma_Q(rhobeta+L_X, rhoalpha+(double)SumN_X[0], 1.0);
      LogPosterior -= log(IntConst1);
      transform(alpha[0].begin(), alpha[0].end(), SumLocusAncestry_X, alphaparams0.begin(), std::plus<double>());
      LogPosterior += getDirichletLogDensity(alphaparams0, theta);
    }
    else if( Loci->isX_data() ){
      for( unsigned int g = 0; g < 2; g++ ){
	LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN[g], rhobeta + L, rho[g] );
	if(!options->RhoFlatPrior() && !options->logRhoFlatPrior() )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN[g], rhobeta+L, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[g], 1.0);
	LogPosterior -= log(IntConst1);
	transform(alpha[g].begin(), alpha[g].end(), SumLocusAncestry +g*Populations, alphaparams0.begin(), std::plus<double>());
	LogPosterior += getDirichletLogDensity(alphaparams0, theta + g*Populations);
      }
      for( unsigned int g = 0; g < gametes[X_posn]; g++ ){
	LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN_X[g], rhobeta + L_X, rhoX[g] );
	if( options->RhoFlatPrior() || options->logRhoFlatPrior() )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN_X[g], rhobeta+L_X, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L_X, rhoalpha+(double)SumN_X[g], 1.0);
	LogPosterior -= log(IntConst1);
	transform(alpha[g].begin(), alpha[g].end(), SumLocusAncestry_X +g*Populations, alphaparams0.begin(), std::plus<double>());
	LogPosterior += getDirichletLogDensity(alphaparams0, thetaX + g*Populations);
      }
    }
    else if( options->isSymmetric() ){//both gametes have same admixture
      vector<double> x(2,0.0);
      x[0] += getGammaLogDensity( rhoalpha + (double)SumN[0], rhobeta + L, rho[0] );
      x[1] += getGammaLogDensity( rhoalpha + (double)SumN[1], rhobeta + L, rho[0] );
      x[0] += getGammaLogDensity( rhoalpha + (double)SumN[1], rhobeta + L, rho[1] );
      x[1] += getGammaLogDensity( rhoalpha + (double)SumN[0], rhobeta + L, rho[1] );
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
      x[0] += getDirichletLogDensity(alphaparams0, theta);
      x[1] += getDirichletLogDensity(alphaparams0, theta+Populations);
      transform(alpha[1].begin(), alpha[1].end(), SumLocusAncestry+Populations, alphaparams1.begin(), std::plus<double>());
      x[0] += getDirichletLogDensity(alphaparams1, theta+Populations);
      x[1] += getDirichletLogDensity(alphaparams1, theta);
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
	LogPosterior = getGammaLogDensity( rhoalpha + (double)SumN[0], rhobeta + L, rho[0] );
	if( options->RhoFlatPrior() || options->logRhoFlatPrior() )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN[0], rhobeta+L, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[0], 1.0);
           LogPosterior -= log( IntConst1 );
	   transform(alpha[0].begin(), alpha[0].end(), SumLocusAncestry, alphaparams0.begin(), std::plus<double>());
	   LogPosterior += getDirichletLogDensity(alphaparams0, theta);
      }
      if(  options->isAdmixed(1) ){//admixed second gamete
	LogPosterior += getGammaLogDensity( rhoalpha + (double)SumN[1], rhobeta + L, rho[1] );
	if( options->RhoFlatPrior() || options->logRhoFlatPrior() )
	  IntConst1 = IntegratingConst(rhoalpha+(double)SumN[1], rhobeta+L, 1.0, options->getTruncPt() );
	else
	  IntConst1 = gsl_cdf_gamma_Q(rhobeta+L, rhoalpha+(double)SumN[1], 1.0);
	LogPosterior -= log( IntConst1 );
	transform(alpha[1].begin(), alpha[1].end(), SumLocusAncestry+Populations, alphaparams1.begin(), std::plus<double>());
	LogPosterior += getDirichletLogDensity(alphaparams1, theta+Populations);
      }
    }
  }
 return LogPosterior;
}

double Individual::IntegratingConst( double alpha, double beta, double a, double b )const
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

void Individual::OutputLikRatios(const char* const filename, int iterations, const std::string* const PopLabels){
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
