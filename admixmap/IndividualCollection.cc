/** 
 *   ADMIXMAP
 *   IndividualCollection.cc 
 *   Class to hold an array of Individuals
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
#include "IndividualCollection.h"
#include "StringSplitter.h"
#include "Regression.h"

using namespace std;

// **** CONSTRUCTORS  ****
IndividualCollection::IndividualCollection()
{
  SumLogTheta = 0;
  OutcomeType = 0;
  CovariateLabels = 0;
  NumCovariates = 0;
  NumberOfInputCovariates = 0;
  ReportedAncestry = 0;
  SumDeviance = SumDevianceSq = 0.0;
}

IndividualCollection::IndividualCollection(const AdmixOptions* const options, const InputData* const Data, const Genome& Loci, 
					   const Chromosome* const* chrm)
{
  OutcomeType = 0;
  NumOutcomes = 0;
  NumCovariates = 0;
  NumberOfInputCovariates = 0;
  indadmixoutput = 0;
  SumLogLikelihood = 0.0;
  SumDeviance = SumDevianceSq = 0.0;
  CovariateLabels = 0;
  ExpectedY = 0;
  SumResiduals = 0;
  SumLogTheta = 0;
  ReportedAncestry = 0;
  NumInd = Data->getNumberOfIndividuals();
  size = NumInd;
  NumCompLoci = Loci.GetNumberOfCompositeLoci();
  sigma.resize(2);
  sigma[0] = sigma[1] = 1.0;
  TestInd = 0;

  Individual::SetStaticMembers(&Loci, options);
  unsigned i0 = 0;
    // Fill separate individuals.
  if(options->getAnnealIndicator()==2){

    TestInd = new Individual(1, options, Data, Loci, chrm); 
    i0 = 1;
    --size;
  }

  _child = new Individual*[size];
  for (unsigned int i = 0; i < size; ++i) {
    _child[i] = new Individual(i+i0+1, options, Data, Loci, chrm);
    }
}

// ************** DESTRUCTOR **************
IndividualCollection::~IndividualCollection()
{
  Individual::DeleteStaticMembers();
  for(unsigned int i = 0; i < size; i++){
    delete _child[i];
  }
  delete[] _child;
  delete TestInd;
  delete indadmixoutput;
  delete[] OutcomeType;
  free_matrix(ExpectedY, NumOutcomes);
  free_matrix(SumResiduals, NumOutcomes);
  delete[] SumLogTheta;
  delete[] CovariateLabels;
  delete[] ReportedAncestry;
}

// ************** INITIALISATION AND LOADING OF DATA **************

void IndividualCollection::Initialise(const AdmixOptions* const options, const Genome* const Loci, const string* const PopulationLabels,
				      double rhoalpha, double rhobeta, LogWriter &Log, const DataMatrix &MLEMatrix){
  Log.setDisplayMode(Quiet);
  //Open indadmixture file  
  if ( strlen( options->getIndAdmixtureFilename() ) ){
    Log << "Writing individual-level parameters to " << options->getIndAdmixtureFilename() <<"\n";
    indadmixoutput = new IndAdmixOutputter(options, Loci, PopulationLabels);
  }
  else {
    Log << "No indadmixturefile given\n";
  }
  //Set locusfortest if specified
 if( options->getLocusForTestIndicator() )
     _locusfortest = Loci->GetChrmAndLocus( options->getLocusForTest() );

 //Initialise Admixture Proportions
 //admix_null holds initial values of admixture props, assigned to each Individual
 //indexed first by gamete, then by population
  double *admix_null;
  size_t size_admix;
  int K = options->getPopulations();
   if( options->isRandomMatingModel() )
     size_admix = K*2; // double the size for 2 gametes in RMM
   else//assortative mating
     size_admix = K;

   admix_null = new double[size_admix];
   for(unsigned i = 0; i< size_admix; ++i) admix_null[i] = (double)1.0;
   
   vector<double >alphatemp(K);
   
   int KK0=K, KK1=K;
   if( options->sizeInitAlpha() == 0 ){
     fill(alphatemp.begin(), alphatemp.end(), 1.0 );
   }
   else if( options->sizeInitAlpha() == 1 ){
    alphatemp = options->getInitAlpha(0);
   }
   else  if( !options->getIndAdmixHierIndicator() ){
     alphatemp = options->getInitAlpha(0);
     
     for( int k = 0; k < K; k++ ){
       if( alphatemp[k] == 0 ) {
	 admix_null[k] = 0.0;
	 KK0--;
       }
     }
     
     alphatemp = options->getInitAlpha(1);
     
     //possible problem if not randommating model ?
     for( int k = 0; k < K; k++ ){
       if( alphatemp[k] == 0 ) {
	 admix_null[K + k] = 0.0;
	 KK1--;
       }
     }
   }

   //set initial values of admixture to 1 / (number of admixed pops)
   // eg 1 1 0 -> 1/2, 1/2, 1/2
   // should draw initial values randomly from the prior so that different chains will have 
   // different starting values
   for( int k = 0; k < K; k++ )admix_null[k] /= KK0;//gamete 1
   if( options->isRandomMatingModel() )   for( int k = 0; k < K; k++ )admix_null[K + k] /= KK1;//gamete 2
   
   setAdmixtureProps(admix_null, size_admix);
   if( Loci->isX_data() )setAdmixturePropsX(admix_null, size_admix);
   delete[] admix_null;
   
  //Regression stuff
   if(options->getNumberOfOutcomes() > 0){
    ExpectedY = alloc2D_d(NumOutcomes, NumInd);
    SumResiduals = alloc2D_d(NumOutcomes, NumInd);
    for(int j = 0; j < NumOutcomes; ++j)fill(SumResiduals[j], SumResiduals[j] + NumInd, 0.0);
   }

  SumLogTheta = new double[ options->getPopulations()];
  if( options->getMLIndicator() )
    InitialiseMLEs(rhoalpha,rhobeta,options, MLEMatrix);
  //set to very large negative value (effectively -Inf) so the first value is guaranteed to be greater
  MaxLogLikelihood.assign(NumInd, -9999999 );
}

// ** this function needs debugging
void IndividualCollection::InitialiseMLEs(double rhoalpha, double rhobeta, const AdmixOptions* const options, 
					  const DataMatrix &MLEMatrix){
  //set thetahat and rhohat, estimates of individual admixture and sumintensities

   size_t size_admix;
   int K = options->getPopulations();
   if( options->isRandomMatingModel() )
     size_admix = K*2; // double the size for 2 gametes in RMM
   else//assortative mating
     size_admix = K;

   thetahat = new double[size_admix];
   thetahatX = new double[size_admix];

   //initialise thetahat at initial values of individual admixture
   for(unsigned k = 0; k < size_admix; ++k)
     thetahat[k] = _child[0]->getAdmixtureProps()[k];

   if(!options->isGlobalRho())
     //initialise rhohat at initial value of Individual 1's sumintensities
      rhohat = _child[0]->getRho();
   else{
   //initialise rhohat at initial value of globalsumintensities ie prior mean
     vector<double> r(2, rhoalpha/rhobeta );
     rhohat = r;
     rhohatX = r;
   }
      //TODO: X objects


   //use previously read values from file, if available
   if( NumInd == 1 && strlen(options->getMLEFilename())>0 ){
      rhohat[0] = MLEMatrix.get( options->getPopulations(), 0 );
      if( options->isXOnlyAnalysis() ){
	for(int k = 0; k < options->getPopulations(); ++k) thetahat[k] = MLEMatrix.get(k,0);
      }
      else{
	for(int k = 0; k < options->getPopulations(); ++k) {
	  thetahat[k] = MLEMatrix.get(k,0);
	  thetahat[k+ options->getPopulations()] = MLEMatrix.get(k,1);
	}
	rhohat[1] = MLEMatrix.get(options->getPopulations(), 1 );
      }
      setAdmixtureProps(thetahat, size_admix);
   }
}

void IndividualCollection::LoadData(const AdmixOptions* const options, const InputData* const data_){
  LoadCovariates(data_, options);
  
  if ( strlen( options->getOutcomeVarFilename() ) != 0 ){
    LoadOutcomeVar(data_);
  }
  if ( strlen( options->getReportedAncestryFilename() ) != 0 ){
    LoadRepAncestry(data_);
  }
}

void IndividualCollection::LoadCovariates(const InputData* const data_, const AdmixOptions* const options){
  if ( strlen( options->getCovariatesFilename() ) > 0 ){
    DataMatrix& CovData = (DataMatrix&)data_->getCovariatesMatrix();
    NumberOfInputCovariates = CovData.nCols();

    if( !options->getTestForAdmixtureAssociation() && options->getPopulations() > 1 ){
      Covariates.setDimensions(NumInd, CovData.nCols() + options->getPopulations());
      for(unsigned i = 0; i < NumInd; ++i)for(int j = 0; j < options->getPopulations()-1; ++j)
         Covariates.set(i, j + CovData.nCols() + 1,  1.0 / (double)options->getPopulations() );
      }
    else
      Covariates.setDimensions(NumInd, CovData.nCols()+1);
    for(unsigned i = 0; i < NumInd; ++i)for(unsigned j = 0; j < CovData.nCols(); ++j)
      {
      Covariates.set(i,0, 1.0); //for intercept
      Covariates.set(i,j+1, CovData.get(i+1, j) );
      Covariates.isMissing(i,j+1, CovData.isMissing(i+1,j));
      }
    if ( Covariates.hasMissing() ) Covariates.SetMissingValuesToColumnMeans();
    CovariateLabels = new string[ NumberOfInputCovariates ];
    getLabels(data_->getInputData()[0], CovariateLabels);

    //vector<double> mean(NumberOfInputCovariates);

    //centre covariates about their means
//    for( int j = 0; j < NumberOfInputCovariates; j++ ){
  //    int count = 0;
    //  mean[j] = 0.0;
      //for(unsigned int i = 0; i < NumInd; i++ )
        //if(!Covariates.IsMissingValue(i,j+1)){
          // mean[j] += Covariates(i,j+1);
           //++count;
        //}
      //mean[j] /= (double)count;
      //}

    //for(unsigned int i = 0; i < NumInd; i++ )
  	 // for( int j = 0; j < NumberOfInputCovariates; j++ )
	    //Covariates( i, j+1 ) -= mean[j];
  }
  else {
    if( !options->getTestForAdmixtureAssociation() && options->getPopulations() > 1 ){
      Covariates.setDimensions(NumInd, options->getPopulations());
        for(unsigned i = 0; i < NumInd; ++i)for(int j = 1; j < options->getPopulations(); ++j)
          Covariates.set(i, j, 1.0 / (double)options->getPopulations() );
    }
    else
      Covariates.setDimensions(NumInd, 1);
    for(unsigned i = 0; i < NumInd; ++i)
    Covariates.set(i,0, 1.0 );
    }

  if( options->getTestForAdmixtureAssociation() )
       NumCovariates = NumberOfInputCovariates + 1;
  else
       NumCovariates = NumberOfInputCovariates + options->getPopulations();

}

void IndividualCollection::LoadOutcomeVar(const InputData* const data_){
  DataMatrix& OutcomeVarData = (DataMatrix&)data_->getOutcomeVarMatrix();
  NumOutcomes = OutcomeVarData.nCols();
  Outcome = data_->getOutcomeVarMatrix();;
 
  delete[] OutcomeType;
  OutcomeType = new DataType[ NumOutcomes ];
  data_->getOutcomeTypes(OutcomeType);
}

void IndividualCollection::LoadRepAncestry(const InputData* const data_){
  ReportedAncestry = new DataMatrix[NumInd];
  DataMatrix& temporary = (DataMatrix&)data_->getReportedAncestryMatrix();
  for( unsigned i = 0; i < temporary.nRows() / 2; i++ )
    ReportedAncestry[i] = temporary.SubMatrix( 2*i, 2*i + 1, 0, temporary.nCols() - 1 );
 
}
void IndividualCollection::getLabels(const Vector_s& data, string *labels)
{
  for (size_t i = 0, index = 0; i < data.size(); ++i) {
    labels[index++] = data[i];
  }
}
void IndividualCollection::setAdmixtureProps(const double* const a, size_t thetasize)
{
  if(TestInd)
    TestInd->setAdmixtureProps(a, thetasize);
  for(unsigned int i = 0; i < size; i++){
    _child[i]->setAdmixtureProps(a, thetasize);
  }
}

void IndividualCollection::setAdmixturePropsX(const double* const a, size_t thetasize)
{
  if(TestInd)
    TestInd->setAdmixtureProps(a, thetasize);
  for(unsigned int i = 0; i < size; i++){
    _child[i]->setAdmixturePropsX(a, thetasize);
  }
}

void IndividualCollection::SetExpectedY(int k, const double* const beta){
  //sets ExpectedY = X * Beta
  if(ExpectedY){
    matrix_product(Covariates.getData(), beta, ExpectedY[k], Covariates.nRows(), Covariates.nCols(), 1);
    if(OutcomeType[k] == Binary)
      //for binary outcome sets EY as logit^-1(X*beta)
      for(unsigned int i = 0; i < NumInd; i++ )
	ExpectedY[k][i] = 1 / ( 1 + exp( -ExpectedY[k][i] ) );
    }
}

// void IndividualCollection::calculateExpectedY(int k)
// {
//   if(ExpectedY)
//     for(unsigned int i = 0; i < NumInd; i++ )
//       ExpectedY[k][i] = 1 / ( 1 + exp( -ExpectedY[k][i] ) );
// }

void IndividualCollection::UpdateSumResiduals(){
  if(SumResiduals)
    for(int k = 0; k < NumOutcomes; ++k)
      for(unsigned i = 0; i < NumInd; ++i)
	SumResiduals[k][i] += Outcome.get(i,k) - ExpectedY[k][i];
}

void IndividualCollection::HMMIsBad(bool b){
  if(TestInd)TestInd->HMMIsBad(b);
  for(unsigned i = 0; i < size; ++i)
    _child[i]->HMMIsBad(b);
}

// ************** UPDATING **************
void IndividualCollection::Update(int iteration, const AdmixOptions* const options, Chromosome **chrm, AlleleFreqs *A,
				  const Regression* const R, const double* const poptheta,
				  const vector<vector<double> > &alpha, double globalrho,
				  double rhoalpha, double rhobeta, LogWriter &Log, double coolness, bool anneal=false){
  fill(SumLogTheta, SumLogTheta+options->getPopulations(), 0.0);//reset to 0
  if(iteration > options->getBurnIn())Individual::ResetScores(options);

  vector<double> lambda;
  vector<const double*> beta;
  for(int i = 0; i < options->getNumberOfOutcomes(); ++i){
    lambda.push_back( R[i].getlambda());
    beta.push_back( R[i].getbeta());
  }
  int i0 = 0;
  if(options->getAnnealIndicator()==2){
    i0 = 1;
    Chromosome::setCoolness(coolness);
    TestInd->SampleParameters(SumLogTheta, A, iteration , &Outcome, OutcomeType, ExpectedY,
			      lambda, NumCovariates, &Covariates, beta, poptheta, options,
			      chrm, alpha, rhoalpha, rhobeta, sigma,  
			      DerivativeInverseLinkFunction(0),
			      R[0].getDispersion(), anneal );
    if(options->getPopulations() > 1 && (iteration %2))//conjugate update of theta on even-numbered iterations
      TestInd->SampleTheta(iteration, SumLogTheta, &Outcome, chrm, OutcomeType, ExpectedY, lambda, NumCovariates,
			   & Covariates, beta, poptheta, options, alpha, sigma,
			   DerivativeInverseLinkFunction(0), 
			   R[0].getDispersion(), false, anneal);
    
    _child[size-1]->HMMIsBad(false);//The HMMs are shared between individuals so if there are two or more individuals
    //an update of one will overwrite the HMM and Chromosome values for the other. However, the stored value of loglikelihood will
    //still be valid as the allelefreqs, global rho and that individual's have not changed.
    
    Chromosome::setCoolness(1.0); 
    TestInd->HMMIsBad(false);  //TODO:
  }

  for(unsigned int i = 0; i < size; i++ ){
    int prev = i-1;
    if(i==0)prev = size-1;
    
    _child[i]->SampleParameters(SumLogTheta, A, iteration , &Outcome, OutcomeType, ExpectedY,
				lambda, NumCovariates, &Covariates, beta, poptheta, options,
				chrm, alpha, rhoalpha, rhobeta, sigma,  
				DerivativeInverseLinkFunction(i+i0),
				R[0].getDispersion(), anneal );
    if(options->getPopulations() > 1 && (iteration %2))//conjugate update of theta on even-numbered iterations
      _child[i]->SampleTheta(iteration, SumLogTheta, &Outcome, chrm, OutcomeType, ExpectedY, lambda, NumCovariates,
			     & Covariates, beta, poptheta, options, alpha, sigma,
			     DerivativeInverseLinkFunction(i+i0), 
			     R[0].getDispersion(), false, anneal);
    
    if(size > 1)_child[prev]->HMMIsBad(false);//The HMMs are shared between individuals so if there are two or more individuals
    //an update of one will overwrite the HMM and Chromosome values for the other. However, the stored value of loglikelihood will
    //still be valid until allelefreqs, global rho or individual admixture are updated.
    
    if( options->getMLIndicator() && (i == 0) )//compute marginal likelihood for first individual
      _child[i]->Chib(iteration, &SumLogLikelihood, &(MaxLogLikelihood[i]),
				options, chrm, alpha, globalrho, rhoalpha, rhobeta,
				thetahat, thetahatX, rhohat, rhohatX, Log, &MargLikelihood, A);
  }
}

void IndividualCollection::setGenotypeProbs(Chromosome** C, unsigned nchr){
  if(TestInd)
    for(unsigned j = 0; j < nchr; ++j)
      TestInd->SetGenotypeProbs(j, C[j],false);
  for(unsigned int i = 0; i < size; i++ ){
    for(unsigned j = 0; j < nchr; ++j)
      _child[i]->SetGenotypeProbs(j, C[j],false);
  }

}
// ************** ACCESSORS **************
int IndividualCollection::getSize()const
{
  return size;
}

double IndividualCollection::GetSumrho()const
{
   double Sumrho = 0;
   for( unsigned int i = 0; i < size; i++ )
      Sumrho += (*_child[i]).getSumrho();
   return Sumrho;
}

vector<double> IndividualCollection::getOutcome(int j)const{
  return Outcome.getCol(j);
}

double IndividualCollection::getOutcome(int j, int ind)const{
    return Outcome.get(ind, j);
}

int IndividualCollection::getNumberOfOutcomeVars()const{
  return NumOutcomes;
}

DataType IndividualCollection::getOutcomeType(int i)const{
  return OutcomeType[i];
}

Individual* IndividualCollection::getIndividual(int num)const
{
  if (num < (int)size){
    return _child[num];
  } else {
    return 0;
  }
}

int IndividualCollection::GetNumberOfInputCovariates()const{
  return NumberOfInputCovariates;
}
int IndividualCollection::GetNumCovariates() const{
  return NumCovariates;
}

const double* IndividualCollection::getCovariates()const{
  return Covariates.getData();
}

const std::string IndividualCollection::getCovariateLabels(int i)const{
  return CovariateLabels[i];
  }
const std::string *IndividualCollection::getCovariateLabels()const{
  return CovariateLabels;
  }

double IndividualCollection::getExpectedY(int i)const{
  if(ExpectedY)
    return ExpectedY[0][i];
  else
    return 0.0;
 }
double IndividualCollection::getSumLogTheta(int i)const{
  return SumLogTheta[i];
}
const double* IndividualCollection::getSumLogTheta()const{
  return SumLogTheta;
}
const chib* IndividualCollection::getChib()const{
  return &MargLikelihood;
}
// ************** OUTPUT **************

void IndividualCollection::OutputDeviance(const AdmixOptions* const options, Chromosome** C, Regression *R, 
					  LogWriter &Log, double SumLogRho, unsigned numChromosomes){
  //SumRho = ergodic sum of global sumintensities

  int iterations = options->getTotalSamples()-options->getBurnIn();
  double E, V;
  E = SumDeviance / (double) iterations;//ergodic average of deviance
  V = SumDevianceSq / (double)iterations - E*E;//ergodic variance of deviance 

  if(options->getDisplayLevel()>0)Log.setDisplayMode(On);
  else Log.setDisplayMode(Off);
  Log << "\nMeanDeviance(D_bar)\t" << E << "\n"
      << "VarDeviance(V)\t" << V << "\n"
      << "PritchardStat(D_bar+0.25V)\t" << E + 0.25 *V << "\n";

  //update chromosomes using globalrho, for globalrho model
  if(options->isGlobalRho())
    for( unsigned int j = 0; j < numChromosomes; j++ )
      C[j]->SetLociCorr(exp(SumLogRho / (double)iterations));

  //set haplotype pair probs to posterior means
  for( unsigned int j = 0; j < numChromosomes; j++ )
    for(unsigned jj = 0; jj < C[j]->GetNumberOfCompositeLoci(); ++jj)
      (*C[j])(jj)->SetHapPairProbsToPosteriorMeans(iterations);

  //set genotype probs using happair probs calculated at posterior means of allele freqs 
  setGenotypeProbs(C, numChromosomes);

  //accumulate deviance at posterior means for each individual
  double Lhat = 0.0; // Lhat = loglikelihood at estimates
  for(unsigned int i = 0; i < size; i++ ){
    Lhat += _child[i]->getLogLikelihoodAtPosteriorMeans(options, C);
  }

  Log << "DevianceAtPosteriorMean(IndAdmixture)" << -2.0*Lhat << "\n";
  for(int c = 0; c < options->getNumberOfOutcomes(); ++c){
    double RegressionLogL = R[c].getLogLikelihoodAtPosteriorMeans(this, iterations);
    Lhat += RegressionLogL;
    Log << "DevianceAtPosteriorMean(Regression " << c << ")"
	<< -2.0*RegressionLogL << "\n";
  }

  double pD = E + 2.0*Lhat;
  double DIC = E + pD;

  Log << "DevianceAtPosteriorMean(D_hat)\t" << -2.0*Lhat << "\n"
      << "EffectiveNumParameters(pD)\t" << pD << "\n"
      << "DevianceInformationCriterion\t" << DIC << "\n\n"; 
}

void IndividualCollection::OutputIndAdmixture()
{
  indadmixoutput->visitIndividualCollection(*this);
  if(TestInd)
    indadmixoutput->visitIndividual(*TestInd, _locusfortest);
  for(unsigned int i = 0; i < size; i++){
    indadmixoutput->visitIndividual(*_child[i], _locusfortest);
  }
}

void IndividualCollection::OutputChibEstimates(LogWriter &Log, int Populations)const{
  //Used only if marglikelihood = 1
  Log.setDisplayMode(Off);
  Log << "EstimatesUsedForChibAlgorithm\t";

  for(int k = 0; k < Populations; ++k){
    Log << thetahat[k] << "\t";
  }
  for(int k = 0; k < Populations; ++k){
    Log << thetahat[Populations +k] << "\t";
  }
  Log << rhohat[0] <<"\t" << rhohat[1] << "\n";
}
void IndividualCollection::OutputChibResults(LogWriter& Log)const{
  Log.setDisplayMode(On);
  Log << "\nChib values at estimates:"
      << "\nDeviance\t" << -2.0*MargLikelihood.getLogLikelihood()
      << "\nLogLikelihood\t" << MargLikelihood.getLogLikelihood()
      << "\nLogPrior\t" << MargLikelihood.getLogPrior()
      << "\nLogPosterior\t" << MargLikelihood.getLogPosterior()
      << "\n\nLogMarginalLikelihood\t" << MargLikelihood.getLogMarginalLikelihood()
      << "\n\n";
}

//for single individual analysis
void IndividualCollection::OutputErgodicAvg(int samples, bool ML, std::ofstream *avgstream){
  double E, V;
  E = SumDeviance / (double) samples;//ergodic average of deviance
  V = SumDevianceSq / (double)samples - E*E;//ergodic variance of deviance 
  *avgstream << E<< " "<<V<<" ";
  if(ML)
    *avgstream //<< SumLogLikelihood / samples << " "
      << MargLikelihood.getLogPrior()<< " " << MargLikelihood.getLogPosterior() << " "
      << _child[0]->getLogPosteriorTheta() << " " << _child[0]->getLogPosteriorRho()<< " " 
      << _child[0]->getLogPosteriorAlleleFreqs() << " "
      << MargLikelihood.getLogMarginalLikelihood();
}

void IndividualCollection::OutputResiduals(const char* ResidualFilename, const Vector_s Labels, int iterations){
  std::ofstream ResidualStream(ResidualFilename, ios::out);
  if( !ResidualStream )
    {
      cerr<< "WARNING: Couldn't open residualfile\n";
    }
  else{
    for(int j = 0; j < NumOutcomes; ++j)
      ResidualStream << Labels[j]<< "\t";
    ResidualStream << endl;
    for(unsigned i = 0; i < size; ++i){
      for(int j = 0; j < NumOutcomes; ++j)
	ResidualStream << SumResiduals[j][i] / (double) iterations << "\t";
      ResidualStream << endl;
    }
    ResidualStream.close();
  }
}

void IndividualCollection::getOnePopOneIndLogLikelihood(LogWriter &Log, const string* const PopulationLabels)
{
  Log.setDisplayMode(On);
  Log << "Log-likelihood for unadmixed "  << (*PopulationLabels)[0] << ": "
      << _child[0]->getLogLikelihoodOnePop() << "\n";

}
double IndividualCollection::getLogLikelihood(const AdmixOptions* const options, Chromosome **C, const Regression* R, 
					      bool sumdeviance=false){
  double LogL = 0.0;
  for(unsigned i = 0; i < size; ++i){
    int prev = i-1;
    if(i==0)prev = size-1;
    LogL += _child[i]->getLogLikelihood(options, C);
    if(size > 1)_child[prev]->HMMIsBad(false);
  }
  double LogLikelihood = LogL;
  //accumulate deviance
  if(sumdeviance){
    for(int c = 0; c < options->getNumberOfOutcomes(); ++c)LogLikelihood += R[c].getLogLikelihood(this);
    SumDeviance += -2.0*LogLikelihood;
    SumDevianceSq += 4.0*LogLikelihood*LogLikelihood;
  }

  return LogL;
}
double IndividualCollection::getModifiedLogLikelihood(const AdmixOptions* const options, Chromosome **C, double coolness){
  double LogL = 0.0;
  Chromosome::setCoolness(1.0);

  if(options->getAnnealIndicator()==2){
    TestInd->HMMIsBad(true);//to force HMM update and recalculation of log likelihood
    LogL += TestInd->getLogLikelihood(options, C);
    _child[size-1]->HMMIsBad(false);
    //if(coolness < 1.0)
    TestInd->HMMIsBad(true);
  }
  else if(options->getAnnealIndicator()==1){
    for(unsigned i = 0; i < size; ++i){
      int prev = i-1;
      if(i==0)prev = size-1;
      _child[i]->HMMIsBad(true);//to force HMM update and recalculation of log likelihood
      LogL += _child[i]->getLogLikelihood(options, C);
      if(size > 1)_child[prev]->HMMIsBad(false);
      if(coolness < 1.0)_child[i]->HMMIsBad(true);
    }
    Chromosome::setCoolness(coolness);
  }

  return LogL;
}

//returns Derivative of Inverse Link Function for individual i
double IndividualCollection::DerivativeInverseLinkFunction(int i)const{
  double DInvLink = 1.0;
  
  if(OutcomeType){//in case no regression model
    if(OutcomeType[0] == Binary)
      DInvLink = ExpectedY[0][i] * (1.0 - ExpectedY[0][i]);
    else if(OutcomeType[0] == Continuous)DInvLink = 1.0;
  }
  return DInvLink;    
}
void IndividualCollection::ResetChib(){
  MargLikelihood.Reset();
}
