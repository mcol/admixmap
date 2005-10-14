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

IndividualCollection::IndividualCollection(AdmixOptions* options,InputData *Data, Genome& Loci, Chromosome **chrm)
{
  OutcomeType = 0;
  NumOutcomes = 0;
  NumCovariates = 0;
  NumberOfInputCovariates = 0;
  indadmixoutput = 0;
  LogLikelihood=0.0;
  SumLogLikelihood = 0.0;
  SumDeviance = SumDevianceSq = 0.0;
  CovariateLabels = 0;
  ExpectedY = 0;
  SumLogTheta = 0;
  ReportedAncestry = 0;
  NumInd = Data->getNumberOfIndividuals();
  NumCompLoci = Loci.GetNumberOfCompositeLoci();
  sigma.resize(2);
  sigma[0] = sigma[1] = 1.0;

  _child = new Individual*[NumInd];
  Individual::SetStaticMembers(&Loci, options);
    // Fill separate individuals.
  for (unsigned int i = 0; i < NumInd; ++i) {
    _child[i] = new Individual(i+1,options, Data, Loci, chrm);
    }
}

// ************** DESTRUCTOR **************
IndividualCollection::~IndividualCollection()
{
  Individual::DeleteStaticMembers();
  for(unsigned int i = 0; i < NumInd; i++){
    delete _child[i];
  }
  delete[] _child;
  delete indadmixoutput;
  delete[] OutcomeType;
  free_matrix(ExpectedY, NumOutcomes);
  delete[] SumLogTheta;
  delete[] CovariateLabels;
  delete[] ReportedAncestry;
}

// ************** INITIALISATION AND LOADING OF DATA **************

void IndividualCollection::Initialise(AdmixOptions *options, Genome *Loci, std::string *PopulationLabels,
				      double rhoalpha, double rhobeta, LogWriter *Log, const DataMatrix &MLEMatrix){
  //Open indadmixture file  
  if ( strlen( options->getIndAdmixtureFilename() ) ){
    Log->logmsg(true,"Writing individual-level parameters to ");
    Log->logmsg(true,options->getIndAdmixtureFilename());
    Log->logmsg(true,"\n");
    indadmixoutput = new IndAdmixOutputter(options,Loci,PopulationLabels);
  }
  else {
    Log->logmsg(true,"No indadmixturefile given\n");
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
   for( int k = 0; k < K; k++ )admix_null[k] /= KK0;//gamete 1
   if( options->isRandomMatingModel() )   for( int k = 0; k < K; k++ )admix_null[K + k] /= KK1;//gamete 2
   
   setAdmixtureProps(admix_null, size_admix);
   if( Loci->isX_data() )setAdmixturePropsX(admix_null, size_admix);
   delete[] admix_null;
   
  //Regression stuff
   if(options->getNumberOfOutcomes() > 0){
    ExpectedY = alloc2D_d(NumOutcomes, NumInd);
   }

  SumLogTheta = new double[ options->getPopulations()];
  if( options->getMLIndicator() )
    InitialiseMLEs(rhoalpha,rhobeta,options, MLEMatrix);
  //set to very large negative value (effectively -Inf) so the first value is guaranteed to be greater
  MaxLogLikelihood.assign(NumInd, -9999999 );
}

// ** this function needs debugging
void IndividualCollection::InitialiseMLEs(double rhoalpha, double rhobeta, AdmixOptions * options, const DataMatrix &MLEMatrix){
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
      if( options->getXOnlyAnalysis() ){
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

void IndividualCollection::LoadData(AdmixOptions *options, InputData *data_){
  LoadCovariates(data_, options);
  
  if ( strlen( options->getOutcomeVarFilename() ) != 0 ){
    LoadOutcomeVar(data_);
  }
  if ( strlen( options->getReportedAncestryFilename() ) != 0 ){
    LoadRepAncestry(data_);
  }
}

void IndividualCollection::LoadCovariates(InputData *data_, AdmixOptions *options){
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

void IndividualCollection::LoadOutcomeVar(InputData *data_){
  DataMatrix& OutcomeVarData = (DataMatrix&)data_->getOutcomeVarMatrix();
  NumOutcomes = OutcomeVarData.nCols();
  Outcome = data_->getOutcomeVarMatrix();;
 
  delete[] OutcomeType;
  OutcomeType = new DataType[ NumOutcomes ];
  data_->getOutcomeTypes(OutcomeType);
}

void IndividualCollection::LoadRepAncestry(InputData *data_){
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
void IndividualCollection::setAdmixtureProps(double *a, size_t size)
{
  for(unsigned int i=0; i<NumInd; i++){
    _child[i]->setAdmixtureProps(a, size);
  }
}

void IndividualCollection::setAdmixturePropsX(double *a, size_t size)
{
  for(unsigned int i = 0; i < NumInd; i++){
    _child[i]->setAdmixturePropsX(a, size);
  }
}

void IndividualCollection::SetExpectedY(int k, const double* const beta){
  //sets ExpectedY = X * Beta
  if(ExpectedY){
    matrix_product(Covariates.getData(), beta, ExpectedY[k], Covariates.nRows(), Covariates.nCols(), 1);
    }
}

void IndividualCollection::calculateExpectedY( int k)
{
  if(ExpectedY)
    for(unsigned int i = 0; i < NumInd; i++ )
      ExpectedY[k][i] = 1 / ( 1 + exp( -ExpectedY[k][i] ) );
}

// ************** UPDATING **************
void IndividualCollection::Update(int iteration, AdmixOptions *options, Chromosome **chrm, AlleleFreqs *A,
Regression *R0, Regression *R1, const double *poptheta,
				   vector<vector<double> > &alpha, double globalrho,
				  double rhoalpha, double rhobeta, LogWriter *Log){
  fill(SumLogTheta, SumLogTheta+options->getPopulations(), 0.0);//reset to 0
  if(iteration > options->getBurnIn())Individual::ResetScores(options);

  double lambda[] = {R0->getlambda(), R1->getlambda()};
  double *beta[] = {R0->getbeta(), R1->getbeta()};
  LogLikelihood = 0.0; 

  for(unsigned int i = 0; i < NumInd; i++ ){
    int prev = i-1;
    if(i==0)prev = NumInd-1;

    if( options->getPopulations() > 1 ){
      _child[i]->SampleParameters(i, SumLogTheta, &LogLikelihood, A, iteration , &Outcome, NumOutcomes, OutcomeType, ExpectedY,
				  lambda, NumCovariates, &Covariates, beta, poptheta, options,
				  chrm, alpha, rhoalpha, rhobeta, sigma,  
				  DerivativeInverseLinkFunction(i),
				  R0->getDispersion()
				  );
      if((iteration %2))//conjugate update of theta on even-numbered iterations
 	_child[i]->SampleTheta(i, iteration, SumLogTheta, &Outcome, chrm, NumOutcomes, OutcomeType, ExpectedY, lambda, NumCovariates,
 			      & Covariates, beta, poptheta, options, alpha, sigma,
 			       DerivativeInverseLinkFunction(i), 
 			       R0->getDispersion(), false);

      if(NumInd > 1)_child[prev]->HMMIsBad(false);//The HMMs are shared between individuals so if there are two or more individuals
      //an update of one will overwrite the HMM and Chromosome values for the other. However, the stored value of loglikelihood will
      //still be valid as the allelefreqs, global rho and that individual's have not changed.
    }
    
    else{//single population 
      _child[i]->OnePopulationUpdate(i, (bool)(iteration > options->getBurnIn()), &Outcome, NumOutcomes, OutcomeType, ExpectedY, lambda,
				     chrm, A);
      LogLikelihood += _child[i]->getLogLikelihoodOnePop(false);
    }   

    if( options->getMLIndicator() && (i == 0) )//compute marginal likelihood for first individual
      _child[i]->Chib(iteration, &SumLogLikelihood, &(MaxLogLikelihood[i]),
				options, chrm, alpha, globalrho, rhoalpha, rhobeta,
				thetahat, thetahatX, rhohat, rhohatX, Log, &MargLikelihood, A);
  }
  if(iteration > options->getBurnIn()){
    SumDeviance += -2.0*LogLikelihood;
    SumDevianceSq += 4.0*LogLikelihood*LogLikelihood;
  }

}

void IndividualCollection::ConjugateUpdateIndAdmixture(int iteration, Regression *R0, Regression *R1, const double *poptheta, 
						       AdmixOptions *options, Chromosome **chrm, vector<vector<double> > &alpha){
  if( options->getPopulations() > 1 ){
    double lambda[] = {R0->getlambda(), R1->getlambda()};
    double *beta[] = {R0->getbeta(), R1->getbeta()};

    for(unsigned int i = 0; i < NumInd; i++ )
      _child[i]->SampleTheta(i, iteration, SumLogTheta, &Outcome, chrm, NumOutcomes, OutcomeType, ExpectedY, lambda, NumCovariates,
			    &Covariates, beta, poptheta, options, alpha, sigma,
			     DerivativeInverseLinkFunction(i), 
			     R0->getDispersion(), false);
  }
}

// ************** ACCESSORS **************
int IndividualCollection::getSize()
{
  return NumInd;
}

double IndividualCollection::GetSumrho()
{
   double Sumrho = 0;
   for( unsigned int i = 0; i < NumInd; i++ )
      Sumrho += (*_child[i]).getSumrho();
   return Sumrho;
}

vector<double> IndividualCollection::getOutcome(int j){
  return Outcome.getCol(j);
}

double IndividualCollection::getOutcome(int j, int ind){
    return Outcome.get(ind, j);
}

int IndividualCollection::getNumberOfOutcomeVars(){
  return NumOutcomes;
}

DataType IndividualCollection::getOutcomeType(int i){
  return OutcomeType[i];
}

Individual* IndividualCollection::getIndividual(int num)
{
  if (num < (int)NumInd){
    return _child[num];
  } else {
    return 0;
  }
}

int IndividualCollection::GetNumberOfInputCovariates(){
  return NumberOfInputCovariates;
}
int IndividualCollection::GetNumCovariates() const{
  return NumCovariates;
}

const double* IndividualCollection::getCovariates()const{
  return Covariates.getData();
}

std::string IndividualCollection::getCovariateLabels(int i){
  return CovariateLabels[i];
  }
std::string *IndividualCollection::getCovariateLabels(){
  return CovariateLabels;
  }

double IndividualCollection::getExpectedY(int i){
  if(ExpectedY)
    return ExpectedY[0][i];
  else
    return 0.0;
 }

// ************** OUTPUT **************

void IndividualCollection::OutputDeviance(AdmixOptions *options, Chromosome** C, LogWriter *Log, double SumRho, unsigned numChromosomes){
  //SumRho = ergodic sum of global sumintensities

  int iterations = options->getTotalSamples()-options->getBurnIn();
  double E, V;
  E = SumDeviance / (double) iterations;//ergodic average of deviance
  V = SumDevianceSq / (double)iterations - E*E;//ergodic variance of deviance 

  Log->logmsg(true, "MeanDeviance    VarDeviance     GOFStat         pD      DIC\n");
  Log->logmsg(true, E);Log->logmsg(true,"    ");
  Log->logmsg(true, V);Log->logmsg(true,"    ");
  Log->logmsg(true, E + 0.25 *V);Log->logmsg(true, "    ");

  //update chromosomes using globalrho, for globalrho model
  if(options->isGlobalRho())
    for( unsigned int j = 0; j < numChromosomes; j++ )
      C[j]->SetLociCorr(SumRho / (double)iterations);

  //accumulate deviance at posterior means for each individual
  double D = 0.0; // D = deviance at estimates
  if(options->getPopulations() > 1)
    for(unsigned int i = 0; i < NumInd; i++ ){
      D += -2.0*_child[i]->getLogLikelihoodAtPosteriorMeans(options, C);
    }
  else//single population
    for(unsigned int i = 0; i < NumInd; i++ ){
      D += -2.0*_child[i]->getLogLikelihoodAtPosteriorMeansOnePop(iterations);
    }
  double pD = E - D;
  double DIC = E + pD;
  
  Log->logmsg(true, pD);Log->logmsg(true, "    ");
  Log->logmsg(true, DIC);Log->logmsg(true, "\n\n");
}

void IndividualCollection::OutputIndAdmixture()
{
  indadmixoutput->visitIndividualCollection(*this);
  for(unsigned int i = 0; i < NumInd; i++){
    indadmixoutput->visitIndividual(*_child[i], _locusfortest, LogLikelihood);
  }
}

void IndividualCollection::OutputChibEstimates(LogWriter *Log, int Populations){
  //Used only if marglikelihood = 1
  Log->write("Estimates used in Chib algorithm to estimate marginal likelihood for Individual 1:\n");

  for(int k = 0; k < Populations; ++k)
    Log->write(thetahat[k]);
  for(int k = 0; k < Populations; ++k)
    Log->write(thetahat[Populations +k]);
  Log->write( rhohat[0]);Log->write( rhohat[1]);Log->write("\n");

  Log->logmsg(true, "Individual 1:\n");
  Log->logmsg(true,   "Log likelihood (at estimates): ");Log->logmsg(true, MargLikelihood.getLogLikelihood());
  Log->logmsg(true, "\nLog prior      (at estimates): ");Log->logmsg(true, MargLikelihood.getLogPrior());
  Log->logmsg(true, "\nLog posterior  (at estimates): ");Log->logmsg(true, MargLikelihood.getLogPosterior());
  Log->logmsg(true, "\nLog marginal likelihood : ");Log->logmsg(true, MargLikelihood.getLogMarginalLikelihood());
  Log->logmsg(true, "\n\n");
}

//for single individual analysis
void IndividualCollection::OutputErgodicAvg(int samples, std::ofstream *avgstream){
     *avgstream << SumLogLikelihood / samples << " "
               << MargLikelihood.getLogPosterior();
}

void IndividualCollection::getOnePopOneIndLogLikelihood(LogWriter *Log, std::string *PopulationLabels)
{
   Log->logmsg(true,"Log-likelihood for unadmixed ");
   Log->logmsg(true, (*PopulationLabels)[0]);
   Log->logmsg(true, ": ");
   Log->logmsg(true, _child[0]->getLogLikelihoodOnePop(false));
   Log->logmsg(true,"\n");

}
double IndividualCollection::getSumLogTheta(int i){
  return SumLogTheta[i];
}
double *IndividualCollection::getSumLogTheta(){
  return SumLogTheta;
}
double IndividualCollection::getLogLikelihood(AdmixOptions *options, Chromosome **C){
  double LogL = 0.0;
  for(unsigned i = 0; i < NumInd; ++i){
    int prev = i-1;
    if(i==0)prev = NumInd-1;
    LogL += _child[i]->getLogLikelihood(options, C);
    if(NumInd > 1)_child[prev]->HMMIsBad(false);
  }

  return LogL;
}

//returns Derivative of Inverse Link Function for individual i
double IndividualCollection::DerivativeInverseLinkFunction(int i){
  double DInvLink = 1.0;
  
  if(OutcomeType){//in case no regression model
    if(OutcomeType[0] == Binary)
      DInvLink = ExpectedY[0][i] * (1.0 - ExpectedY[0][i]);
    else if(OutcomeType[0] == Continuous)DInvLink = 1.0;
  }
  return DInvLink;    
}
