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

IndividualCollection::IndividualCollection()
{
  SumLogTheta = 0;
  OutcomeType = 0;
  CovariateLabels = 0;
  NumCovariates = 0;
  ReportedAncestry = 0;
}

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

IndividualCollection::IndividualCollection(AdmixOptions* options,InputData *Data, Genome& Loci, Chromosome **chrm)
{
  Matrix_d nullMatrix(1,1);
  OutcomeType = 0;
  NumOutcomes = 0;
  NumCovariates = 0;
  indadmixoutput = 0;
  OutcomeVarLabels = 0;
  LogLikelihood=0.0;
  SumLogLikelihood = 0.0;
  Covariates = nullMatrix;
  CovariateLabels = 0;
  ExpectedY = 0;
  SumLogTheta = 0;
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

  OutcomeType = new int[1];
  OutcomeType[0] = 0; 
  if (options->getAnalysisTypeIndicator() == 3 || options->getAnalysisTypeIndicator() == 4)
    {
      OutcomeType[0] = 1;
    }
 
}

void
IndividualCollection::OutputIndAdmixture()
{
  indadmixoutput->visitIndividualCollection(*this);
  for(unsigned int i=0; i<NumInd; i++){
    indadmixoutput->visitIndividual(*_child[i], _locusfortest, LogLikelihood);
  }
}

int
IndividualCollection::getSize()
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

int IndividualCollection::getOutcomeType(int i){
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

void IndividualCollection::setAdmixtureProps(double *a, size_t size)
{
  for(unsigned int i=0; i<NumInd; i++){
    _child[i]->setAdmixtureProps(a, size);
  }
}

void IndividualCollection::setAdmixturePropsX(double *a, size_t size)
{
  for(unsigned int i=0; i<NumInd; i++){
    _child[i]->setAdmixturePropsX(a, size);
  }
}

int IndividualCollection::GetNumberOfInputCovariates(){
  return InputCovariates.GetNumberOfCols();
}
int IndividualCollection::GetNumCovariates() const{
  return NumCovariates;
}

Matrix_d IndividualCollection::getCovariates(){
  return Covariates;
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

std::string IndividualCollection::getTargetLabels(int k){
  return OutcomeVarLabels[k];
}

void IndividualCollection::SetExpectedY(int k, double *beta){
  if(ExpectedY){
    Matrix_d Beta(Covariates.GetNumberOfCols(), 1);
    Beta.SetElements(0.0);
    for(int i = 0; i < Covariates.GetNumberOfCols(); ++i)Beta(i, 0) = beta[i];
      Beta = Covariates * Beta;  //possible memory leak
      for(int j = 0; j < Covariates.GetNumberOfRows(); ++j)
        ExpectedY[k][j] = Beta(j,0);
    }
}

void IndividualCollection::calculateExpectedY( int k)
{
  if(ExpectedY)
    for(unsigned int i = 0; i < NumInd; i++ )
      ExpectedY[k][i] = 1 / ( 1 + exp( -ExpectedY[k][i] ) );
}

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
   else  if( options->getAnalysisTypeIndicator() < 0 ){
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
   if(options->getAnalysisTypeIndicator() >=2){
    ExpectedY = alloc2D_d(NumOutcomes, NumInd);
    //Covariates.SetNumberOfElements(1);
    Matrix_d temporary( NumInd, 1 );
    temporary.SetElements(1);
      
    if(strlen(options->getCovariatesFilename()) > 0){//if covariatesfile specified
      Covariates = ConcatenateHorizontally( temporary, InputCovariates );
      Vector_d mean;
      mean = InputCovariates.ColumnMean();
      for(unsigned int i = 0; i < NumInd; i++ )
	for( int j = 0; j < InputCovariates.GetNumberOfCols(); j++ )
	  InputCovariates( i, j ) -= mean(j);
    } else {
      Covariates = temporary;
    }

    
    if( !options->getTestForAdmixtureAssociation() && options->getPopulations() > 1 ){
      temporary.SetNumberOfElements( NumInd, options->getPopulations() - 1 );
      temporary.SetElements( 1 / options->getPopulations() );
      Covariates = ConcatenateHorizontally( Covariates, temporary );
    }

   }

  SumLogTheta = new double[ options->getPopulations()];
  InitialiseMLEs(rhoalpha,rhobeta,options, MLEMatrix);
  //set to very large negative value (effectively -Inf) so the first value is guaranteed to be greater
  MaxLogLikelihood.assign(NumInd, -9999999 );
}

// ** this function needs debugging
void IndividualCollection::InitialiseMLEs(double rhoalpha, double rhobeta, AdmixOptions * options, const DataMatrix &MLEMatrix){
  //set thetahat and rhohat, estimates of individual admixture and sumintensities
  //NB NumInd = 1 here

  thetahat = new double *[NumInd];
  thetahatX = new double *[NumInd];

   size_t size_admix;
   int K = options->getPopulations();
   if( options->isRandomMatingModel() )
     size_admix = K*2; // double the size for 2 gametes in RMM
   else//assortative mating
     size_admix = K;

   for(unsigned int i = 0; i < NumInd; ++i){
   thetahat[i] = new double[size_admix];
   thetahatX[i] = new double[size_admix];
   }

   vector<double> r(2, rhoalpha/rhobeta );
   rhohat.resize( NumInd, r );
   rhohatX.resize( NumInd, r );

   //use previously read values from file, if available
   if( options->getAnalysisTypeIndicator() == -2 ){
      rhohat[0][0] = MLEMatrix.get( options->getPopulations(), 0 );
      if( options->getXOnlyAnalysis() ){
	for(int k = 0; k < options->getPopulations(); ++k) thetahat[0][k] = MLEMatrix.get(k,0);
      }
      else{
	for(int k = 0; k < options->getPopulations(); ++k) {
	  thetahat[0][k] = MLEMatrix.get(k,0);
	  thetahat[0][k+ options->getPopulations()] = MLEMatrix.get(k,1);
	}
	rhohat[0][1] = MLEMatrix.get(options->getPopulations(), 1 );
      }
      setAdmixtureProps(thetahat[0], size_admix);
   }
   else if( options->getAnalysisTypeIndicator() == -1 ){
    for(unsigned int i = 0; i < NumInd; ++i){
      for(unsigned k = 0; k < size_admix; ++k)
       thetahat[i][k] = _child[i]->getAdmixtureProps()[k];
       }
   }
 
}

void IndividualCollection::LoadData(AdmixOptions *options, InputData *data_, LogWriter *Log){
   if ( strlen( options->getCovariatesFilename() ) != 0 ){
     if( options->getTextIndicator() ){  
       LoadCovariates(data_);
     }
     if( options->getTestForAdmixtureAssociation() )
       NumCovariates = InputCovariates.GetNumberOfCols() + 1;
     else
       NumCovariates = InputCovariates.GetNumberOfCols() + options->getPopulations();
   }
   else{
     if( options->getTestForAdmixtureAssociation() )
       NumCovariates = 1;
     else
       NumCovariates = options->getPopulations();
   }

  if ( InputCovariates.GetNumberOfMissingValues() ) InputCovariates.SetMissingValuesToColumnMeans();

  if ( strlen( options->getOutcomeVarFilename() ) != 0 ){
    LoadOutcomeVar(options, data_, Log);
  }
  if ( strlen( options->getReportedAncestryFilename() ) != 0 ){
    LoadRepAncestry(data_);
  }
}

void IndividualCollection::LoadCovariates(InputData *data_){
  //LOAD INPUT (COVARIATES)

  InputCovariates = data_->getInputMatrix();

  InputCovariates.SubMatrix2( 1, NumInd, 0, InputCovariates.GetNumberOfCols() - 1 );
  CovariateLabels = new string[ InputCovariates.GetNumberOfCols() ];
  //Vector_i vtemp( InputCovariates.GetNumberOfCols() );
  //vtemp.SetElements( 1 );
  getLabels(data_->getInputData()[0], CovariateLabels);
}

void IndividualCollection::LoadOutcomeVar(AdmixOptions *options, InputData *data_, LogWriter *Log){
  //LOAD TARGET (OUTCOME VARIABLE)
  DataMatrix TempTarget, temporary;
  string *TempLabels = 0;

  //conversion necessary because LoadTarget is changed further down
  DataMatrix& OutcomeVarData = (DataMatrix&)data_->getOutcomeVarMatrix();
  TempLabels = new string[ OutcomeVarData.nCols() ];

  //Vector_i vtemp( LoadTarget.GetNumberOfCols() );
  //vtemp.SetElements(1);
  getLabels(data_->getOutcomeVarData()[0], TempLabels);

  if( options->getAnalysisTypeIndicator() == 5 ){
    OutcomeVarLabels = new string[ OutcomeVarData.nCols() ];
    NumOutcomes = OutcomeVarData.nCols();
    Outcome.setDimensions(NumOutcomes, 1);
    delete[] OutcomeType;
    OutcomeType = new int[ OutcomeVarData.nCols() ];
    Outcome = OutcomeVarData.SubMatrix( 1, NumInd, 0, OutcomeVarData.nCols()-1 );

    for( unsigned j = 0; j < OutcomeVarData.nCols(); j++ ){
      OutcomeVarLabels[j] = TempLabels[j];
      TempTarget = OutcomeVarData.SubMatrix( 1, NumInd, j, j );//jth col of outcomevarfile, excluding header

      for(unsigned int i = 0; i < NumInd; i++ ){
	if( !TempTarget.isMissing( i, 0 ) &&
	    (TempTarget.get( i, 0 ) == 0 || TempTarget.get( i, 0 ) == 1) )//binary outcome
	  OutcomeType[j] = 1;// 1 => binary outcome
	else OutcomeType[j] = 0;// 0 => continuous outcome
      }
      //Outcome[j] = TempTarget;

      if( OutcomeType[j] )
	{
	  Log->logmsg(true,"Binary variable: ");
	  Log->logmsg(true,OutcomeVarLabels[j]);
	  Log->logmsg(true,".\n");
	}
      else
	{
	  Log->logmsg(true,"Continuous variable: ");
	  Log->logmsg(true,OutcomeVarLabels[j]);
	  Log->logmsg(true,".\n");
	}
    }
  }//end analysistypeindicator == 5
  else{
    OutcomeVarLabels = new string[ 1 ];
    OutcomeVarLabels[0] = TempLabels[ options->getTargetIndicator() ];
    NumOutcomes = 1;
    Log->logmsg(true,"Regressing on: ");
    Log->logmsg(true, OutcomeVarLabels[0]);
    Log->logmsg(true,".\n");
    Outcome = OutcomeVarData.SubMatrix( 1, NumInd, options->getTargetIndicator(), options->getTargetIndicator() );
    //Outcome[0] = OutcomeVarData.SubMatrix( 1, NumInd, options->getTargetIndicator(), options->getTargetIndicator() );
  }

  delete [] TempLabels;
}

void IndividualCollection::LoadRepAncestry(InputData *data_){
  //LOAD REPORTED ANCESTRY IF GIVEN   

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

void IndividualCollection::Update(int iteration, AlleleFreqs *A, Regression *R0, Regression *R1, const double *poptheta, 
				  AdmixOptions *options, Chromosome **chrm, vector<vector<double> > &alpha, 
				  double rhoalpha, double rhobeta, LogWriter *Log, chib *MargLikelihood){
  fill(SumLogTheta, SumLogTheta+options->getPopulations(), 0.0);//reset to 0
  if(iteration > options->getBurnIn())Individual::ResetScores(options);

  double lambda[] = {R0->getlambda(), R1->getlambda()};
  double *beta[] = {R0->getbeta(), R1->getbeta()};
 
  for(unsigned int i = 0; i < NumInd; i++ ){
    
    if( options->getPopulations() > 1 ){
      _child[i]->SampleParameters(i, SumLogTheta, A, iteration , &Outcome, NumOutcomes, OutcomeType, ExpectedY,
				  lambda, NumCovariates, Covariates, beta, poptheta, options, 
				  chrm, alpha, rhoalpha, rhobeta, sigma,  
				  DerivativeInverseLinkFunction(options->getAnalysisTypeIndicator(), i),
				  R0->getDispersion()
				  );
      if((iteration %2))//conjugate update of theta
	_child[i]->SampleTheta(i, iteration, SumLogTheta, &Outcome, chrm, NumOutcomes, OutcomeType, ExpectedY, lambda, NumCovariates,
			       Covariates, beta, poptheta, options, alpha, sigma,
			       DerivativeInverseLinkFunction(options->getAnalysisTypeIndicator(), i), 
			       R0->getDispersion(), false);

    }
    //?? possible error, only using dispersion parameter for first regression model
    
    else{//single population 
      _child[i]->OnePopulationUpdate(i, &Outcome, NumOutcomes, OutcomeType, ExpectedY, lambda,
				     options->getAnalysisTypeIndicator(), chrm, A);
    }   

    //calculate log posterior if necessary 
    if( options->getMLIndicator() && i == 0 && iteration > options->getBurnIn() )
      _child[i]->CalculateLogPosterior(options, alpha, rhoalpha, rhobeta);
    
    if( (options->getAnalysisTypeIndicator() < 0) &&  options->getMLIndicator() && (i == 0) )//check if this condition is correct
      _child[i]->ChibLikelihood(iteration, &LogLikelihood, &SumLogLikelihood, &(MaxLogLikelihood[i]),
				options, chrm, alpha, rhoalpha, rhobeta,
				thetahat[i], thetahatX[i], rhohat[i], rhohatX[i], Log, MargLikelihood, A);
  }

}

void IndividualCollection::ConjugateUpdateIndAdmixture(int iteration, Regression *R0, Regression *R1, const double *poptheta, 
						       AdmixOptions *options, Chromosome **chrm, vector<vector<double> > &alpha){
  if( options->getPopulations() > 1 ){
    double lambda[] = {R0->getlambda(), R1->getlambda()};
    double *beta[] = {R0->getbeta(), R1->getbeta()};

    for(unsigned int i = 0; i < NumInd; i++ )
      _child[i]->SampleTheta(i, iteration, SumLogTheta, &Outcome, chrm, NumOutcomes, OutcomeType, ExpectedY, lambda, NumCovariates,
			     Covariates, beta, poptheta, options, alpha, sigma,
			     DerivativeInverseLinkFunction(options->getAnalysisTypeIndicator(), i), 
			     R0->getDispersion(), false);
  }
}

void IndividualCollection::OutputChibEstimates(LogWriter *Log, int Populations){
  //Used only if marglikelihood = 1
  Log->write("Estimates used in Chib algorithm to estimate marginal likelihood:\n");
  for(unsigned  int i = 0; i < NumInd; i++ ){
    for(int k = 0; k < Populations; ++k)
      Log->write(thetahat[i][k]);
    for(int k = 0; k < Populations; ++k)
      Log->write(thetahat[i][Populations +k]);
    Log->write( rhohat[i][0]);Log->write( rhohat[i][1]);Log->write("\n");
  }
}

void IndividualCollection::OutputErgodicAvg(int samples, chib *MargLikelihood, std::ofstream *avgstream){
     *avgstream << SumLogLikelihood / samples << " "
               << MargLikelihood->getLogPosterior();
}

void
IndividualCollection::getOnePopOneIndLogLikelihood(LogWriter *Log, std::string *PopulationLabels)
{
   Log->logmsg(true,"Log-likelihood for unadmixed ");
   Log->logmsg(true, (*PopulationLabels)[0]);
   Log->logmsg(true, ": ");
   Log->logmsg(true, _child[0]->getLogLikelihoodOnePop());
   Log->logmsg(true,"\n");

}
double IndividualCollection::getSumLogTheta(int i){
  return SumLogTheta[i];
}
double *IndividualCollection::getSumLogTheta(){
  return SumLogTheta;
}
double IndividualCollection::getLL(){
  //not currently used
  return SumLogLikelihood;
}

//returns Derivative of Inverse Link Function for individual i
double IndividualCollection::DerivativeInverseLinkFunction(int AnalysisType,int i){
  double DInvLink = 1.0;
  int mOutcomeType = OutcomeType[0];

    //Linear regression
    if(AnalysisType == 2 ){
      DInvLink = 1.0;
      }
    //Logistic Regression
    else if( AnalysisType == 3 || AnalysisType == 4 ){
      DInvLink = ExpectedY[0][i] * (1.0 - ExpectedY[0][i]);
    }
    else if( AnalysisType == 5 ){
      DInvLink = mOutcomeType ? ExpectedY[0][i] * (1.0 - ExpectedY[0][i]) : 1.0;
    }
 
  return DInvLink;    
}
