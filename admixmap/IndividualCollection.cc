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
   
}

IndividualCollection::~IndividualCollection()
{
  Individual::DeleteStaticMembers();
  for(unsigned int i = 0; i < NumInd; i++){
    delete _child[i];
  }
  delete[] _child;
  delete indadmixoutput;
  delete[] Outcome;
  delete[] OutcomeType;
  free_matrix(ExpectedY, NumOutcomes);
  delete[] SumLogTheta;
}

IndividualCollection::IndividualCollection(AdmixOptions* options,InputData *Data, Genome& Loci, Chromosome **chrm)
{
  Matrix_d nullMatrix(1,1);
  OutcomeType = 0;
  NumOutcomes = 0;
  indadmixoutput = 0;
  TargetLabels = 0;
  LogLikelihood=0.0;
  SumLogLikelihood = 0.0;
  Covariates = nullMatrix;
  Outcome = 0;
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

//not needed
int *IndividualCollection::GetSumXi()
{
  return Individual::getSumXi();
}

//not needed
int IndividualCollection::GetSumXi(int j){
  return Individual::getSumXi(j);
}

//not needed 
double IndividualCollection::GetSumrho0()
{
 //   double Sumrho0 = 0;
//    for( unsigned int i = 0; i < NumInd; i++ )
//       Sumrho0 += (*_child[i]).getSumrho0();
//    return Sumrho0;
  return Individual::getSumrho0();
 }

double IndividualCollection::GetSumrho()
{
   double Sumrho = 0;
   for( unsigned int i = 0; i < NumInd; i++ )
      Sumrho += (*_child[i]).getSumrho();
   return Sumrho;
}

Matrix_d IndividualCollection::getOutcome(int j){
  return Outcome[j];
}

double IndividualCollection::getOutcome(int j, int ind){
  if(Outcome){
    return Outcome[j](ind, 0);
  }
  else return 0.0;
}

Vector_d IndividualCollection::getTargetCol(int j, int k){
  return Outcome[j].GetColumn(k);
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

int IndividualCollection::GetNumberOfInputRows(){
  return Input.GetNumberOfRows();
}
int IndividualCollection::GetNumberOfInputCols(){
  return Input.GetNumberOfCols();
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
  return TargetLabels[k];
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

void IndividualCollection::Initialise(AdmixOptions *options, double **beta, Genome *Loci, std::string *PopulationLabels,
				      double rhoalpha, double rhobeta, LogWriter *Log, const Matrix_d &MLEMatrix){
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
      
    if( Input.GetNumberOfRows() == (int)NumInd ){
      Covariates = ConcatenateHorizontally( temporary, Input );
      Vector_d mean;
      mean = Input.ColumnMean();
      for(unsigned int i = 0; i < NumInd; i++ )
	for( int j = 0; j < Input.GetNumberOfCols(); j++ )
	  Input( i, j ) -= mean(j);
    } else {
      Covariates = temporary;
    }

    
    if( !options->getTestForAdmixtureAssociation() && options->getPopulations() > 1 ){
      temporary.SetNumberOfElements( NumInd, options->getPopulations() - 1 );
      temporary.SetElements( 1 / options->getPopulations() );
      Covariates = ConcatenateHorizontally( Covariates, temporary );
    }

    //should be in Regression
    for( int k = 0; k < NumOutcomes; k++ ){
      SetExpectedY(k,beta[k]);
//for logistic regression
      if( OutcomeType[k] )calculateExpectedY(k);
    }
  }

  SumLogTheta = new double[ options->getPopulations()];
  InitialiseMLEs(rhoalpha,rhobeta,options, MLEMatrix);
  //set to very large negative value (effectively -Inf) so the first value is guaranteed to be greater
  MaxLogLikelihood.assign(NumInd, -9999999 );
}

// ** this function needs debugging
void IndividualCollection::InitialiseMLEs(double rhoalpha, double rhobeta, AdmixOptions * options, const Matrix_d &MLEMatrix){
  //set thetahat and rhohat, estimates of individual admixture and sumintensities
  //NB NumInd = 1 here

  thetahat = new double *[NumInd];
  thetahatX = new double *[NumInd];

   Matrix_d temp;
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
      rhohat[0][0] = MLEMatrix( options->getPopulations(), 0 );
      if( options->getXOnlyAnalysis() ){
	temp = MLEMatrix.SubMatrix( 0, options->getPopulations() - 1, 0, 0 );
	for(int k = 0; k < options->getPopulations(); ++k) thetahat[0][k] = temp(k,0);
      }
      else{
	temp = MLEMatrix.SubMatrix( 0, options->getPopulations() - 1, 0, 1 );
	for(int k = 0; k < options->getPopulations(); ++k) {
	  thetahat[0][k] = temp(k,0);
	  thetahat[0][k+ options->getPopulations()] = temp(k,1);
	}
	rhohat[0][1] = MLEMatrix(options->getPopulations(), 1 );
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
  }
  if ( Input.GetNumberOfMissingValues() ) Input.SetMissingValuesToColumnMeans();

  if ( strlen( options->getOutcomeVarFilename() ) != 0 ){
    LoadOutcomeVar(options, data_, Log);
  }
  if ( strlen( options->getReportedAncestryFilename() ) != 0 ){
    LoadRepAncestry(data_);
  }
}

void IndividualCollection::LoadCovariates(InputData *data_){
  //LOAD INPUT (COVARIATES)

  Input = data_->getInputMatrix();

  Input.SubMatrix2( 1, NumInd, 0, Input.GetNumberOfCols() - 1 );
  CovariateLabels = new string[ Input.GetNumberOfCols() ];
  Vector_i vtemp( Input.GetNumberOfCols() );
  vtemp.SetElements( 1 );
  getLabels(data_->getInputData()[0], vtemp, CovariateLabels);//?getCovariatelabels()
}

void IndividualCollection::LoadOutcomeVar(AdmixOptions *options, InputData *data_, LogWriter *Log){
  //LOAD TARGET (OUTCOME VARIABLE)
  Matrix_d TempTarget, temporary;
  string *TempLabels = 0;

  //conversion necessary because LoadTarget is changed further down
  Matrix_d& LoadTarget = (Matrix_d&)data_->getTargetMatrix();
  TempLabels = new string[ LoadTarget.GetNumberOfCols() ];

  Vector_i vtemp( LoadTarget.GetNumberOfCols() );
  vtemp.SetElements(1);
  getLabels(data_->getTargetData()[0], vtemp, TempLabels);

  if( options->getAnalysisTypeIndicator() == 5 ){
    TargetLabels = new string[ LoadTarget.GetNumberOfCols() ];
    NumOutcomes = LoadTarget.GetNumberOfCols();
    Outcome = new Matrix_d[NumOutcomes];
    delete[] OutcomeType;
    OutcomeType = new int[ LoadTarget.GetNumberOfCols() ];

    for( int j = 0; j < LoadTarget.GetNumberOfCols(); j++ ){
      TargetLabels[j] = TempLabels[j];
      TempTarget = LoadTarget;
      TempTarget.SubMatrix2( 1, NumInd, j, j );
      for(unsigned int i = 0; i < NumInd; i++ ){
	if( !TempTarget.IsMissingValue( i, 0 ) &&
	    (TempTarget( i, 0 ) == 0 || TempTarget( i, 0 ) == 1) )//binary outcome
	  OutcomeType[j] = 1;// 1 => binary outcome
	else OutcomeType[j] = 0;// 0 => continuous outcome
      }
      Outcome[j] = TempTarget;

      if( OutcomeType[j] )
	{
	  Log->logmsg(true,"Binary variable: ");
	  Log->logmsg(true,getTargetLabels(j));
	  Log->logmsg(true,".\n");
	}
      else
	{
	  Log->logmsg(true,"Continuous variable: ");
	  Log->logmsg(true,getTargetLabels(j));
	  Log->logmsg(true,".\n");
	}
    }
  }
  else{
    TargetLabels = new string[ 1 ];
    TargetLabels[0] = TempLabels[ options->getTargetIndicator() ];
    NumOutcomes = 1;
    Outcome = new Matrix_d[1];
    Log->logmsg(true,"Regressing on: ");
    Log->logmsg(true,getTargetLabels(0));
    Log->logmsg(true,".\n");

    LoadTarget.SubMatrix2( 1, NumInd, options->getTargetIndicator(), options->getTargetIndicator() );
    Outcome[0] = LoadTarget;
  }

  delete [] TempLabels;
}

void IndividualCollection::LoadRepAncestry(InputData *data_){
  //LOAD REPORTED ANCESTRY IF GIVEN   

  ReportedAncestry = new Matrix_d[NumInd];
  const Matrix_d& temporary = data_->getReportedAncestryMatrix();
  for( int i = 0; i < temporary.GetNumberOfRows() / 2; i++ )
    ReportedAncestry[i] = temporary.SubMatrix( 2*i, 2*i + 1, 0, temporary.GetNumberOfCols() - 1 );
 
}

void IndividualCollection::getLabels( const string buffer, Vector_i temporary, string *labels )
{
  StringSplitter splitter;
  const Vector_s& labels_tmp = splitter.split(buffer);

  for (size_t i = 0, index = 0; i < labels_tmp.size(); ++i) {
    if (temporary.GetNumberOfElements() == 1 || temporary(i)) {            
      labels[index++] = labels_tmp[i];
    }
  }
}
void IndividualCollection::getLabels(const Vector_s& data, Vector_i temporary, string *labels)
{
    for (size_t i = 0, index = 0; i < data.size(); ++i) {
        if (temporary.GetNumberOfElements() == 1 || temporary(i)) {            
            labels[index++] = data[i];
        }
    }
}

void IndividualCollection::Update(int iteration, AlleleFreqs *A, Regression *R, const double *poptheta, AdmixOptions *options,
				  Chromosome **chrm, vector<vector<double> > &alpha, double rhoalpha, double rhobeta,
				  LogWriter *Log, chib *MargLikelihood){
  fill(SumLogTheta, SumLogTheta+options->getPopulations(), 0.0);//reset to 0
  if(iteration > options->getBurnIn())Individual::ResetScores(options);
  Individual::ResetStaticSums();

  for(unsigned int i = 0; i < NumInd; i++ ){
    
    if( options->getPopulations() > 1 ){
      _child[i]->SampleParameters(i, SumLogTheta, A, iteration , Outcome, NumOutcomes, OutcomeType, ExpectedY, R->getlambda(),
				  R->getNoCovariates(),  Covariates, R->getbeta(), poptheta, options, 
				  chrm, alpha, rhoalpha, rhobeta, sigma,  
				  DerivativeInverseLinkFunction(options->getAnalysisTypeIndicator(), i),
				  R->getDispersion(OutcomeType[0]));}
    
    else{
      _child[i]->OnePopulationUpdate(i, Outcome, NumOutcomes, OutcomeType, ExpectedY, R->getlambda(), options->getAnalysisTypeIndicator(), 
				     chrm, A);
    }   
    
    if( (options->getAnalysisTypeIndicator() < 0) &&  options->getMLIndicator() && (i == 0) )//check if this condition is correct
      _child[i]->ChibLikelihood(iteration, &LogLikelihood, &SumLogLikelihood, &(MaxLogLikelihood[i]),
				options, chrm, alpha, rhoalpha, rhobeta,
				thetahat[i], thetahatX[i], rhohat[i], rhohatX[i], Log, MargLikelihood, A);
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
