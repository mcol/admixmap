#include "Regression.h"
#include "matrix_d.h"

std::ofstream Regression::outputstream;

Regression::Regression(){
  lambda = 0; 
  lambda0 = 0.0;
  lambda1 = 0.0;
  SumLambda = 0;
  beta = 0;
  SumBeta = 0;
  NumOutcomeVars = 0;
  NumCovariates = 0;
  BetaDrawArray = 0;
  acceptbeta = 0;
  BetaParameters = 0;
  RegType = None;

  aCovariates = 0;
  dims = 0;
}

Regression::~Regression(){
  if(RegType == Logistic){
    for(int i = 0; i < NumCovariates; i++){
      delete BetaDrawArray[i];
    }
    delete[] BetaDrawArray;
    delete[] dims;
    delete[] aCovariates;
  }

  delete[] beta;
  delete[] SumBeta;
  delete[] BetaParameters;

}

void Regression::OpenOutputFile(AdmixOptions *options, IndividualCollection *individuals,std::string *PopulationLabels, LogWriter *Log){
  //Open paramfile 
  if ( strlen( options->getRegressionOutputFilename() ) ){
    outputstream.open( options->getRegressionOutputFilename(), ios::out );
    if( !outputstream )
      {
	Log->logmsg(true,"ERROR: Couldn't open regparamfile\n");
	exit( 1 );
      }
    else{
      Log->logmsg(true,"Writing regression parameters to ");
      Log->logmsg(true, options->getRegressionOutputFilename());
      Log->logmsg(true,"\n");
      if( options->getTextIndicator() )InitializeOutputFile(options, individuals, PopulationLabels);
    }
  }
  else{
    Log->logmsg(true,"No regparamfile given\n");
  }
}

void Regression::InitializeOutputFile(AdmixOptions *options, IndividualCollection *individuals,std::string *PopulationLabels)
{
  // Header line of paramfile
  if( options->getAnalysisTypeIndicator() == 2 || options->getAnalysisTypeIndicator() == 3 ){
    outputstream << "\"intercept\" ";
   
    if(strlen(options->getCovariatesFilename()) > 0){//if covariatesfile specified
      for( int i = 0; i < individuals->GetNumberOfInputCovariates(); i++ ){
	outputstream << individuals->getCovariateLabels(i) << " ";
      }
    } 
    if( !options->getTestForAdmixtureAssociation() ){
      for( int k = 1; k < options->getPopulations(); k++ ){
	outputstream << "\"slope." << PopulationLabels[k].substr(1) << " ";
      }
    }
    if( options->getAnalysisTypeIndicator() == 2 )//AnalysisType = 2 => linear regression
      outputstream << setprecision(6) << "\"precision\"";
  }
  else if( options->getAnalysisTypeIndicator() == 5 ){
    for( int kk = 0; kk < 2; kk++ ){
      outputstream << "\"intercept\" ";
      for( int i = 0; i < individuals->GetNumberOfInputCovariates(); i++ ){
	outputstream << individuals->getCovariateLabels(i) << " ";
      }
      for( int k = 1; k < options->getPopulations(); k++ ){
	outputstream << "\"slope." << PopulationLabels[k].substr(1) << " ";
      }
      if( ! individuals->getOutcomeType(kk) ){
	outputstream<< setprecision(6) << "\"precision\"";
      }
    }
  }
  outputstream << endl;
}

void Regression::Initialise(unsigned Number, IndividualCollection *individuals, AdmixOptions *options, 
			    LogWriter *Log){
  //set regression number for this object
  RegNumber = Number;

  //?? next 2 lines may be unecessary, not initalising unless there is a regression model
  if( options->getAnalysisTypeIndicator() < 2 ){//no regression model
    NumCovariates = 0;
    RegType = None;
  }
  else{

    //determine regression type
    if( individuals->getOutcomeType(RegNumber)==1 ) RegType = Logistic;
    else if( individuals->getOutcomeType(RegNumber)==0 )RegType = Linear;

    //** Objects common to both regression types
    NumCovariates = individuals->GetNumCovariates();

    beta = new double[ NumCovariates ];
    SumBeta = new double[ NumCovariates ];
    fill(beta, beta + NumCovariates, 0.0);
    fill(SumBeta, SumBeta + NumCovariates, 0.0);

    double p;
    beta0.SetNumberOfElements(NumCovariates, 1 );
    
    p = (individuals->getTargetCol(RegNumber,0)).Mean();
    if(RegType == Logistic )
      beta0(0,0) = log( p / ( 1 - p ) );
    else if(RegType == Linear)
      beta0(0,0) = p;
      
    for(int j = 0; j < NumCovariates; ++j){
      beta[j] = beta0(j,0);
    }

    // if linear regression, n0*lambda is the prior precision matrix of regression coefficients given lambda   
    // if logistic regression, lambda is the prior precision for regression coefficients 
    lambda = 0.01; 
    SumLambda = lambda;

    // ** Initialise Linear Regression objects    
    if(RegType == Linear){
       
      n0.SetNumberOfElements( NumCovariates, NumCovariates );
      n0.SetDiagonal( 1 ); // identity matrix
    
      lambda0 = 0.01;//shape parameter for prior on lambda
      lambda1 = 0.01;//rate parameter for prior on lambda
      
      Log->logmsg(true,"\nNormal-inverse-gamma prior for linear regression model with gamma shape parameter ");
      Log->logmsg(true, lambda0);
      Log->logmsg(true, " and rate parameter "); Log->logmsg(true, lambda1); Log->logmsg(true, "\n");
      
      DrawBeta.SetDimension( NumCovariates );
    }

    // ** Initialise Logistic Regression objects
  else if( RegType == Logistic) {
    Log->logmsg(true,"\nGaussian priors on logistic regression parameters with precision ");
    Log->logmsg(true, lambda); Log->logmsg(true, "\n");
    
    //  ** initialize sampler for logistic regression **
    acceptbeta = 0;
    BetaParameters = new double[ NumCovariates + 4 ]; // array elements consist of 
    // one beta param for each covariate followed by 
    // precision of priors on logistic regression params for each outcome var, 
    // followed by one intercept param for each outcome var
    BetaDrawArray = new GaussianProposalMH*[NumCovariates];
    
    for( int i = 0; i < NumCovariates; i++ ){
      BetaDrawArray[i] = 0;
    }

    Matrix_d Cov = individuals->getCovariates();
    aCovariates = new double[individuals->getSize() * NumCovariates];
    for(int i = 0; i < Cov.GetNumberOfRows(); ++i)
      for(int j = 0; j < Cov.GetNumberOfCols(); ++j)
	aCovariates[i*NumCovariates +j] = Cov(i,j);

    dims = new int[2];
    dims[0] = individuals->getSize();
    dims[1] = NumCovariates;
    fill(BetaParameters, BetaParameters+NumCovariates+4, 0.0);
    BetaParameters[ NumCovariates + 1 ] = lambda;
    BetaParameters[ NumCovariates + 3 ] = beta0(0,0);
    for( int i = 0; i < NumCovariates; i++ ){
      BetaDrawArray[i] = new GaussianProposalMH( BetaParameters, lr, dlr, ddlr, dims, aCovariates );
    }
  }
  }//end else, AnalysisTypeIndicator >2
}

void Regression::SetExpectedY(IndividualCollection *IC){
  IC->SetExpectedY(RegNumber, beta);
  if( RegType == Logistic ) IC->calculateExpectedY(RegNumber);
}

void Regression::Update(bool afterBurnIn, IndividualCollection *individuals){
  // Sample for regression model parameters beta
  // should make sure that matrix returned by getCovariates contains updated values of indiv admixture

    if( RegType == Linear ){
      // full conditional for regression parameters has mean Beta_n = (n0 + X'X)^-1 (n0 Beta0 + X'Y)
      // where Beta0  and n0*lambda are prior mean and prior precision matrix for regression parameters
      // calculate (n0 + X'X)^-1
      Matrix_d temporary = individuals->getCovariates().Transpose() * individuals->getCovariates() + n0;
      temporary.InvertUsingLUDecomposition();
      temporary.Symmetrize();
      // postmultiply by (n0*beta0 + X'*Y) to obtain means of full conditional distribution
      betan = ( n0 * beta0 + individuals->getCovariates().Transpose() * individuals->getOutcome(RegNumber) );
      betan = temporary * betan;
      // lambda_n is rate parameter of gamma full conditional distribution for dispersion parameter, given by  
      // lambda1 + 0.5[(Y - X Beta)' Y + (Beta0 - Beta_n)' n0 Beta0]
      double lambdan = lambda1 + 0.5 *
	( (individuals->getOutcome(RegNumber) - individuals->getCovariates() * betan ).Transpose() * individuals->getOutcome(RegNumber)
	  + (beta0 - betan).Transpose() * n0 * beta0 )(0,0);
      lambda = gengam( lambdan, lambda0 + 0.5 * individuals->getSize() );
      DrawBeta.SetMean( betan );
      DrawBeta.SetCovariance( temporary / lambda);
      DrawBeta.Draw(beta);
    }

    else if( RegType == Logistic ){
      sum = (individuals->getOutcome(RegNumber)).Transpose() * individuals->getCovariates();

      Matrix_d Cov = individuals->getCovariates();
      for(int i = 0; i < Cov.GetNumberOfRows(); ++i)
	for(int j = 0; j < Cov.GetNumberOfCols(); ++j)
	  aCovariates[i*NumCovariates +j] = Cov(i,j);

      for( int j = 0; j < NumCovariates; j++ ){
	BetaDrawArray[j]->UpdateDoubleData( aCovariates );
	BetaParameters[ NumCovariates + 2 ] = j;
	BetaParameters[ NumCovariates ] = sum( 0, j );
	BetaDrawArray[j]->UpdateParameters( BetaParameters );
	acceptbeta = BetaDrawArray[j]->Sample( &( beta[j] ) );
	BetaParameters[j] = beta[j];
      }
    }
    individuals->SetExpectedY(RegNumber,beta);
    if( individuals->getOutcomeType(RegNumber) )
      individuals->calculateExpectedY(RegNumber);

  if(afterBurnIn)
    SumParameters();
}//end Update

void Regression::Output(int iteration, AdmixOptions *options, LogWriter *Log){
  //output to logfile
  if( !options->useCOUT() || iteration == 0 )
    {
      if( RegType != None )
	{
	  if(options->getAnalysisTypeIndicator()==5)Log->write("\nRegression ");Log->write((int)RegNumber);
          for( int j = 0; j < NumCovariates; j++ )
	    {
	      Log->width(9);
	      Log->write(beta[j],6);
	    }
          //Log->width(9);
          if( RegType == Linear )
	    {
	      Log->write(lambda,6);
	      //Log->write( lambda, 6);
	    }
	}
    }
  //output to screen
  if( options->useCOUT() )
    {
      if( RegType != None ){
	if(options->getAnalysisTypeIndicator()==5)cout << "\nRegression " << RegNumber << " ";
	for( int j = 0; j < NumCovariates; j++ ){
	  (cout).width(9);
	  cout << setprecision(6) << beta[j] << " ";
	}
	(cout).width(9);
	if( RegType == Linear )
	  cout << setprecision(6)
	       << lambda<<" ";

      }
    }
  //Output to paramfile after BurnIn
  if( iteration > options->getBurnIn() ){
	if( RegType != None ){
	  for( int j = 0; j < NumCovariates; j++ ){
	    outputstream.width(9);
	    outputstream << setprecision(6) << beta[j] << " ";
	  }
	  outputstream.width(9);
     if( RegType == Linear )
       outputstream << setprecision(6) << lambda << " ";
	}
	if(options->getAnalysisTypeIndicator()<5 || RegNumber==1)outputstream << endl;
	//output new line in paramfile when last regression model
  }
}
void Regression::OutputErgodicAvg(int samples, std::ofstream *avgstream){
 //output to ergodicaveragefile
  if( RegType != None ){
    for( int j = 0; j < NumCovariates; j++ ){
      avgstream->width(9);
      *avgstream << setprecision(6) << SumBeta[j] / samples << " ";
    }
    avgstream->width(9);
    if( RegType == Linear )
      *avgstream << setprecision(6) << SumLambda / samples << " ";
  }
}

void Regression::SumParameters(){
  // accumulate sum of parameters after burnin.
  if( NumCovariates > 0 )
      for(int j = 0; j < NumCovariates; ++j)
        SumBeta[j] += beta[j];
  SumLambda += lambda;
}
double *Regression::getbeta(){
  if(beta)//in case beta not allocated (happens if no regression model); may be unnecessary if beta initialised to 0
    return beta;
  else return NULL;
}
double Regression::getlambda(){
  return lambda;
}
int Regression::getNumCovariates(){
  return NumCovariates;
}

double Regression::lr( const double* const parameters, const int* const dims, const double* const data, const double beta )
{
  //dims is an array of length 2 containing the dimensions of data
  int n = dims[0];
  int d = dims[1];
  int index = (int)parameters[ d + 2 ];
  double beta0 = 0;
  if( index == 0 )
    beta0 = parameters[ d + 3 ];
  double *Xbeta, *beta1;

  beta1 = new double[ d ];
  Xbeta = new double[ n ];

  double f = parameters[ d ] * beta
    - 0.5 * parameters[ d + 1 ] * (beta - beta0) * (beta - beta0);
  
  for( int i = 0; i < d; i++ )
    {
      if( i != index )
	beta1[ i ] = parameters[i];
      else
	beta1[ i ] = beta;
    }
  
  //Xbeta = data * beta1;
  matrix_product(data, beta1, Xbeta, n, d, 1);

  for( int i = 0; i < n; i++ ){
    f -= log( 1. + exp( Xbeta[ i ] ) );}
  
  delete[] beta1;
  delete[] Xbeta;
  return( f );
}

//a lot of duplicated code in the next 2 functions
//can we find a way to use a single function to compute Xbeta?
double Regression::dlr( const double* const parameters, const int* const dims, const double* const data, const double beta )
{
  //dims is an array of length 2 containing the dimensions of data
  int n = dims[0];
  int d = dims[1];
  int index = (int)parameters[ d + 2 ];
  double beta0 = 0;
  if( index == 0 )
    beta0 = parameters[ d + 3 ];
  double f = parameters[ d ] - parameters[ d + 1 ] * (beta - beta0);
  double *Xbeta, *beta1;

  beta1 = new double[ d ];
  Xbeta = new double[ n ];
  
  for( int i = 0; i < d; i++ )
    {
      if( i != index )
	beta1[ i ] = parameters[i];
      else
	beta1[ i ] = beta;
    }
  
  //Xbeta = data * beta1;
  matrix_product(data, beta1, Xbeta, n, d, 1);
  for( int i = 0; i < n; i++ )
    {
      f -= data[ i*d + index ] / ( 1. + exp( -Xbeta[ i ] ) );
    }
  delete[] beta1;
  delete[] Xbeta;
  return( f );
}

double Regression::ddlr( const double* const parameters, const int* const dims, const double* const data, const double beta )
{
  //dims is an array of length 2 containing the dimensions of data
  int n = dims[0];
  int d = dims[1];
  int index = (int)parameters[ d + 2 ];
  double f = -parameters[ d + 1 ];
   double *Xbeta, *beta1;

  beta1 = new double[ d ];
  Xbeta = new double[ n ];
  
  for( int i = 0; i < d; i++ )
    {
      if( i != index )
	beta1[ i ] = parameters[i];
      else
	beta1[ i ] = beta;
    }
  
  //Xbeta = data * beta1;
  matrix_product(data, beta1, Xbeta, n, d, 1);
  for( int i = 0; i < n; i++ )
    {
      f -= data[ i*d + index ] * data[ i*d + index ] / ( 2. + exp( -Xbeta[ i ] ) + exp( Xbeta[ i ] ) );
    }
  delete[] beta1;
  delete[] Xbeta;
  return( f );
}

double Regression::getDispersion(){
  //returns dispersion parameter
  double dispersion = 1.0;
  if( RegType == Linear ) dispersion = lambda;//linear regression
  else if(RegType == Logistic) dispersion = 1.0;
  return dispersion;
}

