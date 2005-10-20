#include "Regression.h"
#include <numeric>

std::ofstream Regression::outputstream;

Regression::Regression(){
  lambda = 0; 
  lambda0 = 0.0;
  lambda1 = 0.0;
  SumLambda = 0;
  beta = 0;
  beta0 = 0;
  betan = 0;
  n0 = 0;
  SumBeta = 0;
  NumOutcomeVars = 0;
  NumCovariates = 0;
  BetaDrawArray = 0;
  acceptbeta = 0;
  BetaParameters = 0;
  RegType = None;

  X = 0;
  XtY = 0;
  dims = 0;
}

Regression::~Regression(){
  if(RegType == Logistic){
    for(int i = 0; i < NumCovariates; i++){
      delete BetaDrawArray[i];
    }
    delete[] BetaDrawArray;
    delete[] dims;
  }
  delete[] XtY;
  delete[] beta0;
  delete[] betan;
  delete[] n0;
  delete[] beta;
  delete[] SumBeta;
  delete[] BetaParameters;

}

void Regression::OpenOutputFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
				const std::string* const PopulationLabels, LogWriter *Log){
  //Open paramfile
  if ( options->getIndAdmixHierIndicator()){ 
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
      InitializeOutputFile(options, individuals, PopulationLabels);
      }
    }
    else{
      Log->logmsg(true,"No regparamfile given\n");
    }
  }
}

void Regression::InitializeOutputFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
				      const std::string* const PopulationLabels)
{
  // Header line of paramfile
  for( int kk = 0; kk < options->getNumberOfOutcomes(); kk++ ){
      outputstream << "\"intercept\" ";
      for( int i = 0; i < individuals->GetNumberOfInputCovariates(); i++ ){
	outputstream << individuals->getCovariateLabels(i) << " ";
      }
      if( !options->getTestForAdmixtureAssociation() )
	for( int k = 1; k < options->getPopulations(); k++ ){
	  outputstream << "\"slope." << PopulationLabels[k] << "\" ";
	}
      if( individuals->getOutcomeType(kk) == Continuous ){
	outputstream<< setprecision(6) << "\"precision\"";
      }
    }
  outputstream << endl;
}

void Regression::Initialise(unsigned Number, const IndividualCollection* const individuals, LogWriter *Log){
  //set regression number for this object
  RegNumber = Number;

    //determine regression type
    if( individuals->getOutcomeType(RegNumber)== Binary ) RegType = Logistic;
    else if( individuals->getOutcomeType(RegNumber)== Continuous )RegType = Linear;

    //** Objects common to both regression types
    NumCovariates = individuals->GetNumCovariates();

    beta = new double[ NumCovariates ];
    SumBeta = new double[ NumCovariates ];
    fill(beta, beta + NumCovariates, 0.0);
    fill(SumBeta, SumBeta + NumCovariates, 0.0);

    double p;
    beta0 = new double[ NumCovariates ];
    fill(beta0, beta0 + NumCovariates, 0.0);
    
    std::vector<double> v = individuals->getOutcome(RegNumber);
    p = accumulate(v.begin(), v.end(), 0.0, std::plus<double>()) / (double)v.size();
    if(RegType == Logistic )
      beta0[0] = log( p / ( 1 - p ) );
    else if(RegType == Linear)
      beta0[0] = p;
      
    for(int j = 0; j < NumCovariates; ++j){
      beta[j] = beta0[j];
    }

    X = individuals->getCovariates();
    XtY = new double[NumCovariates];

    // if linear regression, n0*lambda is the prior precision matrix of regression coefficients given lambda   
    // if logistic regression, lambda is the prior precision for regression coefficients 
    lambda = 0.01; 
    SumLambda = lambda;

    // ** Initialise Linear Regression objects    
    if(RegType == Linear){
      betan = new double[NumCovariates];
      n0 = new double[ NumCovariates * NumCovariates ];
      fill(n0, n0+NumCovariates*NumCovariates, 0.0);
      for(int i = 0; i < NumCovariates; ++i)n0[i*NumCovariates +i] = 1.0; // identity matrix
    
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

    dims = new int[2];
    dims[0] = individuals->getSize();
    dims[1] = NumCovariates;
    fill(BetaParameters, BetaParameters+NumCovariates+4, 0.0);
    BetaParameters[ NumCovariates + 1 ] = lambda;
    BetaParameters[ NumCovariates + 3 ] = beta0[0];
    for( int i = 0; i < NumCovariates; i++ ){
      BetaDrawArray[i] = new GaussianProposalMH( BetaParameters, lr, dlr, ddlr, dims, X );
    }
  }
}

void Regression::SetExpectedY(IndividualCollection *IC)const{
  IC->SetExpectedY(RegNumber, beta);
  if( RegType == Logistic ) IC->calculateExpectedY(RegNumber);
}

void Regression::Update(bool afterBurnIn, IndividualCollection* individuals){
  // Sample for regression model parameters beta
  // should make sure that matrix returned by getCovariates contains updated values of indiv admixture
  std::vector<double> Outcome = individuals->getOutcome(RegNumber);
  double Y[Outcome.size()];
  for(unsigned i = 0; i < Outcome.size(); ++i)
    Y[i] = Outcome[i];
  
  X = individuals->getCovariates();
  
  int NumIndividuals = individuals->getSize();
  matrix_product(Y, X, XtY, 1, NumIndividuals, NumCovariates);//XtY = Xt * Y
  
  
  if( RegType == Linear ){
    // full conditional for regression parameters has mean Beta_n = (n0 + X'X)^-1 (n0 Beta0 + X'Y)
    // where Beta0  and n0*lambda are prior mean and prior precision matrix for regression parameters
    // calculate (n0 + X'X)^-1

    //compute (X'X + n0)^-1
    double Covariance[NumCovariates*NumCovariates];

    //TODO: this block should be hidden in matrix functions
    gsl_matrix_view XX, C;
    XX = gsl_matrix_view_array(const_cast<double *>(X), NumIndividuals, NumCovariates);
    C = gsl_matrix_view_array(Covariance, NumCovariates, NumCovariates);
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, &XX.matrix, &XX.matrix, 0, &C.matrix); 
    add_matrix(Covariance, n0, NumCovariates, NumCovariates);
    matrix_inverse(Covariance, NumCovariates);

    // compute betan = temp*(n0*beta0 + X'*Y) to obtain means of full conditional distribution
    double temp2[NumCovariates];
    matrix_product(n0, beta0, temp2, NumCovariates, NumCovariates, 1);
    add_matrix(temp2, XtY, NumCovariates, 1);
    matrix_product(Covariance, temp2, betan, NumCovariates, NumCovariates, 1);

    // lambda_n is rate parameter of gamma full conditional distribution for dispersion parameter, given by  
    // lambda1 + 0.5[(Y - X Beta)' Y + (Beta0 - Beta_n)' n0 Beta0]
    double Xbetan[NumIndividuals];
    matrix_product(X, betan, Xbetan, NumIndividuals, NumCovariates, 1);//Xbetan = X * betan
    scale_matrix(Xbetan, -1.0, NumIndividuals, 1);
    add_matrix(Xbetan, Y, NumIndividuals, 1);//Xbetan = Y - X* betan
    double lambdan[1];
    matrix_product(Xbetan, Y, lambdan, 1, NumIndividuals, 1);//lambdan = Xbetan' * Y
    lambdan[0] *= 0.5;
    lambdan[0] += lambda1; //lambdan = lambda1  + 0.5 * Xbetan' * Y

    lambda = gengam( lambdan[0], lambda0 + 0.5 * individuals->getSize() );
    scale_matrix(Covariance, lambda, NumCovariates, NumCovariates);

    DrawBeta.SetMean( betan );
    DrawBeta.SetCovariance( Covariance );
    DrawBeta.Draw(beta);
  }
  
  else if( RegType == Logistic ){
   
    for( int j = 0; j < NumCovariates; j++ ){
      BetaDrawArray[j]->UpdateDoubleData( X );
      BetaParameters[ NumCovariates + 2 ] = j;
      BetaParameters[ NumCovariates ] = XtY[ j ];
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

void Regression::Output(int iteration, AdmixOptions *options, LogWriter *Log)const{
  //output to logfile
  if( !options->useCOUT() || iteration == 0 )
    {
      if( RegType != None )
	{
	  if(options->getNumberOfOutcomes()==2)Log->write("\nRegression ");Log->write((int)RegNumber);
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
	if(options->getNumberOfOutcomes()==2)cout << "\nRegression " << RegNumber << " ";
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
	if(options->getNumberOfOutcomes()< 2 || RegNumber==1)outputstream << endl;
	//output new line in paramfile when last regression model
  }
}
void Regression::OutputErgodicAvg(int samples, std::ofstream *avgstream)const{
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
double *Regression::getbeta()const{
  if(beta)//in case beta not allocated (happens if no regression model); may be unnecessary if beta initialised to 0
    return beta;
  else return NULL;
}
double Regression::getlambda()const{
  return lambda;
}
int Regression::getNumCovariates()const{
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

double Regression::getDispersion()const{
  //returns dispersion parameter
  double dispersion = 1.0;
  if( RegType == Linear ) dispersion = lambda;//linear regression
  else if(RegType == Logistic) dispersion = 1.0;
  return dispersion;
}

double Regression::getLogLikelihood(IndividualCollection *IC)const{
  int NumIndividuals = IC->getSize();
  double dev[NumIndividuals];
  double devsq[1];
  for(int i = 0; i < NumIndividuals; ++i)
    dev[i] = IC->getOutcome(RegNumber, i) - IC->getExpectedY(i);
  matrix_product(dev, dev, devsq, 1, NumIndividuals, 1);
  return -0.5* (NumIndividuals * NumCovariates*log(2.0*3.14159) - NumIndividuals*log(lambda) + lambda*devsq[0]);

}
double Regression::getLogLikelihoodAtPosteriorMeans(IndividualCollection *IC, int iterations){
  int NumIndividuals = IC->getSize();
  double logL = 0.0;

  //set expected outcome at posterior means of regression parameters
  for(int i = 0; i < NumCovariates; ++i)SumBeta[i] /= (double)iterations; 
  IC->SetExpectedY(RegNumber,SumBeta);//computes X * BetaBar
  for(int i = 0; i < NumCovariates; ++i)SumBeta[i] *= (double)iterations; 

  double Y;

  if(RegType ==Linear){
    //multivariate Gaussian density
    double dev[NumIndividuals];
    double devsq[1];
    for(int i = 0; i < NumIndividuals; ++i){
      Y = IC->getOutcome(RegNumber, i);
      dev[i] = Y - IC->getExpectedY(i);
    }
    matrix_product(dev, dev, devsq, 1, NumIndividuals, 1);//precision has form lambda * I so
    //quadratic form in density is lambda * dev' * dev
    logL = -0.5* (NumIndividuals * NumCovariates*log(2.0*3.14159) - NumIndividuals*log(lambda) + lambda*devsq[0]);
  }
  else if(RegType == Logistic){
    //loglikelihood is sum of logs of bernoulli probabilities, given by EY
    IC->calculateExpectedY(RegNumber);
    double pi;
    for(int i = 0; i < NumIndividuals; ++i){
      Y = IC->getOutcome(RegNumber, i);
      pi = IC->getExpectedY(i);
      if(Y == 1) logL += log(pi);
      else if(Y == 0)logL += log(1.0-pi);
    }
  }

  return logL;
}
