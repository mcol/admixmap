/** 
 *   ADMIXMAP
 *   Regression.cc 
 *   Class to represent and update parameters of a regression model as used in admixmap
 *   Copyright (c) 2002-2006 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "Regression.h"
#include <numeric>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

std::ofstream Regression::outputstream;

Regression::Regression(){
  lambda = 0; 
  lambda0 = 0.0;
  lambda1 = 0.0;
  SumLambda = 0;
  beta = 0;
  betamean = 0;
  betaprecision = 0;
  R = 0;
  QY = 0;
  QX = 0;
  V = 0;
  betahat = 0;

  SumBeta = 0;
  NumOutcomeVars = 0;
  NumCovariates = 0;
  BetaDrawArray = 0;
  acceptbeta = 0;
  RegType = None;

  X = 0;
  Y = 0;
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
  delete[] beta;
  delete[] SumBeta;

  delete[] betamean;
  delete[] betaprecision;
  delete[] R;
  delete[] QY;
  delete[] QX;
  delete[] V;
  delete[] betahat;
}

void Regression::OpenOutputFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
				const std::string* const PopulationLabels, LogWriter &Log){
  //Open paramfile
  if ( options->getIndAdmixHierIndicator() && options->getNumberOfOutcomes()>0){ 
    if ( strlen( options->getRegressionOutputFilename() ) ){
      outputstream.open( options->getRegressionOutputFilename(), ios::out );
      if( !outputstream )
	{
	  Log.setDisplayMode(On);
	  Log << "ERROR: Couldn't open regparamfile\n";
	  exit( 1 );
	}
      else{
	Log.setDisplayMode(Quiet);
	Log << "Writing regression parameters to " << options->getRegressionOutputFilename() << "\n";
	InitializeOutputFile(options, individuals, PopulationLabels);
      }
    }
    else{
      Log << "No regparamfile given\n";
    }
  }
}

void Regression::InitializeOutputFile(const AdmixOptions* const options, const IndividualCollection* const individuals, 
				      const std::string* const PopulationLabels)
{
  // Header line of paramfile
  for( int kk = 0; kk < options->getNumberOfOutcomes(); kk++ ){
      outputstream << "intercept\t";
      const Vector_s& labels = individuals->getCovariateLabels();
      for( unsigned i = 0; i < labels.size(); i++ ){
	outputstream << labels[i] << "\t";
      }
      if( !options->getTestForAdmixtureAssociation() && !options->getHapMixModelIndicator() )
	for( int k = 1; k < options->getPopulations(); k++ ){
	  outputstream << "slope." << PopulationLabels[k] << "\t";
	}
      if( individuals->getOutcomeType(kk) == Continuous ){
	outputstream<< setprecision(6) << "precision\t";
      }
    }
  outputstream << endl;
}

void Regression::Initialise(unsigned Number, const IndividualCollection* const individuals, LogWriter &Log){
  Log.setDisplayMode(Quiet);
  //set regression number for this object
  RegNumber = Number;

    //determine regression type
    if( individuals->getOutcomeType(RegNumber)== Binary ) RegType = Logistic;
    else if( individuals->getOutcomeType(RegNumber)== Continuous )RegType = Linear;

    //** Objects common to both regression types
    NumCovariates = individuals->GetNumCovariates();
    NumIndividuals = individuals->getSize();

    beta = new double[ NumCovariates ];
    SumBeta = new double[ NumCovariates ];
    fill(beta, beta + NumCovariates, 0.0);
    fill(SumBeta, SumBeta + NumCovariates, 0.0);

    betamean = new double[ NumCovariates ];
    fill(betamean, betamean + NumCovariates, 0.0);
    
    std::vector<double> v = individuals->getOutcome(RegNumber);
    double p = accumulate(v.begin(), v.end(), 0.0, std::plus<double>()) / (double)v.size();
    if(RegType == Logistic )
      betamean[0] = log( p / ( 1 - p ) );
    //     else if(RegType == Linear)
    //       betamean[0] = p;
    
    //initialise regression params at prior mean
    for(int j = 0; j < NumCovariates; ++j){
      beta[j] = betamean[j];
    }
    
    X = individuals->getCovariates();
    XtY = new double[NumCovariates];

    // if linear regression, n0*lambda is the prior precision matrix of regression coefficients given lambda   
    // if logistic regression, lambda is the prior precision for regression coefficients 
    lambda = 1.0; 
    SumLambda = lambda;

    // ** Initialise Linear Regression objects    
    if(RegType == Linear){
      betaprecision = new double[NumCovariates];//prior precision for beta (prior covariances are zero)
    
      lambda0 = 0.01;//shape parameter for prior on lambda
      lambda1 = 0.01;//rate parameter for prior on lambda
      lambda = lambda1 / lambda0;//initialise to prior mean

      betaprecision = new double[NumCovariates];
      fill(betaprecision, betaprecision + NumCovariates, 0.0001);

      Log << "\nGaussian priors on linear regression parameters with zero mean and precision "<< betaprecision[0] << "\n";
      Log << "Gamma("<< lambda0 << ", " << lambda1 << ") prior on data precision.\n";
      
      betahat = new double[NumCovariates];
      V = new double[NumCovariates*NumCovariates];
      QX = new double[(NumIndividuals+NumCovariates)*NumCovariates];
      QY = new double[NumIndividuals+NumCovariates];
      R = new double[NumCovariates*NumCovariates];
      // Augment data matrices
      for(int i = 0; i < NumCovariates; ++i){
	QY[i+NumIndividuals] = betamean[i] * sqrt(betaprecision[i]);//append prior means to Y
	QX[(i+NumIndividuals)*NumCovariates + i] = sqrt(betaprecision[i]);//prior precision on beta
	for(int j = i+1; j < NumCovariates; ++j)
	  QX[(i+NumIndividuals)*NumCovariates + j] = QX[(j+NumIndividuals)*NumCovariates + i] = 0.0;
      }

      DrawBeta.SetDimension( NumCovariates );
    }

    // ** Initialise Logistic Regression objects
  else if( RegType == Logistic) {
    Log << "\nGaussian priors on logistic regression parameters with precision " << lambda << "\n";
    
    //  ** initialize sampler for logistic regression **
    acceptbeta = 0;
    BetaDrawArray = new GaussianProposalMH*[NumCovariates];
    
    for( int i = 0; i < NumCovariates; i++ ){
      BetaDrawArray[i] = 0;
    }

    dims = new int[2];
    BetaParameters.n = NumIndividuals;
    BetaParameters.d = NumCovariates;
    BetaParameters.lambda = lambda;
    BetaParameters.beta0 = betamean[0];
    for( int i = 0; i < NumCovariates; i++ ){
      BetaDrawArray[i] = new GaussianProposalMH( lr, dlr, ddlr);
    }
  }
}

void Regression::SetExpectedY(IndividualCollection *IC)const{
  IC->SetExpectedY(RegNumber, beta);
}

void Regression::Update(bool sumbeta, IndividualCollection* individuals, double coolness){
  // Sample for regression model parameters beta
  //and precision in linear regression
  std::vector<double> Outcome = individuals->getOutcome(RegNumber);

  Y = &(Outcome[0]);
  X = individuals->getCovariates();
//   double sumNAm = 0.0;
//   for(int i = 0; i < NumIndividuals; ++i)
//     sumNAm += X[i*NumCovariates + 4];
//   cout<< "SumNAm ="<<sumNAm<<endl;
  
  if( RegType == Linear ){
    SampleLinearRegressionParametersWithAnnealing(Y, X, beta, &lambda, coolness);
  }
  
  else if( RegType == Logistic ){
    BetaParameters.coolness = coolness; 
    matrix_product(Y, X, XtY, 1, NumIndividuals, NumCovariates);//XtY = X' * Y
  
    for( int j = 0; j < NumCovariates; j++ ){
      BetaParameters.Covariates = X;
      BetaParameters.beta = beta;
      BetaParameters.index = j;
      BetaParameters.XtY = XtY[ j ];

      acceptbeta = BetaDrawArray[j]->Sample( &( beta[j] ), &BetaParameters );
    }
  }
  individuals->SetExpectedY(RegNumber,beta);

  if(sumbeta){
    SumParameters();
    individuals->UpdateSumResiduals();
  }
}//end Update

void Regression::OutputParams(ostream* out){
  if( RegType != None ){
    for( int j = 0; j < NumCovariates; j++ ){
      out->width(9);
      (*out) << setprecision(6) << beta[j] << "\t";
    }
    out->width(9);
    if( RegType == Linear )
      (*out) << setprecision(6) << lambda << "\t";
  }
}

void Regression::Output(int iteration, const AdmixOptions *options, LogWriter &Log){
  //output to logfile
  if( iteration == -1 )
    {
      if( RegType != None )
	{
	  Log.setDisplayMode(Off);
	  Log.setPrecision(6);
	  if(options->getNumberOfOutcomes()==2)Log <<"\nRegression " <<(int)RegNumber << ": ";
          for( int j = 0; j < NumCovariates; j++ )
	    {
	      Log << beta[j] << "\t";
	    }
          if( RegType == Linear )
	    {
	      Log << lambda;
	    }
	}
    }
  //output to screen
  if( options->getDisplayLevel()>2 )
    {
      if(options->getNumberOfOutcomes()==2)cout << "\nRegression " << RegNumber << "\t";
      OutputParams(&cout);
      cout << endl;
    }
  //Output to paramfile after BurnIn
  if( iteration > options->getBurnIn() ){
    OutputParams(&outputstream);
    if(options->getNumberOfOutcomes()< 2 || RegNumber==1) outputstream << endl;
    //output new line in paramfile when last regression model
  }
}
void Regression::OutputErgodicAvg(int samples, std::ofstream *avgstream)const{
 //output to ergodicaveragefile
  if( RegType != None ){
    for( int j = 0; j < NumCovariates; j++ ){
      avgstream->width(9);
      *avgstream << setprecision(6) << SumBeta[j] / samples << "\t";
    }
    avgstream->width(9);
    if( RegType == Linear )
      *avgstream << setprecision(6) << SumLambda / samples << "\t";
  }
}

void Regression::SumParameters(){
  // accumulate sum of parameters after burnin.
  if( NumCovariates > 0 )
      for(int j = 0; j < NumCovariates; ++j)
        SumBeta[j] += beta[j];
  SumLambda += lambda;
}
const double* Regression::getbeta()const{
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

void ExpectedOutcome(const double* const beta, const double* const X, double* Y, int n, int dim, int index, double betaj){
  //given an array of regression parameters beta and covariates X, computes expected outcome Y = X * beta
  double* beta1 = new double[dim];

  for( int i = 0; i < dim; i++ )
    {
      if( i != index )
	beta1[ i ] = beta[i];
      else
	beta1[ i ] = betaj;
    }
  //Xbeta = X * beta1;
  matrix_product(X, beta1, Y, n, dim, 1);
  delete[] beta1;
}

//(unnormalised)log posterior for a single regression parameter
double Regression::lr( const double beta, const void* const vargs )
{
  const BetaArgs* args = (const BetaArgs*)vargs;

  int n = args->n;
  int index = args->index ;
  double beta0 = 0;
  if( index == 0 )
    beta0 = args->beta0;
  double f = 0.0; 

  if(args->coolness > 0.0){

    //log likelihood contribution
    f += args->XtY * beta;

    double *Xbeta = new double[ n ];    
    ExpectedOutcome(args->beta, args->Covariates, Xbeta, n, args->d, index, beta);
    
    for( int i = 0; i < n; i++ ){
      f -= log( 1.0 + exp( Xbeta[ i ] ) );}
    
    //anneal likelihood
    f *= args->coolness;
    
    delete[] Xbeta;
  }

  f -= 0.5 * args->lambda * (beta - beta0) * (beta - beta0); //log prior contribution
  return( f );
}

//first and second derivatives of log posterior for a single regression parameter; 
//used in Newton-Raphson algorithm in MH sampler
double Regression::dlr( const double beta, const void* const vargs )
{
  const BetaArgs* args = (const BetaArgs*)vargs;

  int n = args->n;
  int d = args->d;
  int index = args->index ;
  double beta0 = 0;
  if( index == 0 )
    beta0 = args->beta0;
  double f = 0.0;

  if(args->coolness > 0.0){
    f += args->XtY;
    
    double *Xbeta = new double[ n ];
    
    ExpectedOutcome(args->beta, args->Covariates, Xbeta, n, d, index, beta);
    for( int i = 0; i < n; i++ )
      {
	f -= args->Covariates[ i*d + index ] / ( 1.0 + exp( -Xbeta[ i ] ) );
      }
    delete[] Xbeta;
    f *= args->coolness;
  }

  f -= args->lambda * (beta - beta0);//log prior contribution
  return( f );
}

double Regression::ddlr( const double beta, const void* const vargs )
{
  const BetaArgs* args = (const BetaArgs*)vargs;

  int n = args->n;
  int d = args->d;
  int index = args->index ;

  double f = 0.0;

  if(args->coolness > 0.0){
    double *Xbeta = new double[ n ];
    
    ExpectedOutcome(args->beta, args->Covariates, Xbeta, n, d, index, beta);
    
    for( int i = 0; i < n; i++ )
      {
	f -= args->Covariates[ i*d + index ] * args->Covariates[ i*d + index ] / ( 2.0 + exp( -Xbeta[ i ] ) + exp( Xbeta[ i ] ) );
      }
    delete[] Xbeta;
    f *= args->coolness;
  }

  f -= args->lambda;//log prior contribution
  return( f );
}

double Regression::getDispersion()const{
  //returns dispersion parameter
  double dispersion = 1.0;
  if( RegType == Linear ) dispersion = lambda;//linear regression
  else if(RegType == Logistic) dispersion = 1.0;
  return dispersion;
}

double Regression::getLogLikelihood(const IndividualCollection* const IC)const{
  int NumIndividuals = IC->getSize();
  double loglikelihood = 0.0;
  if(RegType==Linear){    //univariate Gaussian likelihood
    double* dev = new double[NumIndividuals];
    double devsq[1];
    for(int i = 0; i < NumIndividuals; ++i)
      dev[i] = IC->getOutcome(RegNumber, i) - IC->getExpectedY(i, RegNumber);
    matrix_product(dev, dev, devsq, 1, NumIndividuals, 1);
    delete[] dev;
    loglikelihood = -0.5* ( NumIndividuals * (log(2.0*3.14159) - log(lambda)) + lambda * devsq[0] );
  }
  else if(RegType==Logistic){//loglikelihood is sum of logs of bernoulli probabilities, given by EY
    for(int i = 0; i < NumIndividuals; ++i){
      if(IC->getOutcome(RegNumber, i))loglikelihood += log( IC->getExpectedY(i, RegNumber) );
      else loglikelihood += log(1.0 - IC->getExpectedY(i, RegNumber));
    }
  }
  return loglikelihood;

}
double Regression::getLogLikelihoodAtPosteriorMeans(IndividualCollection *IC, int iterations){
  double logL = 0.0;

  //set expected outcome at posterior means of regression parameters
  for(int i = 0; i < NumCovariates; ++i)SumBeta[i] /= (double)iterations; 
  IC->SetExpectedY(RegNumber,SumBeta);//computes X * BetaBar
  for(int i = 0; i < NumCovariates; ++i)SumBeta[i] *= (double)iterations; //restore sumbeta

  //set precision to posterior mean
  double temp = lambda;
  lambda = SumLambda/(double)iterations;

  logL = getLogLikelihood(IC);
  lambda = temp;//restore precision
  IC->SetExpectedY(RegNumber,beta);//restore Xbeta

  return logL;
}

//solves Ax = b by QR decomposition
//on exit R is the inverse of the R matrix in decomposition and V is RR', ready for use later
void Regression::QRSolve(int dim1, int dim2, const double* a, const double* b, double* x){

  //note that views do not allocate memory
  //copy A into QR as QR decomp function destroys argument
  gsl_matrix_view Aview = gsl_matrix_view_array(const_cast<double *>(a), dim1, dim2);
  gsl_matrix* QR = gsl_matrix_alloc(dim1, dim2);
  for(int i1 = 0; i1 < dim1; ++i1)
    for(int i2 = 0; i2 < dim2; ++i2)
      gsl_matrix_set(QR, i1, i2, a[i1*dim2 + i2]); 
  
  gsl_vector* tau = gsl_vector_alloc(dim2);
  gsl_vector_view bview = gsl_vector_view_array(const_cast<double *>(b), dim1);
  gsl_vector_view xview = gsl_vector_view_array(x, dim2);
  
  gsl_vector* residuals = gsl_vector_alloc(dim1);   
  int status = 0;
  std::string errstring;
  try{    
    if(gsl_linalg_QR_decomp(QR, tau))throw ("QR decomp failed...");

    status = gsl_linalg_QR_lssolve(QR, tau, &bview.vector, &xview.vector, residuals);    
    if(status){
      errstring = "QRsolve failed, "; errstring.append( gsl_strerror (status));
      throw(errstring);
    }
    gsl_vector_free(tau);
    gsl_vector_free(residuals);
    
    gsl_permutation* p = gsl_permutation_alloc(dim2);
    gsl_permutation_init(p);//sets to identity permutation
    //copy R into V
    for(int i = 0; i < dim2; ++i){
      V[i*dim2+i] = gsl_matrix_get(QR, i,i);
      for(int j = i+1; j < dim2; ++j){
	V[i*dim2+j] = gsl_matrix_get(QR, i,j);
	V[j*dim2+i] = 0.0;
      }
    }
    gsl_matrix_free(QR);
    
    //V is R, which is its own LU decomposition with identity permutation
    gsl_matrix_view Rview = gsl_matrix_view_array(R, dim2, dim2);
    gsl_matrix_view Vview = gsl_matrix_view_array(V, dim2, dim2);
    if(gsl_linalg_LU_invert(&Vview.matrix, p, &Rview.matrix))throw("Inversion of R failed") ;
    gsl_permutation_free(p);
    //R now contains R^-1
    if(gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1, &Rview.matrix, &Rview.matrix, 0, &Vview.matrix))
      throw("X'X multiplication failed...");
    //V now contains R * R' = inv(X'X) 
  }
  catch(const std::string s){
    //tidy up and throw error message up
    gsl_vector_free(tau);
    gsl_vector_free(residuals);
    gsl_matrix_free(QR);
    std::string error_string = "Regression::QRSolve failed\n";
    error_string.append(s);
    throw (error_string);
    exit(1);
  }
}

void Regression::SamplePrecision(double* lambda, const double* Y, const double* X, int NumIndivs, int NumCovars, double coolness){

  double *Xbeta = new double[NumIndivs];
  //compute X * beta
  matrix_product(X, beta, Xbeta, NumIndivs, NumCovars, 1);

  //compute s^2
  double s2 = 0.0;
  for(int i = 0; i < NumIndivs; ++i){
    double dev = (Y[i] - Xbeta[i]);
    s2 += dev*dev;
  }

  //draw lambda given beta from conjugate Gamma update
  // should replace this with an update that marginalizes over beta
  *lambda = gengam( lambda0 + coolness*0.5*NumIndivs, lambda1 + coolness*0.5*s2);
  
  //cout << "sampled " << *lambda << " from Gamma( " << lambda0 + coolness * NumIndivs << ", " << lambda1 + coolness*0.5*s2 << ")" << endl;
  delete[] Xbeta;
}

void Regression::SampleLinearRegressionParams(double* beta, /*double* lambda, */const double* Y, const double* X, 
					      int NumIndivs, int NumCovars){
  /*
  Samples regression parameters in a linear regression model as described in "Bayesian Data Analysis"
  by Gelman et al (Ch8)
  supply:
  beta - regression parameters
  lambda - precision parameter
  Y - vector of NumIndividuals outcomes
  X - matrix of NumIndividuals * NumCovariates covariates
  */

  // compute betahat using QR decomposition and then V
  try{
    QRSolve(NumIndivs, NumCovars, X, Y, betahat);
  }catch(string s){
    string error_string = "Error occurred while updating Linear Regression parameters:\n";
    error_string.append(s);
    throw(error_string);
  }

  //V currently holds (X'X)^-1
  //scale_matrix(V, 1.0/ *lambda, NumCovars, NumCovars);

  //draw beta from N(betahat, V^-1)
  DrawBeta.SetMean( betahat );
  DrawBeta.SetCovariance( V );
  DrawBeta.Draw(beta);

}

void Regression::SampleLinearRegressionParametersWithAnnealing(const double* Y, const double* X, double* beta, double *lambda, 
							       double coolness){
  //sample precision
  SamplePrecision(lambda, Y, X, NumIndividuals, NumCovariates, coolness);


 //for the case of Var(Y_i) = 1/ (lambda* w_i) 

  //augment X and Y with prior as extra 'data points'
  for(int i = 0; i < NumIndividuals; ++i){
    QY[i] = Y[i] * sqrt(*lambda);
    for(int j = 0; j < NumCovariates; ++j)
      QX[i*NumCovariates +j] = X[i*NumCovariates+j]*sqrt(*lambda);
  }
//   for(int i = 0; i < NumCovariates; ++i){
//     QY[i+NumIndividuals] = betamean[i] * sqrt(betaprecision[i]);//append prior means to Y
//     QX[(i+NumIndividuals)*NumCovariates + i] = sqrt(betaprecision[i]);//prior precision on beta
//     for(int j = i+1; j < NumCovariates; ++j)
//       QX[(i+NumIndividuals)*NumCovariates + j] = QX[(j+NumIndividuals)*NumCovariates + i] = 0.0;
//   }

  //compute Q^{-1/2}Y and Q^{-1/2}X
  for(int i = 0; i < NumIndividuals; ++i){
    QY[i] *= coolness;
    for(int j = 0; j < NumCovariates; ++j)
      QX[i*NumCovariates +j] *= coolness;
  }

  //use standard update algorithm with QX and QY
  SampleLinearRegressionParams(beta, QY, QX, NumIndividuals+NumCovariates, NumCovariates);
}

