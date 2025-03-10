/** 
 *   LinearRegression.cc 
 *   Class to represent and update parameters of a logistic regression model
 *   Copyright (c) 2006-2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "bclib/LinearRegression.h"
#include "bclib/LogWriter.h"
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include "bclib/linalg.h"
#include "bclib/rand.h" // for sampling data precision
#include <cmath>
#include <iomanip>

using namespace::std;
BEGIN_BCLIB_NAMESPACE

LinearRegression::LinearRegression(unsigned Number, double priorPrecision, const bclib::DataMatrix& Covars, const bclib::DataMatrix& Outcome, 
				   LogWriter &Log): Regression(Number, Linear){
  lambda0 = 0.0;
  lambda1 = 0.0;
  R = 0;
  QY = 0;
  QX = 0;
  V = 0;
  //RegType = Linear;
  Initialise(priorPrecision, Covars, Outcome, Log);
}

LinearRegression::~LinearRegression(){
  delete[] R;
  delete[] QY;
  delete[] QX;
  delete[] V;
  delete[] betahat;
}

void LinearRegression::Initialise(double priorPrecision, const bclib::DataMatrix& Covars, const bclib::DataMatrix& Outcome, 
				  LogWriter &Log){
  Regression::Initialise(Covars.nCols(), Covars.nRows(), Covars.getData());
  Log.setDisplayMode(Quiet);
  //set prior precision
  double outcomeSampleVariance = Outcome.getSampleVariance(RegNumber);
  betaprecision[0] = priorPrecision / outcomeSampleVariance;
  Log << "\nGaussian priors on linear regression parameters with zero means and precisions\n ("<< betaprecision[0];
  
  for(int j = 1; j < NumCovariates; ++j){
    //get sample variance of covariate
    double svc = Covars.getSampleVariance(j);
    if (svc < 0.0001) svc = 1.0;// incase covar is sampled and begins constant or near constant over individuals
    
    betaprecision[j] = priorPrecision * svc / outcomeSampleVariance;
    Log << ", " << betaprecision[j];
  }
  Log << ")\n";
  
  SetExpectedY(beta);

  // ** Initialise Linear Regression objects    
  lambda0 = 0.01;//shape parameter for prior on lambda
  lambda1 = 0.01;//rate parameter for prior on lambda
  lambda = lambda1 / lambda0;//initialise to prior mean
  
  //fill(betaprecision, betaprecision + NumCovariates, 0.0001);
  //std::vector<double> v = individuals->getOutcome(RegNumber);
  //double p = accumulate(v.begin(), v.end(), 0.0, std::plus<double>()) / (double)v.size();
  //betamean[0] = p;
  
  
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

void LinearRegression::InitializeOutputFile(const std::vector<std::string>& CovariateLabels, unsigned NumOutcomes)
{
  Regression::InitializeOutputFile(CovariateLabels, 0);//pass with 0 argument to prevent newline
  //label for precision
  outputstream<< setprecision(6) << "precision";
  if(NumOutcomes == RegNumber+1)outputstream << bclib::newline;
}

void LinearRegression::Update(bool sumbeta, const std::vector<double>& Outcome, double coolness  ){
  // Sample for regression model parameters beta
  //and precision in linear regression
  
  SampleLinearRegressionParametersWithAnnealing(&(Outcome[0]), Covariates, beta, &lambda, coolness);

    if(sumbeta){
      SumParameters();
    }
    SetExpectedY(beta);
}//end Update
  
void LinearRegression::SetExpectedY(const double* const _beta){
  //sets ExpectedY = X * Beta
  if(ExpectedY){
    matrix_product(Covariates, _beta, ExpectedY, NumIndividuals, NumCovariates, 1);
  }
}

void LinearRegression::OutputParams(bclib::Delimitedostream& out)const{
  Regression::OutputParams(out);
  out << lambda;
}

void LinearRegression::OutputErgodicAvg(int samples, std::ofstream& avgstream)const{
  //output to ergodicaveragefile
  Regression::OutputErgodicAvg(samples, avgstream);
  avgstream.width(9);
  avgstream << setprecision(6) << SumLambda / samples << "\t";
}

//solves Ax = b by QR decomposition
//on exit R is the inverse of the R matrix in decomposition and V is RR', ready for use later
void LinearRegression::QRSolve(int dim1, int dim2, const double* a, const double* b, double* x){

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
  }
}

void LinearRegression::SamplePrecision(double* lambda, const double* Y, const double* X, int NumIndivs, int NumCovars, double coolness){

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
  *lambda = Rand::gengam( lambda0 + coolness*0.5*NumIndivs, lambda1 + coolness*0.5*s2);
  
  //cout << "sampled " << *lambda << " from Gamma( " << lambda0 + coolness * NumIndivs << ", " << lambda1 + coolness*0.5*s2 << ")" << endl;
  delete[] Xbeta;
}

void LinearRegression::SampleLinearRegressionParams(double* beta, /*double* lambda, */const double* Y, const double* X, 
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

void LinearRegression::SampleLinearRegressionParametersWithAnnealing(const double* Y, const double* X, double* beta, double *lambda, 
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
double LinearRegression::getDispersion()const{
  return lambda;
}
///returns Derivative of Inverse Link Function for individual i
double LinearRegression::DerivativeInverseLinkFunction(unsigned)const{
  return 1.0;    
}

double LinearRegression::getLogLikelihood(const std::vector<double>& Outcome)const{
  double loglikelihood = 0.0;
  //univariate Gaussian likelihood
  double* dev = new double[NumIndividuals];
  double devsq[1];
  for(int i = 0; i < NumIndividuals; ++i)
    dev[i] = Outcome[ i ] - ExpectedY[ i ];
  matrix_product(dev, dev, devsq, 1, NumIndividuals, 1);
  delete[] dev;
  loglikelihood = -0.5* ( NumIndividuals * (log(2.0*3.14159) - log(lambda)) + lambda * devsq[0] );
  return loglikelihood;
}

double LinearRegression::getLogLikelihoodAtPosteriorMeans(int iterations, const std::vector<double>& Outcome){
  double logL = 0.0;

  //set expected outcome at posterior means of regression parameters
  for(int i = 0; i < NumCovariates; ++i)SumBeta[i] /= (double)iterations; 
  //IC->SetExpectedY(RegNumber,SumBeta);//computes X * BetaBar
  SetExpectedY(SumBeta);
  for(int i = 0; i < NumCovariates; ++i)SumBeta[i] *= (double)iterations; //restore sumbeta

  //set precision to posterior mean
  double temp = lambda;
  lambda = SumLambda/(double)iterations;

  logL = getLogLikelihood(Outcome);
  lambda = temp;//restore precision
  //IC->SetExpectedY(RegNumber,beta);//restore Xbeta
  SetExpectedY(beta);

  return logL;
}
END_BCLIB_NAMESPACE
