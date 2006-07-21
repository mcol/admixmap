#include "LinearRegression.h"
#include <numeric>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

LinearRegression::LinearRegression(){
  lambda0 = 0.0;
  lambda1 = 0.0;
  R = 0;
  QY = 0;
  QX = 0;
  V = 0;
  RegType = Linear;
}

LinearRegression::~LinearRegression(){
  delete[] R;
  delete[] QY;
  delete[] QX;
  delete[] V;
  delete[] betahat;
}

void LinearRegression::Initialise(unsigned Number, double priorPrecision, const IndividualCollection* const individuals, LogWriter &Log){
  Regression::Initialise(Number, priorPrecision, individuals, Log);
  
  // ** Initialise Linear Regression objects    
  if(RegType == Linear){
    lambda0 = 0.01;//shape parameter for prior on lambda
    lambda1 = 0.01;//rate parameter for prior on lambda
    lambda = lambda1 / lambda0;//initialise to prior mean
    
    //fill(betaprecision, betaprecision + NumCovariates, 0.0001);
    
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
  
}

void LinearRegression::InitializeOutputFile(const std::vector<std::string>& CovariateLabels, unsigned NumOutcomes)
{
  Regression::InitializeOutputFile(CovariateLabels, 0);//pass with 0 argument to prevent newline
  //label for precision
  outputstream<< setprecision(6) << "precision\t";
  if(NumOutcomes == RegNumber+1)outputstream << endl;
}

void LinearRegression::Update(bool sumbeta, IndividualCollection* individuals, double coolness
#ifdef PARALLEL
			, MPI::Intracomm &Comm){
  if(Comm.Get_rank() == 0){
#else 
  ){
#endif

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
  
#ifdef PARALLEL
  }
  //broadcast parameters to workers
  Comm.Barrier();
  Comm.Bcast(beta, NumCovariates, MPI::DOUBLE, 0);
  if(RegType == Linear)Comm.Bcast(&lambda, 1, MPI::DOUBLE, 0);
#endif

  individuals->SetExpectedY(RegNumber,beta);

  if(sumbeta){
    SumParameters();
  }
}//end Update

void LinearRegression::OutputParams(ostream* out){
  //if( RegType != None ){
    for( int j = 0; j < NumCovariates; j++ ){
      out->width(9);
      (*out) << setprecision(6) << beta[j] << "\t";
    }
    out->width(9);
    (*out) << setprecision(6) << lambda << "\t";
    //}
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

double LinearRegression::getLogLikelihood(const IndividualCollection* const IC)const{
  int NumIndividuals = IC->getSize();
  double loglikelihood = 0.0;
  //univariate Gaussian likelihood
  double* dev = new double[NumIndividuals];
  double devsq[1];
  for(int i = 0; i < NumIndividuals; ++i)
    dev[i] = IC->getOutcome(RegNumber, i) - IC->getExpectedY(i, RegNumber);
  matrix_product(dev, dev, devsq, 1, NumIndividuals, 1);
  delete[] dev;
  loglikelihood = -0.5* ( NumIndividuals * (log(2.0*3.14159) - log(lambda)) + lambda * devsq[0] );
  return loglikelihood;
}

double LinearRegression::getLogLikelihoodAtPosteriorMeans(IndividualCollection *IC, int iterations){
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
