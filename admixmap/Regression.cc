#include "Regression.h"

Regression::Regression(){
  lambda0 = 0.0; //population
  lambda1 = 0.0; //population
}

Regression::Regression(int nocovariates){
  NoCovariates = nocovariates;
  lambda0 = 0.0; //population
  lambda1 = 0.0; //population
}

Regression::~Regression(){
  for(int i=0;i<NoCovariates;i++){
    delete BetaDrawArray[i];
  }
  delete [] BetaDrawArray;

}

void Regression::Initialise(IndividualCollection *individuals,AdmixOptions *options, LogWriter *Log){
  //Open paramfile 
  if ( strlen( options->getRegressionOutputFilename() ) ){
    outputstream.open( options->getRegressionOutputFilename(), ios::out );
    if( !outputstream )
      {
	Log->logmsg(true,"ERROR: Couldn't open regparamfile\n");
	//exit( 1 );
      }
    else{
      Log->logmsg(true,"Writing population-level regression parameters to ");
      Log->logmsg(true,options->getRegressionOutputFilename());
      Log->logmsg(true,"\n");
    }
  }
  else{
    Log->logmsg(true,"No regparamfile given\n");
    //exit(1);
  }

  MatrixArray_i empty_i(1);

  if( options->getAnalysisTypeIndicator() < 2 )
    NoCovariates = 0;
  else{//AnalysisTypeIndicator >= 2
    if( individuals->GetNumberOfInputRows() == individuals->getSize() ){
      if( options->getScoreTestIndicator() )
	NoCovariates = individuals->GetNumberOfInputCols() + 1;
      else
	NoCovariates = individuals->GetNumberOfInputCols() + options->getPopulations();}
    else{
      if( options->getScoreTestIndicator() )
	NoCovariates = 1;
      else
	NoCovariates = options->getPopulations();
    }

    //  Prior parameters for linear regression
    if( options->getAnalysisTypeIndicator() > 1 ){
      SumBeta.SetNumberOfElementsWithDimensions( individuals->getTargetSize(), NoCovariates, 1 );
      beta.SetNumberOfElementsWithDimensions( individuals->getTargetSize(), NoCovariates, 1 );
     }

    double p;
    beta0.SetNumberOfElementsWithDimensions( individuals->getTargetSize(), NoCovariates, 1 );
    for( int k = 0; k < individuals->getTargetSize(); k++ ){
      p = (individuals->getTargetCol(k,0)).Mean();
      if( individuals->getTargetType(k) )
	beta0(k)(0,0) = log( p / ( 1 - p ) );
      else
	beta0(k)(0,0) = p;
    }
    beta = beta0;
    
    n0.SetNumberOfElements( NoCovariates, NoCovariates );
    n0.SetDiagonal( 1 );
  }//end case AnalysisTypeIndicator >= 2

  lambda0 = .1;
  lambda1 = .1;
  lambda.SetNumberOfElements( individuals->getTargetSize() );
  lambda.SetElements( .1 );

  SumLambda.SetNumberOfElements( individuals->getTargetSize() );
  SumLambda.SetElements( .1 );    
  acceptbeta.SetNumberOfElements( options->getPopulations());
  DrawBeta.SetDimension( NoCovariates );
  BetaParameters.SetNumberOfElements( NoCovariates + 4 );
  BetaDrawArray = new MetropolisHastings*[NoCovariates];

  for( int i = 0; i < NoCovariates; i++ ){
    BetaDrawArray[i] = 0;
  }
  if( options->getAnalysisTypeIndicator() == 3 || options->getAnalysisTypeIndicator() == 4 || options->getAnalysisTypeIndicator() == 5){
    for( int k = 0; k < individuals->getTargetSize(); k++ ){
      if( options->getAnalysisTypeIndicator() != 5 || (options->getAnalysisTypeIndicator() == 5 && individuals->getTargetType(k) ) ){
	BetaParameters.SetElements(0);
	BetaParameters( NoCovariates + 1 ) = lambda(k);
	BetaParameters( NoCovariates + 3 ) = beta0(k)(0,0);
	for( int i = 0; i < NoCovariates; i++ ){
	  BetaDrawArray[i] = new MetropolisHastings( BetaParameters, lr, dlr, ddlr, empty_i, individuals->getCovariates() );
	}
      }
    }
  }
}

void Regression::Update(int AnalysisTypeIndicator,IndividualCollection *individuals){
  // Sample for regression model parameters beta
  for( int k = 0; k < individuals->getTargetSize(); k++ ){
    //      linear
    if( AnalysisTypeIndicator == 2 || (AnalysisTypeIndicator == 5 && ! individuals->getTargetType(k)) ){
      Matrix_d temporary = individuals->getCovariates(0).Transpose() * individuals->getCovariates(0) + n0;
      temporary.InvertUsingLUDecomposition();
      temporary.Symmetrize();
      betan = temporary * ( n0 * beta0(k) + individuals->getCovariates(0).Transpose() * individuals->getTarget(k) );
      double lambdan = lambda0 + 0.5 *
	( (individuals->getTarget(k) - individuals->getCovariates(0) * betan ).Transpose() * individuals->getTarget(k)
	  + (beta0(k) - betan).Transpose() * n0 * beta0(k) )(0,0);
      lambda(k) = gengam( lambdan, lambda0 + 0.5 * individuals->getSize() );
      DrawBeta.SetMean( betan );
      DrawBeta.SetCovariance( temporary / lambda(k) );
      beta(k) = DrawBeta.Draw();
    }
    //      logistic
    else if( AnalysisTypeIndicator == 3 ||
	     AnalysisTypeIndicator == 4 ||
	     (AnalysisTypeIndicator == 5 &&  individuals->getTargetType(k) ) ){
      sum = (individuals->getTarget(k)).Transpose() * individuals->getCovariates(0);
      for( int j = 0; j < NoCovariates; j++ ){
	BetaDrawArray[j]->UpdateDoubleData( individuals->getCovariates() );
	BetaParameters( NoCovariates + 2 ) = j;
	BetaParameters( NoCovariates ) = sum( 0, j );
	BetaDrawArray[j]->UpdateParameters( BetaParameters );
	acceptbeta(k) = BetaDrawArray[j]->Sample( &( beta(k)(j,0) ) );
	BetaParameters(j) = beta(k)( j, 0 );
      }
    }
    //ExpectedY(k) = individuals->getCovariates(0) * beta(k);
    individuals->SetExpectedY(k,beta(k));
    if( individuals->getTargetType(k) )
      // calculateExpectedY( ExpectedY(k), individuals->getSize() );
      individuals->calculateExpectedY(k);
  }

}//end Update

void Regression::InitializeOutputFile(AdmixOptions *options, IndividualCollection *individuals,std::string *PopulationLabels)
{
  // Header line of paramfile
  if( options->getAnalysisTypeIndicator() == 2 || options->getAnalysisTypeIndicator() == 3 ){
    outputstream << "\"intercept\" ";
   
    if( individuals->GetNumberOfInputRows() == individuals->getSize() ){
      for( int i = 0; i < individuals->GetNumberOfInputCols(); i++ ){
	outputstream << individuals->getCovariateLabels(i) << " ";
      }
    } 
    if( !options->getScoreTestIndicator() ){
      for( int k = 1; k < options->getPopulations(); k++ ){
	outputstream << "\"slope." << PopulationLabels[k].substr(1) << " ";
      }
    }
    if( options->getAnalysisTypeIndicator() == 2 )
      outputstream << setprecision(6) << "\"precision\"";
  }
  else if( options->getAnalysisTypeIndicator() == 5 ){
    for( int kk = 0; kk < 2; kk++ ){
      outputstream << "\"intercept\" ";
      for( int i = 0; i < individuals->GetNumberOfInputCols(); i++ ){
	outputstream << individuals->getCovariateLabels(i) << " ";
      }
      for( int k = 1; k < options->getPopulations(); k++ ){
	outputstream << "\"slope." << PopulationLabels[k].substr(1) << " ";
      }
      if( ! individuals->getTargetType(kk) ){
	outputstream<< setprecision(6) << "\"precision\"";
      }
    }
  }
  outputstream << endl;
}

void Regression::Output(int iteration, std::ofstream *LogFileStreamPtr, AdmixOptions *options, IndividualCollection *individuals){
  //output to logfile
  if( !options->useCOUT() || iteration == 0 )
    {
      if( options->getAnalysisTypeIndicator() == 2 || options->getAnalysisTypeIndicator() == 3 )
	{
          for( int j = 0; j < NoCovariates; j++ )
	    {
	      LogFileStreamPtr->width(9);
	      (*LogFileStreamPtr) << setprecision(6) << beta(0)( j, 0 ) << " ";
	    }
          LogFileStreamPtr->width(9);
          if( options->getAnalysisTypeIndicator() == 2 )
	    {
	      (*LogFileStreamPtr) << setprecision(6) << lambda;
	    }
	}
      else if( options->getAnalysisTypeIndicator() == 5 )
	{
          for( int k = 0; k < individuals->getTargetSize(); k++ )
	    {
	      *LogFileStreamPtr<< "\nRegression " << k << " ";
	      for( int j = 0; j < NoCovariates; j++ )
		{
		  LogFileStreamPtr->width(9);
		  *LogFileStreamPtr<< setprecision(6) << beta(k)( j, 0 ) << " ";
		}
	      LogFileStreamPtr->width(9);
	      if( ! individuals->getTargetType(k) )
		{
                  *LogFileStreamPtr<< setprecision(6) << lambda(k);
		}
	    }
	}
    }
  //output to screen
  if( options->useCOUT() )
    {
     if( options->getAnalysisTypeIndicator() == 2 || options->getAnalysisTypeIndicator() == 3 ){
	for( int j = 0; j < NoCovariates; j++ ){
	  (cout).width(9);
	  cout << setprecision(6) << beta(0)( j, 0 ) << " ";
	}
	(cout).width(9);
	if( options->getAnalysisTypeIndicator() == 2 )
	  cout << setprecision(6) << lambda;
      }
      else if( options->getAnalysisTypeIndicator() == 5 ){
	for( int k = 0; k < individuals->getTargetSize(); k++ ){
	  cout << "\nRegression " << k << " ";
	  for( int j = 0; j < NoCovariates; j++ ){
	    (cout).width(9);
	    cout << setprecision(6) << beta(k)( j, 0 ) << " ";
	  }
	  (cout).width(9);
	  if( ! individuals->getTargetType(k) )
	    cout << setprecision(6) << lambda(k);
	}
      }
    }
  //Output to paramfile after BurnIn
 if( iteration > options->getBurnIn() ){
   if( options->getAnalysisTypeIndicator() == 2 || options->getAnalysisTypeIndicator() == 3 ){
     for( int j = 0; j < NoCovariates; j++ ){
       outputstream.width(9);
       outputstream << setprecision(6) << beta(0)( j, 0 ) << " ";
     }
     outputstream.width(9);
     if( options->getAnalysisTypeIndicator() == 2 )
       outputstream << setprecision(6) << lambda(0) << " ";
   }
   else if( options->getAnalysisTypeIndicator() == 5 ){
     for( int k = 0; k < individuals->getTargetSize(); k++ ){
       for( int j = 0; j < NoCovariates; j++ ){
	 outputstream.width(9);
	 outputstream << setprecision(6) << beta(k)( j, 0 ) << " ";
       }
       if( ! individuals->getTargetType(k) ){
	 outputstream.width(9);
	 outputstream<< setprecision(6) << lambda(k) << " ";
       }
     }
   }
      outputstream << endl;
 }
}
void Regression::OutputErgodicAvg(int samples, AdmixOptions *options, IndividualCollection *individuals, std::ofstream *avgstream){
 //output to ergodicaveragefile
  if( options->getAnalysisTypeIndicator() == 2 || options->getAnalysisTypeIndicator() == 3 ){
    for( int j = 0; j < NoCovariates; j++ ){
      avgstream->width(9);
      *avgstream << setprecision(6) << SumBeta(0)( j, 0 ) / samples << " ";
    }
    avgstream->width(9);
    if( options->getAnalysisTypeIndicator() == 2 )
      *avgstream << setprecision(6) << SumLambda(0) / samples << " ";
  }
  else if( options->getAnalysisTypeIndicator() == 5 ){
    for( int k = 0; k < individuals->getTargetSize(); k++ ){
      for( int j = 0; j < NoCovariates; j++ ){
	avgstream->width(9);
	*avgstream << setprecision(6) << SumBeta(k)( j, 0 ) / samples << " ";
      }
      if( ! individuals->getTargetType(k) ){
	avgstream->width(9);
	*avgstream << setprecision(6) << SumLambda(k) / samples << " ";
      }
    }
  }
}

void Regression::SumParameters(int AnalysisTypeIndicator){
  // accumulate sum of parameters after burnin.
  if( NoCovariates )
    SumBeta += beta;
  if( AnalysisTypeIndicator > 1 )
    SumLambda += lambda;
}
MatrixArray_d *Regression::getbeta(){
  return &beta;
}
Vector_d *Regression::getlambda(){
  return &lambda;
}
double Regression::getlambda0(){
  return lambda(0);
}
MatrixArray_d *Regression::getSumBeta(){
  return &SumBeta;
}
Vector_d *Regression::getSumLambda(){
  return &SumLambda;
}
int Regression::getNoCovariates(){
  return NoCovariates;
}

double
Regression::lr( Vector_d &parameters, MatrixArray_i &, MatrixArray_d &data, double beta )
{
  int n = (int)data(0).GetNumberOfRows();
  int d = (int)data(0).GetNumberOfCols();
  int index = (int)parameters( d + 2 );
  double beta0 = 0;
  if( index == 0 )
    beta0 = parameters( d + 3 );
  Matrix_d Xbeta, beta1( d, 1 );
  double f = parameters( d ) * beta
    - 0.5 * parameters( d + 1 ) * (beta - beta0) * (beta - beta0);
  
  for( int i = 0; i < d; i++ )
    {
      if( i != index )
	beta1( i , 0 ) = parameters(i);
      else
	beta1( i , 0 ) = beta;
    }
  
  Xbeta = data(0) * beta1;
  for( int i = 0; i < n; i++ ){
    f -= log( 1. + exp( Xbeta( i, 0 ) ) );}
  
  return( f );
}

double
Regression::dlr( Vector_d &parameters, MatrixArray_i&, MatrixArray_d &data, double beta )
{
  int n = (int)data(0).GetNumberOfRows();
  int d = (int)data(0).GetNumberOfCols();
  int index = (int)parameters( d + 2 );
  double beta0 = 0;
  if( index == 0 )
    beta0 = parameters( d + 3 );
  double f = parameters( d ) - parameters( d + 1 ) * (beta - beta0);
  Matrix_d Xbeta, beta1( d, 1 );
  
  for( int i = 0; i < d; i++ )
    {
      if( i != index )
	beta1( i, 0 ) = parameters(i);
      else
	beta1( i, 0 ) = beta;
    }
  
  Xbeta = data(0) * beta1;
  for( int i = 0; i < n; i++ )
    {
      f -= data(0)( i, index ) / ( 1. + exp( -Xbeta( i, 0 ) ) );
    }
  return( f );
}

double
Regression::ddlr( Vector_d &parameters, MatrixArray_i&, MatrixArray_d &data, double beta )
{
  int n = (int)data(0).GetNumberOfRows();
  int d = (int)data(0).GetNumberOfCols();
  int index = (int)parameters( d + 2 );
  double f = -parameters( d + 1 );
  Matrix_d Xbeta, beta1( d, 1 );
  
  for( int i = 0; i < d; i++ )
    {
      if( i != index )
	beta1( i, 0 ) = parameters(i);
      else
	beta1( i, 0 ) = beta;
    }
  
  Xbeta = data(0) * beta1;
  for( int i = 0; i < n; i++ )
    {
      f -= data(0)( i, index ) * data(0)( i, index ) / ( 2. + exp( -Xbeta( i, 0 ) ) + exp( Xbeta( i, 0 ) ) );
    }
  return( f );
}
