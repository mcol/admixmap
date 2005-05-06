#include "Regression.h"
#include "matrix_d.h"

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

void Regression::Initialise(IndividualCollection *individuals,AdmixOptions *options, std::string *PopulationLabels, LogWriter *Log){
  AnalysisTypeIndicator = options->getAnalysisTypeIndicator();
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
      if( options->getTextIndicator() )InitializeOutputFile(options, individuals, PopulationLabels);
    }
  }
  else{
    Log->logmsg(true,"No regparamfile given\n");
    //exit(1);
  }



  if( AnalysisTypeIndicator < 2 )
    NoCovariates = 0;
  else{
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
    if( AnalysisTypeIndicator > 1 ){
      SumBeta.SetNumberOfElementsWithDimensions( individuals->getTargetSize(), NoCovariates, 1 );
      beta.SetNumberOfElementsWithDimensions( individuals->getTargetSize(), NoCovariates, 1 );
     }

    double p;
    beta0.SetNumberOfElementsWithDimensions( individuals->getTargetSize(), NoCovariates, 1 );
    for( int k = 0; k < individuals->getTargetSize(); k++ ){
      p = (individuals->getTargetCol(k,0)).Mean();
      if( individuals->getOutcomeType(k) )
	beta0(k)(0,0) = log( p / ( 1 - p ) );
      else
	beta0(k)(0,0) = p;
    }
    beta = beta0;
    
    n0.SetNumberOfElements( NoCovariates, NoCovariates );
    n0.SetDiagonal( 1 );
  }

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
  if( AnalysisTypeIndicator == 3 || AnalysisTypeIndicator == 4 || AnalysisTypeIndicator == 5){
    Matrix_i empty_i;
    for( int k = 0; k < individuals->getTargetSize(); k++ ){
      if( AnalysisTypeIndicator != 5 || (AnalysisTypeIndicator == 5 && individuals->getOutcomeType(k) ) ){
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

void Regression::Update(IndividualCollection *individuals){
  // Sample for regression model parameters beta
  for( int k = 0; k < individuals->getTargetSize(); k++ ){
    //      linear
    if( AnalysisTypeIndicator == 2 || (AnalysisTypeIndicator == 5 && ! individuals->getOutcomeType(k)) ){
      Matrix_d temporary = individuals->getCovariates().Transpose() * individuals->getCovariates() + n0;
      temporary.InvertUsingLUDecomposition();
      temporary.Symmetrize();
      betan = temporary * ( n0 * beta0(k) + individuals->getCovariates().Transpose() * individuals->getOutcome(k) );
      double lambdan = lambda0 + 0.5 *
	( (individuals->getOutcome(k) - individuals->getCovariates() * betan ).Transpose() * individuals->getOutcome(k)
	  + (beta0(k) - betan).Transpose() * n0 * beta0(k) )(0,0);
      lambda(k) = gengam( lambdan, lambda0 + 0.5 * individuals->getSize() );
      DrawBeta.SetMean( betan );
      DrawBeta.SetCovariance( temporary / lambda(k) );
      beta(k) = DrawBeta.Draw();
    }
    //      logistic
    else if( AnalysisTypeIndicator == 3 ||
	     AnalysisTypeIndicator == 4 ||
	     (AnalysisTypeIndicator == 5 &&  individuals->getOutcomeType(k) ) ){
      sum = (individuals->getOutcome(k)).Transpose() * individuals->getCovariates();
      for( int j = 0; j < NoCovariates; j++ ){
	BetaDrawArray[j]->UpdateDoubleData( individuals->getCovariates() );
	BetaParameters( NoCovariates + 2 ) = j;
	BetaParameters( NoCovariates ) = sum( 0, j );
	BetaDrawArray[j]->UpdateParameters( BetaParameters );
	acceptbeta(k) = BetaDrawArray[j]->Sample( &( beta(k)(j,0) ) );
	BetaParameters(j) = beta(k)( j, 0 );
      }
    }
    //ExpectedY(k) = individuals->getCovariates * beta(k);
    individuals->SetExpectedY(k,beta(k));
    if( individuals->getOutcomeType(k) )
      // calculateExpectedY( ExpectedY(k), individuals->getSize() );
      individuals->calculateExpectedY(k);
  }

}//end Update

void Regression::InitializeOutputFile(AdmixOptions *options, IndividualCollection *individuals,std::string *PopulationLabels)
{
  // Header line of paramfile
  if( AnalysisTypeIndicator == 2 || AnalysisTypeIndicator == 3 ){
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
    if( AnalysisTypeIndicator == 2 )//AnalysisType = 2 => linear regression
      outputstream << setprecision(6) << "\"precision\"";
  }
  else if( AnalysisTypeIndicator == 5 ){
    for( int kk = 0; kk < 2; kk++ ){
      outputstream << "\"intercept\" ";
      for( int i = 0; i < individuals->GetNumberOfInputCols(); i++ ){
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

void Regression::Output(int iteration, std::ofstream *LogFileStreamPtr, AdmixOptions *options, IndividualCollection *individuals){
  //output to logfile
  if( !options->useCOUT() || iteration == 0 )
    {
      if( AnalysisTypeIndicator == 2 || AnalysisTypeIndicator == 3 )
	{
          for( int j = 0; j < NoCovariates; j++ )
	    {
	      LogFileStreamPtr->width(9);
	      (*LogFileStreamPtr) << setprecision(6) << beta(0)( j, 0 ) << " ";
	    }
          LogFileStreamPtr->width(9);
          if( AnalysisTypeIndicator == 2 )
	    {
	      (*LogFileStreamPtr) << setprecision(6) << lambda;
	    }
	}
      else if( AnalysisTypeIndicator == 5 )
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
	      if( ! individuals->getOutcomeType(k) )
		{
                  *LogFileStreamPtr<< setprecision(6) << lambda(k);
		}
	    }
	}
    }
  //output to screen
  if( options->useCOUT() )
    {
     if( AnalysisTypeIndicator == 2 || AnalysisTypeIndicator == 3 ){
	for( int j = 0; j < NoCovariates; j++ ){
	  (cout).width(9);
	  cout << setprecision(6) << beta(0)( j, 0 ) << " ";
	}
	(cout).width(9);
	if( AnalysisTypeIndicator == 2 )
	  cout << setprecision(6) << lambda;
      }
      else if( AnalysisTypeIndicator == 5 ){
	for( int k = 0; k < individuals->getTargetSize(); k++ ){
	  cout << "\nRegression " << k << " ";
	  for( int j = 0; j < NoCovariates; j++ ){
	    (cout).width(9);
	    cout << setprecision(6) << beta(k)( j, 0 ) << " ";
	  }
	  (cout).width(9);
	  if( ! individuals->getOutcomeType(k) )
	    cout << setprecision(6) << lambda(k);
	}
      }
    }
  //Output to paramfile after BurnIn
 if( iteration > options->getBurnIn() ){
   if( AnalysisTypeIndicator == 2 || AnalysisTypeIndicator == 3 ){
     for( int j = 0; j < NoCovariates; j++ ){
       outputstream.width(9);
       outputstream << setprecision(6) << beta(0)( j, 0 ) << " ";
     }
     outputstream.width(9);
     if( AnalysisTypeIndicator == 2 )
       outputstream << setprecision(6) << lambda(0) << " ";
   }
   else if( AnalysisTypeIndicator == 5 ){
     for( int k = 0; k < individuals->getTargetSize(); k++ ){
       for( int j = 0; j < NoCovariates; j++ ){
	 outputstream.width(9);
	 outputstream << setprecision(6) << beta(k)( j, 0 ) << " ";
       }
       if( ! individuals->getOutcomeType(k) ){
	 outputstream.width(9);
	 outputstream<< setprecision(6) << lambda(k) << " ";
       }
     }
   }
      outputstream << endl;
 }
}
void Regression::OutputErgodicAvg(int samples, IndividualCollection *individuals, std::ofstream *avgstream){
 //output to ergodicaveragefile
  if( AnalysisTypeIndicator == 2 || AnalysisTypeIndicator == 3 ){
    for( int j = 0; j < NoCovariates; j++ ){
      avgstream->width(9);
      *avgstream << setprecision(6) << SumBeta(0)( j, 0 ) / samples << " ";
    }
    avgstream->width(9);
    if( AnalysisTypeIndicator == 2 )
      *avgstream << setprecision(6) << SumLambda(0) / samples << " ";
  }
  else if( AnalysisTypeIndicator == 5 ){
    for( int k = 0; k < individuals->getTargetSize(); k++ ){
      for( int j = 0; j < NoCovariates; j++ ){
	avgstream->width(9);
	*avgstream << setprecision(6) << SumBeta(k)( j, 0 ) / samples << " ";
      }
      if( ! individuals->getOutcomeType(k) ){
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

int Regression::getNoCovariates(){
  return NoCovariates;
}

double
Regression::lr( Vector_d &parameters, Matrix_i &, Matrix_d &data, double beta )
{
  int n = (int)data.GetNumberOfRows();
  int d = (int)data.GetNumberOfCols();
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
  
  Xbeta = data * beta1;
  for( int i = 0; i < n; i++ ){
    f -= log( 1. + exp( Xbeta( i, 0 ) ) );}
  
  return( f );
}

double
Regression::dlr( Vector_d &parameters, Matrix_i&, Matrix_d &data, double beta )
{
  int n = (int)data.GetNumberOfRows();
  int d = (int)data.GetNumberOfCols();
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
  
  Xbeta = data * beta1;
  for( int i = 0; i < n; i++ )
    {
      f -= data( i, index ) / ( 1. + exp( -Xbeta( i, 0 ) ) );
    }
  return( f );
}

double
Regression::ddlr( Vector_d &parameters, Matrix_i&, Matrix_d &data, double beta )
{
  int n = (int)data.GetNumberOfRows();
  int d = (int)data.GetNumberOfCols();
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
  
  Xbeta = data * beta1;
  for( int i = 0; i < n; i++ )
    {
      f -= data( i, index ) * data( i, index ) / ( 2. + exp( -Xbeta( i, 0 ) ) + exp( Xbeta( i, 0 ) ) );
    }
  return( f );
}

double Regression::getDispersion(int OutcomeType){
  //return 1.0/dispersion parameter
  //Note: OutcomeType = individuals->getOutcomeType(0)
  //                  = 0 =>linear regression
  //                  = 1=> logistic
  double dispersion = 1.0;
  if( AnalysisTypeIndicator == 2 ) dispersion = lambda(0);//linear regression
  else {
    if( AnalysisTypeIndicator == 3 || AnalysisTypeIndicator == 4 )dispersion = 1.0;//logistic
    else if (AnalysisTypeIndicator == 5)dispersion = OutcomeType ? 1.0 : lambda(0);
    //else cout<<"Invalid call to Regression::getDispersion()."<<endl;
  }
  return dispersion;
}

