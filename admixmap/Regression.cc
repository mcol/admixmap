#include "Regression.h"
#include "matrix_d.h"

Regression::Regression(){
  lambda0 = 0.0; //population
  lambda1 = 0.0; //population
  beta = 0;
  beta0 = 0;
  SumBeta = 0;
  NumOutcomeVars = 0;
  NoCovariates = 0;
  BetaDrawArray = 0;
}

Regression::Regression(int nocovariates){
  NoCovariates = nocovariates;
  lambda0 = 0.0; //population
  lambda1 = 0.0; //population
  beta = 0;
  beta0 = 0;
  SumBeta = 0;
  NumOutcomeVars = 0;
  BetaDrawArray = 0;
}

Regression::~Regression(){
  for(int i = 0; i < NoCovariates; i++){
    delete BetaDrawArray[i];
  }
  delete[] BetaDrawArray;
  //delete[] beta;
  free_matrix(beta,NumOutcomeVars);
  free_matrix(SumBeta,NumOutcomeVars);
  delete[] beta0;
}

void Regression::Initialise(IndividualCollection *individuals,AdmixOptions *options, std::string *PopulationLabels, LogWriter *Log){
  AnalysisTypeIndicator = options->getAnalysisTypeIndicator();
  NumOutcomeVars = individuals->getNumberOfOutcomeVars();
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
      if( options->getTestForAdmixtureAssociation() )
	NoCovariates = individuals->GetNumberOfInputCols() + 1;
      else
	NoCovariates = individuals->GetNumberOfInputCols() + options->getPopulations();}
    else{
      if( options->getTestForAdmixtureAssociation() )
	NoCovariates = 1;
      else
	NoCovariates = options->getPopulations();
    }

    //  Prior parameters for linear regression
    if( AnalysisTypeIndicator > 1 ){
      beta = alloc2D_d(NumOutcomeVars, NoCovariates);
      SumBeta = alloc2D_d(NumOutcomeVars, NoCovariates);
     }
    for(int i = 0; i < NumOutcomeVars; ++i)
      for(int j = 0; j < NoCovariates; ++j)
	SumBeta[i][j] = 0.0;

    double p;
    beta0 = new Matrix_d[NumOutcomeVars];
      for(int i = 0; i < NumOutcomeVars; ++i){
	beta0[i].SetNumberOfElements(NoCovariates, 1 );
      }
    for( int k = 0; k < NumOutcomeVars; k++ ){
      p = (individuals->getTargetCol(k,0)).Mean();
      if( individuals->getOutcomeType(k) )
	beta0[k](0,0) = log( p / ( 1 - p ) );
      else
	beta0[k](0,0) = p;
    }
      for(int i = 0; i < NumOutcomeVars; ++i)for(int j = 0; j < NoCovariates; ++j){
	beta[i][j] = beta0[i](j,0);
      }
    
    n0.SetNumberOfElements( NoCovariates, NoCovariates );
    n0.SetDiagonal( 1 );
  }

  lambda0 = .1;
  lambda1 = .1;
  
  lambda.SetNumberOfElements( NumOutcomeVars );
  lambda.SetElements( .1 );

  if( AnalysisTypeIndicator == 2 || AnalysisTypeIndicator == 5) {
    Log->logmsg(false,"\nNormal-inverse-gamma prior for linear regression model with gamma shape parameter ");
    Log->logmsg(false, lambda0); 
    Log->logmsg(false, " and rate parameter "); Log->logmsg(false, lambda1); Log->logmsg(false, "\n");
  }
  if( AnalysisTypeIndicator == 3 ||  AnalysisTypeIndicator == 4 ||  AnalysisTypeIndicator == 5) {
    Log->logmsg(false,"\nGaussian priors on logistic regression parameters with precision ");
    Log->logmsg(false, "0.1"); Log->logmsg(false, "\n");
  }
  
  SumLambda.SetNumberOfElements( NumOutcomeVars );
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
    for( int k = 0; k < NumOutcomeVars; k++ ){
      if( AnalysisTypeIndicator != 5 || (AnalysisTypeIndicator == 5 && individuals->getOutcomeType(k) ) ){
	BetaParameters.SetElements(0);
	BetaParameters( NoCovariates + 1 ) = lambda(k);
	BetaParameters( NoCovariates + 3 ) = beta0[k](0,0);
	for( int i = 0; i < NoCovariates; i++ ){
	  BetaDrawArray[i] = new MetropolisHastings( BetaParameters, lr, dlr, ddlr, empty_i, individuals->getCovariates() );
	}
      }
    }
  }
}

void Regression::Update(bool afterBurnIn, IndividualCollection *individuals){
  // Sample for regression model parameters beta
  for( int k = 0; k < NumOutcomeVars; k++ ){
    //      linear
    if( AnalysisTypeIndicator == 2 || (AnalysisTypeIndicator == 5 && ! individuals->getOutcomeType(k)) ){
      Matrix_d temporary = individuals->getCovariates().Transpose() * individuals->getCovariates() + n0;
      temporary.InvertUsingLUDecomposition();
      temporary.Symmetrize();
      betan = temporary * ( n0 * beta0[k] + individuals->getCovariates().Transpose() * individuals->getOutcome(k) );
      double lambdan = lambda1 + 0.5 *
	( (individuals->getOutcome(k) - individuals->getCovariates() * betan ).Transpose() * individuals->getOutcome(k)
	  + (beta0[k] - betan).Transpose() * n0 * beta0[k] )(0,0);
      lambda(k) = gengam( lambdan, lambda0 + 0.5 * individuals->getSize() );
      DrawBeta.SetMean( betan );
      DrawBeta.SetCovariance( temporary / lambda(k) );
      DrawBeta.Draw(beta[k]);
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
	acceptbeta(k) = BetaDrawArray[j]->Sample( &( beta[k][j] ) );
	BetaParameters(j) = beta[k][j];
      }
    }
    //ExpectedY(k) = individuals->getCovariates * beta(k);
    individuals->SetExpectedY(k,beta[k]);
    if( individuals->getOutcomeType(k) )
      // calculateExpectedY( ExpectedY(k), individuals->getSize() );
      individuals->calculateExpectedY(k);
  }
  if(afterBurnIn)
    SumParameters();
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
    if( !options->getTestForAdmixtureAssociation() ){
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

void Regression::Output(int iteration, AdmixOptions *options, IndividualCollection *individuals, LogWriter *Log){
  //output to logfile
  if( !options->useCOUT() || iteration == 0 )
    {
      if( AnalysisTypeIndicator == 2 || AnalysisTypeIndicator == 3 )
	{
          for( int j = 0; j < NoCovariates; j++ )
	    {
	      Log->width(9);
	      Log->write(beta[0][j],6);
	    }
          //Log->width(9);
          if( AnalysisTypeIndicator == 2 )
	    {
	      for( int k = 0; k <  NumOutcomeVars; k++ )
		Log->write(lambda(k),6);
		//Log->write( lambda, 6);
	    }
	}
      else if( AnalysisTypeIndicator == 5 )
	{
          for( int k = 0; k <  NumOutcomeVars; k++ )
	    {
	      Log->write("\nRegression ");Log->write(k);
	      for( int j = 0; j < NoCovariates; j++ )
		{
		  Log->width(9);
		  Log->write(beta[k][j], 9);
		}
	      Log->width(9);
	      if( ! individuals->getOutcomeType(k) )
		{
                  Log->write(lambda(k),6);
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
	  cout << setprecision(6) << beta[0][j] << " ";
	}
	(cout).width(9);
	if( AnalysisTypeIndicator == 2 )
	  cout << setprecision(6) << lambda;
      }
      else if( AnalysisTypeIndicator == 5 ){
	for( int k = 0; k < NumOutcomeVars; k++ ){
	  cout << "\nRegression " << k << " ";
	  for( int j = 0; j < NoCovariates; j++ ){
	    (cout).width(9);
	    cout << setprecision(6) << beta[k][j] << " ";
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
       outputstream << setprecision(6) << beta[0][j] << " ";
     }
     outputstream.width(9);
     if( AnalysisTypeIndicator == 2 )
       outputstream << setprecision(6) << lambda(0) << " ";
   }
   else if( AnalysisTypeIndicator == 5 ){
     for( int k = 0; k < NumOutcomeVars; k++ ){
       for( int j = 0; j < NoCovariates; j++ ){
	 outputstream.width(9);
	 outputstream << setprecision(6) << beta[k][j] << " ";
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
      *avgstream << setprecision(6) << SumBeta[0][j] / samples << " ";
    }
    avgstream->width(9);
    if( AnalysisTypeIndicator == 2 )
      *avgstream << setprecision(6) << SumLambda(0) / samples << " ";
  }
  else if( AnalysisTypeIndicator == 5 ){
    for( int k = 0; k < NumOutcomeVars; k++ ){
      for( int j = 0; j < NoCovariates; j++ ){
	avgstream->width(9);
	*avgstream << setprecision(6) << SumBeta[k][j] / samples << " ";
      }
      if( ! individuals->getOutcomeType(k) ){
	avgstream->width(9);
	*avgstream << setprecision(6) << SumLambda(k) / samples << " ";
      }
    }
  }
}

void Regression::SumParameters(){
  // accumulate sum of parameters after burnin.
  if( NoCovariates )
    for(int i = 0; i < NumOutcomeVars; ++i){
      for(int j = 0; j < NoCovariates; ++j)
        SumBeta[i][j] += beta[i][j];
    }
  //if( AnalysisTypeIndicator > 1 )
    SumLambda += lambda;
}
double **Regression::getbeta(){
  if(beta)
    return beta;
  else return NULL;
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

