#include "Latent.h"

using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

Latent::Latent( AdmixOptions * op, Genome *loci, LogWriter *l)
{
  options = 0;
  rho = 0.0;
  rhoalpha = 0.0;
  rhobeta = 0.0;
  SumRho = 0.0;

  options = op;
  Loci = loci;
  Log = l;
}

void Latent::Initialise(IndividualCollection *individuals, std::ofstream *LogFileStreamPtr,
			std::vector<bool> *_admixed, bool *_symmetric, Vector_d *poptheta){

  Vector_d alphatemp;
  SumAlpha.SetNumberOfElements( options->getPopulations() );

  _admixed->resize(2,true);
  *_symmetric = true;
  if( options->sizeInitAlpha() == 0 ){
     alphatemp.SetNumberOfElements( options->getPopulations() );
     alphatemp.SetElements( 1.0 );
     alpha.resize(2,alphatemp);
     *LogFileStreamPtr << "Prior for gamete/individual admixture: "
                      << alphatemp << endl;
     *LogFileStreamPtr << "Shape parameter for rho prior: " <<  options->getRho() << endl;
  }
  else if( options->sizeInitAlpha() == 1 ){
    alphatemp = options->getInitAlpha(0);
    (*_admixed)[0] = CheckInitAlpha( alphatemp );
     alpha.resize(2,alphatemp);
     *LogFileStreamPtr << "Prior for gamete/individual admixture: "
                      << alphatemp << endl;
  }
  else if( options->getAnalysisTypeIndicator() < 0 ){
    alphatemp = options->getInitAlpha(0);
    (*_admixed)[0] = CheckInitAlpha( alphatemp );
     alpha.push_back(alphatemp);

     *LogFileStreamPtr << "Prior for gamete 1 admixture: "
                      << alphatemp << endl;
     alphatemp = options->getInitAlpha(1);
     (*_admixed)[1] = CheckInitAlpha( alphatemp );
     alpha.push_back(alphatemp);

     *LogFileStreamPtr << "Prior for gamete 2 admixture: "
                      << alphatemp << endl;
     *_symmetric = false;
  }
  else{
     Log->logmsg(true,"Can only specify seperate priors for gamete admixture with analysis of single individual.\n");
  }

  if( !options->getRhoIndicator() )
    Log->logmsg(true,"Model with global rho.\n");
  else if( options->getModelIndicator() )
    Log->logmsg(true,"Model with gamete specific rho.\n");
  else
    Log->logmsg(true,"Model with individual specific rho.\n");

  rho = options->getRho();
//   f->SetNumberOfElements( Loci->GetNumberOfCompositeLoci() );
//   for( int j = 1; j < Loci->GetNumberOfCompositeLoci(); j++ )
//     (*f)(j) = strangExp( -Loci->GetDistance( j ) * rho );
  if( options->getRho() == 99 ){
     rhoalpha = 1.0;
     rhobeta = 0.0;
     Log->logmsg(true,"Flat prior on rho.\n");
  }
  else if( options->getRho() == 98 ){
     rhoalpha = 0.0;
     rhobeta = 0.0;
     Log->logmsg(true,"Flat prior on log rho.\n");
  }
  else{
    rhoalpha = options->getRho();
    rhobeta = 1;
  }
  rhobeta0 = 1;
  rhobeta1 = 1;

  //Open paramfile 
  if ( strlen( options->getParameterFilename() ) ){
    outputstream.open( options->getParameterFilename(), ios::out );
    if( !outputstream )
      {
	Log->logmsg(true,"ERROR: Couldn't open paramfile\n");
	//exit( 1 );
      }
    else{
      Log->logmsg(true,"Writing population-level parameters to ");
      Log->logmsg(true,options->getParameterFilename());
      Log->logmsg(true,"\n");
    }
  }
  else{
    Log->logmsg(true,"No paramfile given\n");
    //exit(1);
  }

  //Open IndAdmixtureOutputFile    
  if ( strlen( options->getIndAdmixtureFilename() ) ){
    Log->logmsg(true,"Writing individual-level parameters to ");
    Log->logmsg(true,options->getIndAdmixtureFilename());
    Log->logmsg(true,"\n");
    //indadmixoutput = new IndAdmixOutputter(options,Loci,PopulationLabels);
  } else {
    Log->logmsg(true,"No indadmixturefile given\n");
  }

  //Misc.  
  poptheta->SetNumberOfElements( options->getPopulations() );

  PreUpdate(individuals);
}

Latent::~Latent()
{
  delete RhoDraw;
  for(int i=0; i<options->getPopulations(); i++){
    delete DirParamArray[i];
  }
  /**
   * From Admixbase destructor
   **/
 
  //    if( options->getAnalysisTypeIndicator() > 1 ){
  //      delete [] TargetLabels;
  //    }
}

//FUNCTION TO BE CALLED BEFORE UPDATE TO INITIALISE MEMBERS USED IN LOOP
void Latent::PreUpdate(IndividualCollection *individuals){

  MatrixArray_i empty_i(1);
  MatrixArray_d empty_d(1);
   
  AlphaParameters.SetNumberOfElements(5);
  if( options->getModelIndicator() ){
    AlphaParameters(0) = 2 * individuals->getSize();
  } else {
    AlphaParameters(0) = individuals->getSize();
  }
  //if( options->getAnalysisTypeIndicator() > -1 ){
  AlphaParameters(1) = alpha[0].Sum();
     AlphaParameters(2) = 1;
     AlphaParameters(3) = 1;
     AlphaParameters(4) = 1;
     //}
   
  DirParamArray = new DARS*[ options->getPopulations() ];
  for( int j = 0; j < options->getPopulations(); j++ ){
    DirParamArray[j] = new DARS();
    DirParamArray[j]->SetParameters( 0, 1, 0.1, AlphaParameters,
				     logf, dlogf, ddlogf, empty_i, empty_d );
    DirParamArray[j]->SetLeftTruncation( 0.1 );
  }


  // rho stuff
  rhodata_i.SetNumberOfElementsWithDimensions( 1, Loci->GetNumberOfCompositeLoci(), 1 );
  rhodata_d.SetNumberOfElementsWithDimensions( 1, Loci->GetNumberOfCompositeLoci(), 1 );
  rhodata_d(0).SetColumn( 0, Loci->GetDistances().Double() );
   
   RhoParameters.SetNumberOfElements(4);
   RhoParameters(0) = rhoalpha;
   RhoParameters(1) = rhobeta;
   RhoParameters(2) = Loci->GetNumberOfCompositeLoci();

   RhoDraw = new DARS(0,1,(double)1,RhoParameters,frho,dfrho,ddfrho,
                            rhodata_i,rhodata_d);

 
  
}//END OF PREUPDATE

void
Latent::load_f(double rho,Vector_d *f,Genome *chrm){
  int locus = 0;
  for( int j = 0; j < chrm->size(); j++ ){
    locus++;
    for( int jj = 1; jj < (*chrm)(j)->GetSize(); jj++ ){
      (*f)(locus) = exp( -Loci->GetDistance( locus ) * rho );
      locus++;
    }
  }
}

double Latent::sampleForRho(Vector_d& RhoParameters, DARS* RhoDraw,
                            MatrixArray_i& rhodata_i, MatrixArray_d& rhodata_d)
{
  // Sample for global sum of intensities parameter rho
  RhoDraw->UpdateParameters( RhoParameters );
  RhoDraw->UpdateIntegerData( rhodata_i );
  RhoDraw->UpdateDoubleData( rhodata_d );
  return RhoDraw->Sample();
}     

double
Latent::logf( Vector_d &parameters , MatrixArray_i&, MatrixArray_d&, double x )
{
  double f = parameters(0) * ( gsl_sf_lngamma( x + parameters(1) ) - gsl_sf_lngamma( x ) ) - x * ( parameters(3) - parameters(4) )  + (parameters(2) - 1) * log(x);
  
  return(f);
}

double
Latent::dlogf( Vector_d &parameters, MatrixArray_i&, MatrixArray_d&, double x )
{
  double f,x2,y1,y2;
  
  x2 = x + parameters(1);
   if(x < 0)cout<<"\nError in Latent::dlogf - arg x to ddigam is negative\n";   
  ddigam( &x , &y1 );
  if(x2 < 0)cout<<"\nError in Latent::dlogf - arg x2 to ddigam is negative\n";   
  ddigam( &x2 , &y2 );
  
  f = parameters(0) * ( y2 - y1 ) - ( parameters(3) - parameters(4) ) + (parameters(2) - 1)/x;
  
  return(f);
}

double
Latent::ddlogf( Vector_d &parameters, MatrixArray_i&, MatrixArray_d&, double x )
{
  double f,x2,y1,y2;
  
  x2 = x + parameters(1);
  
  trigam( &x, &y1 );
  trigam( &x2, &y2 );
  
  f = parameters(0) * ( y2 - y1 ) - (parameters(2) - 1)/ (x*x);
  
  return(f);
}

double
Latent::frho( Vector_d &parameters, MatrixArray_i& xi, MatrixArray_d &distance, double x )
{
  int genes = (int)parameters(2);
  double f = -x * ( parameters(1) + parameters(3) ) + ( parameters(0) - 1 ) * log(x);

  for( int j = 1; j < genes; j++ ){
    f += xi(0)( j, 0 ) * log(1 - exp( -x * distance(0)( j, 0 ) ) );
  }
  return(f);
}

double
Latent::dfrho( Vector_d &parameters, MatrixArray_i& xi, MatrixArray_d &distance, double x )
{
  int genes = (int)parameters(2);
  double f = -( parameters(1) + parameters(3) ) + ( parameters(0) - 1 ) / x;

  for( int j = 1; j < genes; j++ ){
    f += xi(0)( j, 0 ) * distance(0)( j, 0 ) * exp( -x * distance(0)( j, 0 ) ) / (1 - exp( -x * distance(0)( j, 0 ) ) );
  }
  return(f);
}

double
Latent::ddfrho( Vector_d &parameters, MatrixArray_i& xi, MatrixArray_d &distance, double x )
{
  float temporary;
  int genes = (int)parameters(2);
  double f = -( parameters(0) - 1 ) / (x * x);
  
  for( int j = 1; j < genes; j++ ){
    temporary = exp( x * distance(0)( j, 0 ) ) - 2 + exp( -x * distance(0)( j, 0 ) );
    f -= xi(0)( j, 0 ) * distance(0)( j, 0 ) * distance(0)( j, 0 ) / temporary;
  }
  return(f);
}

void Latent::InitializeOutputFile(std::string *PopulationLabels)
{
  // Header line of paramfile

//   if(options->getAnalysisTypeIndicator() < 0){//Analysis for a single individual or set of individuals
//     outputstream << "\"Log Likelihood\"\t \"Log Posterior\"\n";
//   }
  if( options->getAnalysisTypeIndicator() >= 0 ){

    //Pop. Admixture
    for( int i = 0; i < options->getPopulations(); i++ ){
      outputstream << PopulationLabels[i] << " ";
    }
    //SumIntensities
    if( !options->getRhoIndicator() )
      outputstream << "\"sumIntensities\" ";
    else
      outputstream << "\"sumIntensities.beta\" ";

    outputstream << endl;
  }

}//end InitializeOutputFile

void Latent::OutputErgodicAvg( int samples, std::ofstream *avgstream)
{
  for( int j = 0; j < options->getPopulations(); j++ ){
    avgstream->width(9);
    *avgstream << setprecision(6) << SumAlpha(j) / samples << " ";
  }
  avgstream->width(9);
  *avgstream << setprecision(6) << SumRho / samples << " ";
}

void Latent::OutputParams(int iteration, std::ofstream *LogFileStreamPtr){
  //output to logfile
  if( !options->useCOUT() || iteration == 0 )
    {
      for( int j = 0; j < options->getPopulations(); j++ ){
	LogFileStreamPtr->width(9);
	(*LogFileStreamPtr) << setprecision(6) << alpha[0]( j ) << " ";
      }
      LogFileStreamPtr->width(9);
      if( options->getRhoIndicator() )
	(*LogFileStreamPtr) << setprecision(6) << rhobeta << " ";
      else
	(*LogFileStreamPtr) << setprecision(6) << rho << " ";
    }
  //output to screen
  if( options->useCOUT() )
    {
      for( int j = 0; j < options->getPopulations(); j++ ){
       (cout).width(9);
       (cout) << setprecision(6) << alpha[0]( j ) << " ";
      }
     (cout).width(9);
      if( options->getRhoIndicator() )
	(cout) << setprecision(6) << rhobeta << " ";
      else
	(cout) << setprecision(6) << rho << " ";
    }
  //Output to paramfile after BurnIn
    //output alpha
  if( iteration > options->getBurnIn() ){
    for( int j = 0; j < options->getPopulations(); j++ ){
      outputstream.width(9);
      outputstream << setprecision(6) << alpha[0]( j ) << " ";}
    //output rho
    outputstream.width(9);
    if( options->getRhoIndicator() )
      (outputstream) << setprecision(6) << rhobeta << " ";
    else
      (outputstream) << setprecision(6) << rho << " ";

    outputstream << endl;
  }

}

double Latent::strangExp( double x )
{
  double y;
  if( x > -700 )
    y = exp(x);
  else
    y = 0;
  return( y );
}

bool
Latent::CheckInitAlpha( Vector_d alphatemp )
{
   bool admixed = true;
   int count=0;
   for( int i = 0; i < alphatemp.GetNumberOfElements(); i++ )
      if( alphatemp(i) != 0.0 )
         count++;
   if( count == 1 )
      admixed = false;
   if( alphatemp.GetNumberOfElements() != options->getPopulations() ){
      cout << "Error in specification of alpha.\n"
           << alphatemp << endl;
      exit(0);
   }
   return admixed;
}

void Latent::Update(int iteration, IndividualCollection *individuals,
		    Vector_d *SumLogTheta,Vector_d *poptheta,std::ofstream *LogFileStreamPtr){

  if( options->getPopulations() > 1 && individuals->getSize() > 1 &&
      options->getAnalysisTypeIndicator() != -3 ){
    if( Loci->GetLengthOfGenome() > 0.0 ){
      // Sample for global rho
      if( !options->getRhoIndicator() ){
	RhoParameters(3) = individuals->GetSumrho0();
	rhodata_i(0).SetColumn( 0, individuals->GetSumXi() );
	rho = sampleForRho(RhoParameters,RhoDraw,rhodata_i,rhodata_d);
      }
      else{
	if( options->getModelIndicator() )
	  rhobeta = gengam( individuals->GetSumrho() + rhobeta1,
			    2*rhoalpha * individuals->getSize() + rhobeta0 );
	else
	  rhobeta = gengam( individuals->GetSumrho() + rhobeta1,
			    rhoalpha* individuals->getSize() + rhobeta0 );
      }
    }
         
// Sample for population admixture distribution Dirichlet parameters alpha
    for( int j = 0; j < options->getPopulations(); j++ ){
      AlphaParameters(1) -= alpha[0]( j );
      AlphaParameters(4) = (*SumLogTheta)( j );
      DirParamArray[j]->UpdateParameters( AlphaParameters );
      alpha[0]( j ) = DirParamArray[j]->Sample();
      AlphaParameters(1) += alpha[0]( j );
    }
    SumAlpha += alpha[0];
  }

  if( iteration == options->getBurnIn() && options->getAnalysisTypeIndicator() > 1 ){
    *LogFileStreamPtr << "Individual admixture centred in regression model around: "
		     << *poptheta << endl;
    Loci->ResetSumAlleleFreqs();
    SumAlpha.SetElements(0);
  }

  if( iteration < options->getBurnIn() && options->getPopulations() > 1
      && options->getAnalysisTypeIndicator() > 0 ){
    *poptheta = SumAlpha / SumAlpha.Sum();
  }
  if( iteration > options->getBurnIn() ){
    // accumulate sum of rho parameters after burnin.
    if( options->getPopulations() > 1 ){
      SumRho += rho;
    }
  }

//       if( !options->getRhoIndicator() ){
// 	load_f(rho,f,chrm);
//       }
}
//end Update

Vector_d *Latent::getalpha0(){
  return &alpha[0];
}
std::vector<Vector_d> Latent::getalpha(){
  return alpha;
}

double Latent::getrhoalpha(){
  return rhoalpha;
}
double Latent::getrhobeta(){
  return rhobeta;
}
double Latent::getrho(){
  return rho;
}
