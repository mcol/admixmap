/** 
 *   ADMIXMAP
 *   Latent.cc 
 *   Class to hold and update population admixture and sumintensities parameters and their priors
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
  mu = 0;
  sumlogtheta = 0;

  options = op;
  Loci = loci;
  Log = l;
}

void Latent::Initialise(IndividualCollection *individuals, std::ofstream *LogFileStreamPtr,
			std::vector<bool> *_admixed, bool *_symmetric, Vector_d *poptheta, std::string *PopulationLabels){
  //Initialise population admixture distribution Dirichlet parameters alpha
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
     *LogFileStreamPtr << "Shape parameter for sumintensities prior: " <<  options->getRho() << endl;
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
  if(!options->getIndAdmixHierIndicator())  SumAlpha = alpha[0];

  eta = alpha[0].Sum();
  mu = new double[ options->getPopulations() ];
  sumlogtheta = new double[ options->getPopulations() ];
  PopAdmixSampler.SetSize( options->getPopulations() );
  for( int i = 0; i < options->getPopulations(); i++ ){
     mu[i] = alpha[0](i)/eta;
  }

//Initialise sum-of-intensities parameter rho and the parameters of its prior, rhoalpha and rhobeta
  rho = options->getRho();

  if( options->getRho() == 99 ){
     rhoalpha = 1.0;
     rhobeta = 0.0;
     Log->logmsg(true,"Flat prior on sumintensities.\n");
  }
  else if( options->getRho() == 98 ){
     rhoalpha = 0.0;
     rhobeta = 0.0;
     Log->logmsg(true,"Flat prior on log sumintensities.\n");
  }
  else{
    rhoalpha = options->getRho();
    rhobeta = 1;
  }
  rhobeta0 = 1;
  rhobeta1 = 1;

  //Open paramfile 
  if ( options->getIndAdmixHierIndicator()){
    if( strlen( options->getParameterFilename() ) ){
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
	if( options->getTextIndicator() )	InitializeOutputFile(PopulationLabels);
      }
    }
    else{
      Log->logmsg(true,"No paramfile given\n");
      //exit(1);
    }
  }

  //Misc.  
  poptheta->SetNumberOfElements( options->getPopulations() );

  // AlphaParameters is an array with 5 elements
  // element 0 is num individuals or num gametes (if random mating model)
  // element 1 is the sum of the Dirichlet parameter vector
  // where are elements 2 and 3 used?
  // element 4 is the sum of log admixture proportions
   
  if( options->isRandomMatingModel() ){
    AlphaParameters[0] = 2 * individuals->getSize();
  } else {
    AlphaParameters[0] = individuals->getSize();
  }
  //if( options->getAnalysisTypeIndicator() > -1 ){
  AlphaParameters[1] = alpha[0].Sum();
  AlphaParameters[2] = 1;
  AlphaParameters[3] = 1;
  AlphaParameters[4] = 1;
  //}

  Matrix_i empty_i(1,1);
  Matrix_d empty_d(1,1);

   
  DirParamArray = new DARS*[ options->getPopulations() ];
  for( int j = 0; j < options->getPopulations(); j++ ){
    DirParamArray[j] = new DARS();
    DirParamArray[j]->SetParameters( 0, 1, 0.1, AlphaParameters,5,
				     logf, dlogf, ddlogf, empty_i, empty_d );
    DirParamArray[j]->SetLeftTruncation( 0.1 );
  }


  // rho stuff
  rhodata_i.SetNumberOfElements(Loci->GetNumberOfCompositeLoci(), 1 );
  rhodata_d.SetNumberOfElements(Loci->GetNumberOfCompositeLoci(), 1 );
  for(unsigned j = 0; j < Loci->GetNumberOfCompositeLoci(); ++j)rhodata_d(j,0) = Loci->GetDistance(j);
     
   RhoParameters[0] = rhoalpha;
   RhoParameters[1] = rhobeta;
   RhoParameters[2] = Loci->GetNumberOfCompositeLoci();
   // RhoParameters[3] is a sum over individuals 

   RhoDraw = new DARS(0,1,(double)1,RhoParameters,4,frho,dfrho,ddfrho,
                      rhodata_i,rhodata_d);
}

Latent::~Latent()
{
  delete RhoDraw;
  delete[] mu;
  delete[] sumlogtheta;
  for(int i=0; i<options->getPopulations(); i++){
    delete DirParamArray[i];
  }
}

double Latent::sampleForRho()
{
  // Sample for global sum of intensities parameter rho
  // this algorithm is unnecessarily complicated
  // to sample conditional on locus ancestry, we can simply sample the total number of 
  // arrivals on each gamete, then update of rho is conjugate with Poisson likelihood and gamma prior.  
  // But sampling conditional on locus ancestry gives poor mixing
  // it would be better to update rho by a Metropolis random walk conditional on the genotype data 
  // and individual admixture proportions, using the HMM likelihood.  
  RhoDraw->UpdateParameters( RhoParameters,4 );
  RhoDraw->UpdateIntegerData( rhodata_i );
  RhoDraw->UpdateDoubleData( rhodata_d );
  return RhoDraw->Sample();
}     

// these 3 functions calculate log-likelihood and derivatives for adaptive rejection sampling of 
// Dirichlet population admixture parameters
double
Latent::logf( Vector_d &parameters , Matrix_i&, Matrix_d&, double x )
{
  double f = parameters(0) * ( gsl_sf_lngamma( x + parameters(1) ) - gsl_sf_lngamma( x ) ) - x * ( parameters(3) - parameters(4) );//  + (parameters(2) - 1) * log(x);
  
  return(f);
}

double
Latent::dlogf( Vector_d &parameters, Matrix_i&, Matrix_d&, double x )
{
  double f,x2,y1,y2;
  
  x2 = x + parameters(1);
   if(x < 0)cout<<"\nError in Latent::dlogf - arg x to ddigam is negative\n";   
  ddigam( &x , &y1 );
  if(x2 < 0)cout<<"\nError in Latent::dlogf - arg x2 to ddigam is negative\n";   
  ddigam( &x2 , &y2 );
  
  f = parameters(0) * ( y2 - y1 ) - ( parameters(3) - parameters(4) );// + (parameters(2) - 1)/x;
  
  return(f);
}

double
Latent::ddlogf( Vector_d &parameters, Matrix_i&, Matrix_d&, double x )
{
  double f,x2,y1,y2;
  
  x2 = x + parameters(1);
  
  trigam( &x, &y1 );
  trigam( &x2, &y2 );
  
  f = parameters(0) * ( y2 - y1 );// - (parameters(2) - 1)/ (x*x);
  
  return(f);
}

double
Latent::frho( Vector_d &parameters, Matrix_i& xi, Matrix_d &distance, double x )
{
  int genes = (int)parameters(2);
  double f = -x * ( parameters(1) + parameters(3) ) + ( parameters(0) - 1 ) * log(x);

  for( int j = 1; j < genes; j++ ){
    f += xi( j, 0 ) * log(1 - exp( -x * distance( j, 0 ) ) );
  }
  return(f);
}

// these 3 functions calculate log-likelihood and derivatives for adaptive rejection sampling of 
// global rho parameter - will not be needed if this algorithm is replaced
double
Latent::dfrho( Vector_d &parameters, Matrix_i& xi, Matrix_d &distance, double x )
{
  int genes = (int)parameters(2);
  double f = -( parameters(1) + parameters(3) ) + ( parameters(0) - 1 ) / x;

  for( int j = 1; j < genes; j++ ){
    f += xi( j, 0 ) * distance( j, 0 ) * exp( -x * distance( j, 0 ) ) / (1 - exp( -x * distance( j, 0 ) ) );
  }
  return(f);
}

double
Latent::ddfrho( Vector_d &parameters, Matrix_i& xi, Matrix_d &distance, double x )
{
  float temporary;
  int genes = (int)parameters(2);
  double f = -( parameters(0) - 1 ) / (x * x);
  
  for( int j = 1; j < genes; j++ ){
    temporary = exp( x * distance( j, 0 ) ) - 2 + exp( -x * distance( j, 0 ) );
    f -= xi( j, 0 ) * distance( j, 0 ) * distance( j, 0 ) / temporary;
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

}

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

//this should be in InputData
bool Latent::CheckInitAlpha( Vector_d alphatemp )
  // check that Dirichlet parameter vector, if specified by user, has correct length  
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
		    Vector_d *poptheta,std::ofstream *LogFileStreamPtr){

  if( options->getPopulations() > 1 && individuals->getSize() > 1 &&
      options->getIndAdmixHierIndicator() ){
     if( Loci->GetLengthOfGenome() > 0.0 ){
      // Sample for global rho
        if( !options->getRhoIndicator() ){
           //RhoParameters[3] = individuals->GetSumrho0();//equivalent to next line
           RhoParameters[3] = Individual::getSumrho0();
           for(unsigned int j = 0; j < Loci->GetNumberOfCompositeLoci(); ++j)
              //rhodata_i(j,0) = individuals->GetSumXi()[j];//equivalent to next line
              rhodata_i(j,0) = Individual::getSumXi(j);
           rho = sampleForRho();
        }
        else{
           // sample for location parameter of gamma distribution of sumintensities parameters 
           // in population 
           if( options->isRandomMatingModel() )
              rhobeta = gengam( individuals->GetSumrho() + rhobeta1,
                                2*rhoalpha * individuals->getSize() + rhobeta0 );
           else
              rhobeta = gengam( individuals->GetSumrho() + rhobeta1,
                                rhoalpha* individuals->getSize() + rhobeta0 );
        }
     }
     
     // For a model in which the distribution of individual admixture in the population is a mixture
     // of components, we will have one Dirichlet parameter vector for each component, 
     // updated only from those individuals who belong to the component
     // Sample for population admixture distribution Dirichlet parameters alpha
     for( int j = 0; j < options->getPopulations(); j++ ){
//         AlphaParameters[1] -= alpha[0]( j );        
//         AlphaParameters[4] = individuals->getSumLogTheta(j);
//         // elements of Dirichlet parameter vector are updated one at a time
//         DirParamArray[j]->UpdateParameters( AlphaParameters, 5 );
//         alpha[0]( j ) = DirParamArray[j]->Sample();
//         AlphaParameters[1] += alpha[0]( j );
        sumlogtheta[j] = individuals->getSumLogTheta(j);
     }
     PopAdmixSampler.Sample( (unsigned)AlphaParameters[0], sumlogtheta, &eta, mu );
     for( int j = 0; j < options->getPopulations(); j++ ){
        alpha[0]( j ) = mu[j] * eta;
     }
     // accumulate sum of Dirichlet parameter vector over iterations 
     SumAlpha += alpha[0];
  }
  
  if( iteration == options->getBurnIn() && options->getAnalysisTypeIndicator() > 1 ){
     *LogFileStreamPtr << "Individual admixture centred in regression model around: "
                       << *poptheta << endl;
     
     SumAlpha.SetElements(0);
  }
  
  if( iteration < options->getBurnIn() && options->getPopulations() > 1
      && options->getAnalysisTypeIndicator() > 1 ){
     // accumulate ergodic average of population admixture, which is used to centre 
     // the values of individual admixture in the regression model
     *poptheta = SumAlpha / SumAlpha.Sum();
     
  }
  if( iteration > options->getBurnIn() ){
     // accumulate sum of rho parameters after burnin.
     if( options->getPopulations() > 1 ){
        SumRho += rho;
     }
  }
  
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
