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
#include "Chromosome.h"
#include <algorithm>
#include <numeric>

using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

#if POPADMIXSAMPLER == 3
static double xlog(double x){
  return log(x);
}
static double xexp(double x){
  return exp(x);
}
#endif


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
  poptheta = 0;
#if POPADMIXSAMPLER == 3
  logalpha = 0;
  initialAlphaStepsize = 0.03;//need a way of setting this without recompiling, or a sensible fixed value
  targetAlphaAcceptRate = 0.5;//need to choose suitable value for this
#endif
}

void Latent::Initialise(int Numindividuals, std::string *PopulationLabels){
  //Initialise population admixture distribution Dirichlet parameters alpha
  alpha = options->getAndCheckInitAlpha(Log);
  SumAlpha.resize( options->getPopulations() );

  if(!options->getIndAdmixHierIndicator())  copy(alpha[0].begin(), alpha[0].end(), SumAlpha.begin());

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
  Log->logmsg(false,"\nShape parameter for gamma prior on sum-of-intensities: ");
  Log->logmsg(false, rhoalpha); Log->logmsg(false,"\n");
  Log->logmsg(false,"Rate (1 / location) parameter for gamma prior on sum-of-intensities: ");
  Log->logmsg(false, rhobeta); Log->logmsg(false,"\n");

  //Open paramfile 
  if ( options->getIndAdmixHierIndicator()){
    if( strlen( options->getParameterFilename() ) ){
      outputstream.open( options->getParameterFilename(), ios::out );
      if( !outputstream )
	{
	  Log->logmsg(true,"ERROR: Couldn't open paramfile\n");
	  exit( 1 );
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

  //ergodic average of population admixture, which is used to centre 
  // the values of individual admixture in the regression model  
  poptheta =new double[ options->getPopulations() ];
  for( int i = 0; i < options->getPopulations(); i++ )poptheta[i] = 0.0;

#if POPADMIXSAMPLER == 1  
  // AlphaParameters is an array with 5 elements
  // element 0 is num individuals or num gametes (if random mating model)
  // element 1 is the sum of the Dirichlet parameter vector
  // elements 2 and 3 are the parameters of the gamma prior
  // element 4 is the sum of log admixture proportions
   
  if( options->isRandomMatingModel() ){
    AlphaParameters[0] = 2 * Numindividuals;
  } else {
    AlphaParameters[0] = Numindividuals;
  }
  //if( options->getAnalysisTypeIndicator() > -1 ){
  AlphaParameters[1] = accumulate(alpha[0].begin(), alpha[0].end(), 0.0, plus<double>());//sum of alpha[0]
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
#elif POPADMIXSAMPLER == 2
  eta = accumulate(alpha[0].begin(), alpha[0].end(), 0.0, std::plus<double>());//eta = sum of alpha[0]
  mu = new double[ options->getPopulations() ];
  PopAdmixSampler.SetSize( options->getPopulations() );
  for( int i = 0; i < options->getPopulations(); i++ ){
    mu[i] = alpha[0][i]/eta;
  }
  if( options->isRandomMatingModel() ){
    obs = 2 * Numindividuals;
  } else {
    obs = Numindividuals;
  }
  
#elif POPADMIXSAMPLER == 3
  logalpha = new double[options->getPopulations()];
  transform(alpha[0].begin(), alpha[0].end(), logalpha, xlog);//logalpha = log(alpha)

  AlphaArgs = new double*[4];
  AlphaArgs[0] = new double[options->getPopulations()];
  for(unsigned i = 1; i < 4;++i)AlphaArgs[i] = new double[1];
  //elem 0 is sum of log admixture props
  AlphaArgs[1][0] = (double)Numindividuals;// elem 1 is num individuals/gametes
  if( options->isRandomMatingModel() )AlphaArgs[1][0] *= 2.0;
  AlphaArgs[2][0] = 1.0;//params of gamma prior
  AlphaArgs[3][0] = 1.0;

  AlphaSampler.SetDimensions(options->getPopulations(), initialAlphaStepsize, 20, targetAlphaAcceptRate, findE, gradE);
#endif

  //set up DARS sampler for global sumintensities
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
  delete[] poptheta;
#if POPADMIXSAMPLER == 1
  for(int i=0; i<options->getPopulations(); i++){
    delete DirParamArray[i];
  }
#elif POPADMIXSAMPLER == 2
  delete[] mu;
#elif POPADMIXSAMPLER == 3
  for(unsigned i = 0; i < 4;++i)delete[] AlphaArgs[i];
  delete[] AlphaArgs;
  delete[] logalpha;
#endif
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

void Latent::UpdateRhoWithRW(IndividualCollection *IC, Chromosome **C){
  //code for generating proposal
  double rhoprop;
  double LogLikelihood = 0.0;
  // ...

  //get log likelihood at current parameter values
  for(int i = 0; i < IC->getSize(); ++i)LogLikelihood -= IC->getIndividual(i)->getLogLikelihood(options, C);
  //get log likelihood at proposal rho and current admixture proportions
  for( unsigned int j = 0; j < Loci->GetNumberOfChromosomes(); j++ ) C[j]->SetLociCorr(rhoprop);
  for(int i = 0; i < IC->getSize(); ++i)LogLikelihood += IC->getIndividual(i)->getLogLikelihood(options, C);
  
  //compute rest of acceptance prob
  // ...
  
  //accept/reject proposal
  // ...
  
  //update TuneRW object
  // ...

  // update f in Chromosomes, must do this regardless of whether proposal is accepted
  for( unsigned int j = 0; j < Loci->GetNumberOfChromosomes(); j++ )
    C[j]->SetLociCorr(rho);
}

#if POPADMIXSAMPLER == 1
// these 3 functions calculate log-likelihood and derivatives for adaptive rejection sampling of 
// Dirichlet population admixture parameters
double
Latent::logf( Vector_d &parameters , Matrix_i&, Matrix_d&, double x )
{
  double f = parameters(0) * ( gsl_sf_lngamma( x + parameters(1) ) - gsl_sf_lngamma( x ) ) - x * ( parameters(3) - parameters(4) ) + (parameters(2) - 1) * log(x);
  
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
  
  f = parameters(0) * ( y2 - y1 ) - ( parameters(3) - parameters(4) ) + (parameters(2) - 1)/x;
  
  return(f);
}

double
Latent::ddlogf( Vector_d &parameters, Matrix_i&, Matrix_d&, double x )
{
  double f,x2,y1,y2;
  
  x2 = x + parameters(1);
  
  trigam( &x, &y1 );
  trigam( &x2, &y2 );
  
  f = parameters(0) * ( y2 - y1 ) - (parameters(2) - 1)/ (x*x);
  
  return(f);
}
#endif

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
    *avgstream << setprecision(6) << SumAlpha[j] / samples << " ";
  }
  avgstream->width(9);
  *avgstream << setprecision(6) << SumRho / samples << " ";
}

void Latent::OutputParams(int iteration){
  //output to logfile for first iteration or every iteration if cout = 0
  if( !options->useCOUT() || iteration == 0 )
    {
      for( int j = 0; j < options->getPopulations(); j++ ){
	Log->width(9);
	Log->write(alpha[0][j], 6);
      }
       if( options->getRhoIndicator() )
	Log->write(rhobeta,6);
      else
	Log->write(rho,6);
    }
  //output to screen
  if( options->useCOUT() )
    {
      for( int j = 0; j < options->getPopulations(); j++ ){
       (cout).width(9);
       (cout) << setprecision(6) << alpha[0][ j ] << " ";
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
      outputstream << setprecision(6) << alpha[0][ j ] << " ";}
    //output rho
    outputstream.width(9);
    if( options->getRhoIndicator() )
      (outputstream) << setprecision(6) << rhobeta << " ";
    else
      (outputstream) << setprecision(6) << rho << " ";

    outputstream << endl;
  }

}

void Latent::Update(int iteration, IndividualCollection *individuals){

  if( options->getPopulations() > 1 && individuals->getSize() > 1 &&
      options->getIndAdmixHierIndicator() ){
    if( Loci->GetLengthOfGenome() > 0.0 ){
      //  ** Sample for global rho **
      if( !options->getRhoIndicator() ){
	//RhoParameters[3] = individuals->GetSumrho0();//equivalent to next line
	RhoParameters[3] = Individual::getSumrho0();
	for(unsigned int j = 0; j < Loci->GetNumberOfCompositeLoci(); ++j)
	  //rhodata_i(j,0) = individuals->GetSumXi()[j];//equivalent to next line
	  rhodata_i(j,0) = Individual::getSumXi(j);
	//rhodata_i.SetColumn( 0, individuals->GetSumXi() );
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
    // ** Sample for population admixture distribution Dirichlet parameters, alpha **
           
    // For a model in which the distribution of individual admixture in the population is a mixture
    // of components, we will have one Dirichlet parameter vector for each component, 
    // updated only from those individuals who belong to the component

#if POPADMIXSAMPLER == 1
    for( int j = 0; j < options->getPopulations(); j++ ){
      AlphaParameters[1] -= alpha[0][ j ];
      AlphaParameters[4] = individuals->getSumLogTheta(j);
      // elements of Dirichlet parameter vector are updated one at a time
      DirParamArray[j]->UpdateParameters( AlphaParameters, 5 );
      alpha[0][ j ] = DirParamArray[j]->Sample();
      AlphaParameters[1] += alpha[0][ j ];
    }
#elif POPADMIXSAMPLER == 2
    PopAdmixSampler.Sample( obs, individuals->getSumLogTheta(), &eta, mu );
    for( int j = 0; j < options->getPopulations(); j++ )
      alpha[0][j] = mu[j]*eta;
    
#elif POPADMIXSAMPLER == 3
    //for( int j = 0; j < options->getPopulations(); j++ ){
    //AlphaArgs[0][j] = individuals->getSumLogTheta(j);
    //}
    copy(individuals->getSumLogTheta(), individuals->getSumLogTheta()+options->getPopulations(), AlphaArgs[0]);
    AlphaSampler.Sample(logalpha, AlphaArgs);//sample new values for logalpha
    transform(logalpha, logalpha+options->getPopulations(), alpha[0].begin(), xexp);//alpha = exp(logalpha)
    if(!((iteration+1) % 10)){
      AlphaSampler.Tune();//tune Hamiltonian sampler every 10 iterations
    }
#endif
    
    // ** accumulate sum of Dirichlet parameter vector over iterations  **
    transform(alpha[0].begin(), alpha[0].end(), SumAlpha.begin(), SumAlpha.begin(), std::plus<double>());//SumAlpha += alpha[0];
  }

  if( iteration == options->getBurnIn() && options->getAnalysisTypeIndicator() > 1 ){
    Log->write("Individual admixture centred in regression model around: ");
    Log->write( poptheta, options->getPopulations());
    Log->write("\n");

    fill(SumAlpha.begin(), SumAlpha.end(), 0.0);
  }

  if( iteration < options->getBurnIn() && options->getPopulations() > 1
      && options->getAnalysisTypeIndicator() > 1 ){
    // accumulate ergodic average of population admixture, which is used to centre 
    // the values of individual admixture in the regression model
    double sum = accumulate(SumAlpha.begin(), SumAlpha.end(), 0.0);
    for( int j = 0; j < options->getPopulations(); j++ )poptheta[j] = SumAlpha[j] / sum;

  }
  if( iteration > options->getBurnIn() ){
    // accumulate sum of rho parameters after burnin.
    if( options->getPopulations() > 1 ){
      SumRho += rho;
    }
  }

}
//end Update

vector<double > &Latent::getalpha0(){
  return alpha[0];
}
std::vector<vector<double> > &Latent::getalpha(){
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
const double *Latent::getpoptheta(){
  return poptheta;
}

#if POPADMIXSAMPLER == 3
float Latent::getAlphaSamplerAcceptanceRate(){
  return AlphaSampler.getAcceptanceRate();
}
float Latent::getAlphaSamplerStepsize(){
  return AlphaSampler.getStepsize();
}
//calculate objective function (-log posterior) for log alpha, used in Hamiltonian Metropolis algorithm
double Latent::findE(unsigned dim, const double* const theta, const double* const* args){
  /*
    theta = log dirichlet parameters (alpha)
    args[1] = n = #individuals/gametes
    args[0] = sumlogtheta (array, length dim) = sums of logs of individual admixture proportions
    args[2] = eps0, args[3] = eps1 = parameters of Gamma prior for alpha
  */

  double E = 0.0;
  double sumalpha = 0.0, sumgamma = 0.0, sumtheta = 0.0, sume = 0.0;
  bool flag = true;
  for(unsigned j = 0; j < dim;++j){
    if(exp(theta[j]) == 0.0){flag = false;break;} //to avoid underflow problems
    sumalpha += exp(theta[j]);
    sumgamma += gsl_sf_lngamma(exp(theta[j]));
    sume += exp(theta[j]) * (args[3][0] - args[0][j]);
    sumtheta += theta[j];
  }
  if(flag){
    E = args[1][0] * (gsl_sf_lngamma(sumalpha) - sumgamma) - sume + args[2][0] * sumtheta;
    return -E;
  }
  else return -1.0;//is there a better return value? possibly use flag pointer
}

//calculate gradient for log alpha
void Latent::gradE(unsigned dim, const double* const theta, const double* const* args, double *g){
  double sumalpha = 0.0, x, y1, y2;
  for(unsigned j = 0; j < dim; ++j) {
    g[j] = 0.0;
    sumalpha += exp(theta[j]);
  }
    ddigam(&sumalpha, &y1);
    for(unsigned j = 0; j < dim; ++j) {
      x = exp(theta[j]);
      ddigam(&x, &y2);
      if(x > 0.0 && gsl_finite(y1) && gsl_finite(y2)){//to avoid over/underflow problems
	g[j] = x *( args[1][0] *(y2 - y1) + (args[3][0] - args[0][j])) - args[2][0];
      }
    }
}

#endif
