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
#include "functions.h"
#include <algorithm>
#include <numeric>
#include "gsl/gsl_math.h"

using namespace std;

#define PR(x) cerr << #x << " = " << x << endl;

Latent::Latent( AdmixOptions* op, const Genome* const loci)
{
  options = 0;
  rho = 0.0;
  rhoalpha = 0.0;
  rhobeta = 0.0;
  SumLogRho = 0.0;
  options = op;
  Loci = loci;
  poptheta = 0;
#if POPADMIXSAMPLER == 2

  SumLocusAncestry = 0;
#elif POPADMIXSAMPLER == 3
  logalpha = 0;
  initialAlphaStepsize = 0.05;//need a way of setting this without recompiling, or a sensible fixed value
  targetAlphaAcceptRate = 0.44;//need to choose suitable value for this
#endif
}

void Latent::Initialise(int Numindividuals, const std::string* const PopulationLabels, LogWriter &Log){
  Log.setDisplayMode(On);
  // ** Initialise population admixture distribution Dirichlet parameters alpha **
  int K = options->getPopulations();
  //ergodic average of population admixture, which is used to centre 
  // the values of individual admixture in the regression model  
  poptheta = new double[ K ];
  for( int i = 0; i < K; i++ )poptheta[i] = 0.0;
  alpha = options->getInitAlpha();
  SumAlpha.resize( K );
  if(!options->getIndAdmixHierIndicator())  copy(alpha[0].begin(), alpha[0].end(), SumAlpha.begin());

  if(K > 1){
    // ** set up sampler for alpha **
    
#if POPADMIXSAMPLER == 1  
    // AlphaParameters is an array with 5 elements
    // element 0 is num individuals or num gametes (if random mating model)
    // element 1 is the sum of the Dirichlet parameter vector
    // elements 2 and 3 are the parameters of the gamma prior
    // element 4 is the sum of log admixture proportions
    
    double alphapriormean = options->getAlphamean();
    double alphapriorvar = options->getAlphavar();
    if( options->isRandomMatingModel() ){
      AlphaParameters[0] = 2 * Numindividuals;
    } else {
      AlphaParameters[0] = Numindividuals;
    }
    //if( NumIndividuals > 1 ){
    AlphaParameters[1] = accumulate(alpha[0].begin(), alpha[0].end(), 0.0, plus<double>());//sum of alpha[0]
    AlphaParameters[2] = alphapriormean*alphapriormean / alphapriorvar; // shape parameter of gamma prior
    AlphaParameters[3] = alphapriormean / alphapriorvar; // rate parameter of gamma prior
    AlphaParameters[4] = 1;
    //}
    
    
    DirParamArray = new AdaptiveRejection*[ K ];
    for( int j = 0; j < K; j++ ){
      DirParamArray[j] = new Adaptiverejection();
      DirParamArray[j]->Initialise(false, true, 0, 1, logf, dlogf);
      DirParamArray[j]->setLowerBound( 0.1 );
    }
#elif POPADMIXSAMPLER == 2
    PopAdmixSampler.SetSize( Numindividuals, K,  options->getAlphamean(), options->getAlphavar());

    if( options->isRandomMatingModel() ){
      obs = 2 * Numindividuals;
    } else {
      obs = Numindividuals;
    }
    SumLocusAncestry = new int[Numindividuals*K];

#endif
    
    // ** Initialise sum-of-intensities parameter rho and the parameters of its prior, rhoalpha and rhobeta **
    rhobeta0 = options->getRhobetaShape();
    rhobeta1 = options->getRhobetaRate();
    rhobeta = options->getRhobeta();
    rhoalpha = options->getRhoalpha();

    if(!options->isGlobalRho() && rhobeta0 > 1)
      rho = rhoalpha * rhobeta1 / (rhobeta0 - 1);//initialise at prior mean
    else if(!options->RhoFlatPrior() && !options->logRhoFlatPrior() )
      rho = rhoalpha / rhobeta;//initialise at prior mean for globalrho
    else rho = 2.0;//initialise at min value if flat prior
  
    
    // ** set up TuneRW object for global rho updates **
    NumberOfUpdates = 0;
    w = 1;
    step0 = 1.0; // sd of proposal distribution for log rho
    //need to choose sensible value for this initial RW sd
    step = step0;
    TuneRhoSampler.SetParameters( step0, 0.01, 10, 0.44);  
    
    // ** Open paramfile **
    if ( options->getIndAdmixHierIndicator()){
      Log.setDisplayMode(Quiet);
      if( strlen( options->getParameterFilename() ) ){
	outputstream.open( options->getParameterFilename(), ios::out );
	if( !outputstream )
	  {
	    Log.setDisplayMode(On);
	    Log << "ERROR: Couldn't open paramfile\n";
	    exit( 1 );
	  }
	else{
	  Log << "Writing population-level parameters to " << options->getParameterFilename() << "\n";
	  InitializeOutputFile(PopulationLabels);
	}
      }
      else{
	Log << "No paramfile given\n";
      }
    }
  }//end if Populations > 1
}

Latent::~Latent()
{
  delete[] poptheta;

#if POPADMIXSAMPLER == 1
  for(int i=0; i<options->getPopulations(); i++){
    delete DirParamArray[i];
  }
#elif POPADMIXSAMPLER == 2
  delete[] SumLocusAncestry;
#endif
}

void Latent::Update(int iteration, const IndividualCollection* const individuals, LogWriter &Log, bool anneal=false)
 {
   if( options->getPopulations() > 1 && individuals->getSize() > 1 &&
       options->getIndAdmixHierIndicator() ){

   // ** Sample for population admixture distribution Dirichlet parameters, alpha **
   // For a model in which the distribution of individual admixture in the population is a mixture
   // of components, we will have one Dirichlet parameter vector for each component, 
   // updated only from those individuals who belong to the component
   
#if POPADMIXSAMPLER == 1
      for( int j = 0; j < options->getPopulations(); j++ ){
         AlphaParameters[1] -= alpha[0][ j ];
         AlphaParameters[4] = individuals->getSumLogTheta(j);
         // elements of Dirichlet parameter vector are updated one at a time
         alpha[0][ j ] = DirParamArray[j]->Sample(&AlphaParameters, ddlogf);
         AlphaParameters[1] += alpha[0][ j ];
      }
#elif POPADMIXSAMPLER == 2
      //sample alpha conditional on individual admixture proportions
      PopAdmixSampler.Sample( obs, individuals->getSumLogTheta(), &alpha[0] );
#endif
   }
   // ** accumulate sum of Dirichlet parameter vector over iterations  **
   transform(alpha[0].begin(), alpha[0].end(), SumAlpha.begin(), SumAlpha.begin(), std::plus<double>());//SumAlpha += alpha[0];
   
   if( iteration < options->getBurnIn() && options->getPopulations() > 1
       && options->getNumberOfOutcomes() > 0 ){
     // accumulate ergodic average of population admixture, which is used to centre 
     // the values of individual admixture in the regression model
     double sum = accumulate(SumAlpha.begin(), SumAlpha.end(), 0.0);
     for( int j = 0; j < options->getPopulations(); j++ )poptheta[j] = SumAlpha[j] / sum;
   }
   
   if( iteration == options->getBurnIn() && options->getNumberOfOutcomes() > 0 
       && options->getPopulations() > 1) {
     Log.setDisplayMode(Off);
     Log << "Individual admixture centred in regression model around: ";
     for(int i = 0; i < options->getPopulations(); ++i)Log << poptheta[i] << "\t";
     Log << "\n";
     fill(SumAlpha.begin(), SumAlpha.end(), 0.0);
   }
   
   if( !anneal && iteration > options->getBurnIn() && options->getPopulations() > 1 ){
     // accumulate sum of log of sumintensities after burnin.
     if(options->isGlobalRho()) SumLogRho += log(rho);
     else SumLogRho += log(rhoalpha) - log(rhobeta);

   }
}
//end Update


void Latent::UpdateRhoWithRW(const IndividualCollection* const IC, Chromosome **C) {
  if( options->isGlobalRho() ){
    double rhoprop;
    double LogLikelihood = 0.0;
    double LogLikelihoodAtProposal = 0.0;
    double LogLikelihoodRatio = 0.0;
    double LogPriorRatio = 0.0;
    double LogAccProbRatio = 0.0;
    bool accept = false;

    NumberOfUpdates++;
    rhoprop = exp(gennor(log(rho), step)); // propose log rho from normal distribution with SD step
    
    //get log likelihood at current parameter values, annealed if this is an annealing run
    // TO DO - IC object should store log HMM likelihood summed over individuals
    // IC or individual objects should use annealed ProbGenotypes where necessary
    for(int i = 0; i < IC->getSize(); ++i) {
      Individual* ind = IC->getIndividual(i);
      ind->HMMIsBad(true);//to force HMM update
      LogLikelihood += ind->getLogLikelihood(options, C, false, true); // don't force update, store result if updated
      ind->HMMIsBad(true); // HMM probs overwritten by next indiv, but stored loglikelihood still ok
   }
    
    // set ancestry correlations using proposed value of sum-intensities 
    for( unsigned int j = 0; j < Loci->GetNumberOfChromosomes(); j++ ) C[j]->SetLociCorr(rhoprop);
    //get log HMM likelihood at proposal rho and current admixture proportions
    for(int i = 0; i < IC->getSize(); ++i) {
      Individual* ind = IC->getIndividual(i);
      LogLikelihoodAtProposal += ind->getLogLikelihood(options, C, true, false); // force update, do not store result 
      ind->HMMIsBad(true); // set HMM probs as bad but stored log-likelihood is still ok
      // line above should not be needed for a forced update with result not stored
    }
    LogLikelihoodRatio = LogLikelihoodAtProposal - LogLikelihood;

    //compute prior ratio
    LogPriorRatio = getGammaLogDensity(rhoalpha, rhobeta, rhoprop) - getGammaLogDensity(rhoalpha, rhobeta, rho);
    LogAccProbRatio = LogLikelihoodRatio + LogPriorRatio; 

    // generic Metropolis step
    if( LogAccProbRatio < 0 ) {
      if( log(myrand()) < LogAccProbRatio ) accept = true;
    } else accept = true;  
    
    if(accept) {
      rho = rhoprop;
      for(int i = 0; i < IC->getSize(); ++i){
	Individual* ind = IC->getIndividual(i);
	ind->storeLogLikelihood(false); // store log-likelihoods calculated at rhoprop, but do not set HMM probs as OK 
      }
    } else { 
      // restore ancestry correlations in Chromosomes using original value of sum-intensities
      for( unsigned int j = 0; j < Loci->GetNumberOfChromosomes(); j++ )
	C[j]->SetLociCorr(rho);
    } // stored loglikelihoods are still ok

    //update sampler object every w updates
    if( !( NumberOfUpdates % w ) ){
      step = TuneRhoSampler.UpdateStepSize( exp(LogAccProbRatio) );  
    }

  }//end if global rho model

  else if(IC->getSize()>1){//non global rho model
    //if a single individual, rhobeta is fixed
    // sample for location parameter of gamma distribution of sumintensities parameters 
    // in population 
    if( options->isRandomMatingModel() )
      rhobeta = gengam( IC->GetSumrho() + rhobeta1,
			2*rhoalpha * IC->getSize() + rhobeta0 );
    else
      rhobeta = gengam( IC->GetSumrho() + rhobeta1,
			rhoalpha* IC->getSize() + rhobeta0 );
  }
}
double Latent::getRhoSamplerAccRate()const{
  return TuneRhoSampler.getExpectedAcceptanceRate();
}
double Latent::getRhoSamplerStepsize()const{
  return step;
}

void Latent::InitializeOutputFile(const std::string* const PopulationLabels)
{
  // Header line of paramfile
  //Pop. Admixture
  for( int i = 0; i < options->getPopulations(); i++ ) {
    outputstream << "\""<<PopulationLabels[i] << "\" ";
  }
  //SumIntensities
  if( options->isGlobalRho() ) outputstream << "sumIntensities\t";
  else outputstream << "sumIntensities.mean\t";
  outputstream << endl;
}

void Latent::OutputErgodicAvg( int samples, std::ofstream *avgstream)
{
  for( int j = 0; j < options->getPopulations(); j++ ){
    avgstream->width(9);
    *avgstream << setprecision(6) << SumAlpha[j] / samples << "\t";
  }
  avgstream->width(9);
  *avgstream << setprecision(6) << exp(SumLogRho / samples) << "\t";
}
void Latent::OutputParams(ostream* out){
  for( int j = 0; j < options->getPopulations(); j++ ){
    out->width(9);
    (*out) << setprecision(6) << alpha[0][ j ] << "\t";
  }
  
  out->width(9);
  if( !options->isGlobalRho() )
    (*out) << setprecision(6) << rhoalpha / rhobeta  << "\t";
  else
    (*out) << setprecision(6) << rho << "\t";
}

void Latent::OutputParams(int iteration, LogWriter &Log){
  //output initial values to logfile
  if( iteration == -1 )
    {
      Log.setDisplayMode(Off);
      Log.setPrecision(6);
      for( int j = 0; j < options->getPopulations(); j++ ){
	//Log->width(9);
	Log << alpha[0][j];
      }
      
      if( !options->isGlobalRho() )
	Log << rhoalpha / rhobeta;
      else
	Log << rho;
    }
  //output to screen
  if( options->getDisplayLevel() > 2 )
    {
      OutputParams(&cout);
    }
  //Output to paramfile after BurnIn
  if( iteration > options->getBurnIn() ){
    OutputParams(&outputstream);
    outputstream << endl;
  }
}



const vector<double > &Latent::getalpha0()const{
  return alpha[0];
}
const std::vector<vector<double> > &Latent::getalpha()const{
  return alpha;
}

double Latent::getrhoalpha()const{
  return rhoalpha;
}
double Latent::getrhobeta()const{
  return rhobeta;
}
double Latent::getrho()const{
  return rho;
}
double Latent::getSumLogRho()const{
  return SumLogRho;
}
const double *Latent::getpoptheta()const{
  return poptheta;
}
#if POPADMIXSAMPLER == 1
// these 3 functions calculate log-likelihood and derivatives for adaptive rejection sampling of 
// Dirichlet population admixture parameters
double Latent::logf( double x, const void* const pars)
{
  const double* parameters = (const double*)pars;
  double f = parameters[0] * ( gsl_sf_lngamma( x + parameters[1] ) - gsl_sf_lngamma( x ) ) - 
    x * ( parameters[3] - parameters[4] ) + (parameters[2] - 1) * log(x);
  
  return(f);
}

double Latent::dlogf( double x, const void* const pars )
{
  const double* parameters = (const double*)pars;
  double f,x2,y1,y2;
  
  x2 = x + parameters[1];
   if(x < 0)cout<<"\nError in Latent::dlogf - arg x to ddigam is negative\n";   
  ddigam( &x , &y1 );
  if(x2 < 0)cout<<"\nError in Latent::dlogf - arg x2 to ddigam is negative\n";   
  ddigam( &x2 , &y2 );
  
  f = parameters[0] * ( y2 - y1 ) - ( parameters[3] - parameters[4] ) + (parameters[2] - 1)/x;
  
  return(f);
}

double Latent::ddlogf( double x, const void* const pars )
{
  const double* parameters = (const double*)pars;
  double f,x2,y1,y2;
  
  x2 = x + parameters[1];
  
  trigam( &x, &y1 );
  trigam( &x2, &y2 );
  
  f = parameters[0] * ( y2 - y1 ) - (parameters[2] - 1)/ (x*x);
  
  return(f);
}
#endif

#if POPADMIXSAMPLER == 2
float Latent::getEtaSamplerAcceptanceRate()const{
  return PopAdmixSampler.getEtaExpectedAcceptanceRate();
}
float Latent::getEtaSamplerStepsize()const{
  return PopAdmixSampler.getEtaStepSize();
}

#endif

#if POPADMIXSAMPLER == 3
float Latent::getAlphaSamplerAcceptanceRate()const{
  return AlphaSampler.getAcceptanceRate();
}
float Latent::getAlphaSamplerStepsize()const{
  return AlphaSampler.getStepsize();
}


#endif
