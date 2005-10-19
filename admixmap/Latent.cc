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

Latent::Latent( AdmixOptions * op, Genome *loci, LogWriter *l)
{
  options = 0;
  rho = 0.0;
  rhoalpha = 0.0;
  rhobeta = 0.0;
  SumLogRho = 0.0;
  options = op;
  Loci = loci;
  Log = l;
  poptheta = 0;
#if POPADMIXSAMPLER == 2
  mu = 0;
  SumLocusAncestry = 0;
#elif POPADMIXSAMPLER == 3
  logalpha = 0;
  initialAlphaStepsize = 0.05;//need a way of setting this without recompiling, or a sensible fixed value
  targetAlphaAcceptRate = 0.44;//need to choose suitable value for this
#endif
}

void Latent::Initialise(int Numindividuals, std::string *PopulationLabels){
  // ** Initialise population admixture distribution Dirichlet parameters alpha **
  int K = options->getPopulations();
  //ergodic average of population admixture, which is used to centre 
  // the values of individual admixture in the regression model  
  poptheta = new double[ K ];
  for( int i = 0; i < K; i++ )poptheta[i] = 0.0;
  alpha = options->getAndCheckInitAlpha(Log);
  SumAlpha.resize( K );
  if(!options->getIndAdmixHierIndicator())  copy(alpha[0].begin(), alpha[0].end(), SumAlpha.begin());

  if(K > 1){
    double alphapriormean = options->getAlphamean();
    double alphapriorvar = options->getAlphavar();
    Log->logmsg(true, "Gamma prior on population admixture Dirichlet parameters with mean ");
    Log->logmsg(true, alphapriormean);
    Log->logmsg(true, " and variance ");
    Log->logmsg(true, alphapriorvar);
    Log->logmsg(true, "\n");
    
    // ** set up sampler for alpha **
    
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
    eta = accumulate(alpha[0].begin(), alpha[0].end(), 0.0, std::plus<double>());//eta = sum of alpha[0]
    mu = new double[ K ];
    PopAdmixSampler.SetSize( Numindividuals, K );
    for( int i = 0; i < K; i++ ){
      mu[i] = alpha[0][i]/eta;
    }
    if( options->isRandomMatingModel() ){
      obs = 2 * Numindividuals;
    } else {
      obs = Numindividuals;
    }
    SumLocusAncestry = new int[Numindividuals*K];
#elif POPADMIXSAMPLER == 3
    logalpha = new double[K];
    transform(alpha[0].begin(), alpha[0].end(), logalpha, xlog);//logalpha = log(alpha)
    
    AlphaArgs = new double*[4];
    AlphaArgs[0] = new double[K];
    for(unsigned i = 1; i < 4;++i)AlphaArgs[i] = new double[1];
    //elem 0 is sum of log admixture props
    AlphaArgs[1][0] = (double)Numindividuals;// elem 1 is num individuals/gametes
    if( options->isRandomMatingModel() )AlphaArgs[1][0] *= 2.0;
    AlphaArgs[2][0] = alphapriormean*alphapriormean / alphapriorvar;//params of gamma prior
    AlphaArgs[3][0] = alphapriormean / alphapriorvar;
    
    AlphaSampler.SetDimensions(K, initialAlphaStepsize, 0.01, 10.0, 20, targetAlphaAcceptRate, findE, gradE);
#endif
    
    // ** Initialise sum-of-intensities parameter rho and the parameters of its prior, rhoalpha and rhobeta **
    rho = options->getRhoalpha()/options->getRhobeta();
    
    if( options->RhoFlatPrior() ){
      rhoalpha = 1.0;
      rhobeta = 0.0;
      Log->logmsg(true,"Flat prior on sumintensities.\n");
    }
    else if( options->logRhoFlatPrior() ){
      rhoalpha = 0.0;
      rhobeta = 0.0;
      Log->logmsg(true,"Flat prior on log sumintensities.\n");
    }
    else{
      rhoalpha = options->getRhoalpha();
      rhobeta = options->getRhobeta();
      Log->logmsg(false,"Gamma prior on sum-of-intensities with shape parameter: ");
      Log->logmsg(false, rhoalpha); Log->logmsg(false,"\n");
      Log->logmsg(false," and rate (1 / location) parameter: ");
      Log->logmsg(false, rhobeta); Log->logmsg(false,"\n");
    }
    rhobeta0 = 1;
    rhobeta1 = 1;
    
    // ** set up TuneRW object for global rho updates **
    NumberOfUpdates = 0;
    w = 1;
    step0 = 1.0; // sd of proposal distribution for log rho
    //need to choose sensible value for this initial RW sd
    step = step0;
    TuneRhoSampler.SetParameters( step0, 0.01, 10, 0.44);  
    
    // ** Open paramfile **
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
	  InitializeOutputFile(PopulationLabels);
	}
      }
      else{
	Log->logmsg(true,"No paramfile given\n");
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
  delete[] mu;
  delete[] SumLocusAncestry;
#elif POPADMIXSAMPLER == 3
  for(unsigned i = 0; i < 4;++i)delete[] AlphaArgs[i];
  delete[] AlphaArgs;
  delete[] logalpha;
#endif
}

void Latent::Update(int iteration, IndividualCollection *individuals)
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
//       if((iteration %2)){//even-numbered iterations
// 	//sample mu conditional on sum of ancestry states where jump indicator==1
// 	unsigned I = individuals->getSize();
// 	unsigned K = options->getPopulations();
// 	for(unsigned i = 0; i < I; ++i){
// 	  int *SLA_ind = individuals->getIndividual(i)->getSumLocusAncestry();
// 	  for(unsigned k = 0; k < K; ++k)
// 	    SumLocusAncestry[k*I + i] = SLA_ind[k];	  
// 	}
// 	PopAdmixSampler.Sample2( obs, individuals->getSumLogTheta(), &eta, mu, SumLocusAncestry );
//       }
//       else//odd-numbered iterations; sample mu conditional on individual admixture proportions
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

#endif
      
      // ** accumulate sum of Dirichlet parameter vector over iterations  **
      transform(alpha[0].begin(), alpha[0].end(), SumAlpha.begin(), SumAlpha.begin(), std::plus<double>());//SumAlpha += alpha[0];
   }
   
   if( iteration < options->getBurnIn() && options->getPopulations() > 1
       && options->getNumberOfOutcomes() > 0 ){
     // accumulate ergodic average of population admixture, which is used to centre 
     // the values of individual admixture in the regression model
     double sum = accumulate(SumAlpha.begin(), SumAlpha.end(), 0.0);
     for( int j = 0; j < options->getPopulations(); j++ )poptheta[j] = SumAlpha[j] / sum;
   }
   
   if( iteration == options->getBurnIn() && options->getNumberOfOutcomes() > 0 ){
     Log->write("Individual admixture centred in regression model around: ");
     Log->write( poptheta, options->getPopulations());
     Log->write("\n");
     
     fill(SumAlpha.begin(), SumAlpha.end(), 0.0);
   }
   
   if( iteration > options->getBurnIn() ){
     // accumulate sum of log of rho parameters after burnin.
     if( options->getPopulations() > 1 ){
       SumLogRho += log(rho);
     }
   }
}
//end Update


void Latent::UpdateRhoWithRW(IndividualCollection *IC, Chromosome **C, double LogL){

  if( options->isGlobalRho() ){
    double rhoprop;
    double LogLikelihoodRatio = 0.0;
    double LogPriorRatio = 0.0;
    double LogAccProb;

    NumberOfUpdates++;
    rhoprop = exp(gennor(log(rho), step)); // propose log rho from normal distribution with SD step
    
    //get log likelihood at current parameter values
    //    for(int i = 0; i < IC->getSize(); ++i)LogLikelihoodRatio -= IC->getIndividual(i)->getLogLikelihood(options, C);
    LogLikelihoodRatio -= LogL;
    
    //get log likelihood at proposal rho and current admixture proportions
    for( unsigned int j = 0; j < Loci->GetNumberOfChromosomes(); j++ ) C[j]->SetLociCorr(rhoprop);
    for(int i = 0; i < IC->getSize(); ++i){
      Individual* ind = IC->getIndividual(i);
      std::vector<double> rhovec(2, rhoprop);//note that the values in here are irrelevant as Chromosome ignores them in a globalrho model
      LogLikelihoodRatio += ind->getLogLikelihood(options, C, ind->getAdmixtureProps(), ind->getAdmixtureProps(),
						  rhovec, rhovec, false);

    }
    
    //compute prior ratio
    LogPriorRatio = getGammaLogDensity(rhoalpha, rhobeta, rhoprop) - getGammaLogDensity(rhoalpha, rhobeta, rho);
    
    LogAccProb = 0.0;
    if(LogLikelihoodRatio + LogPriorRatio < 0.0) 
      LogAccProb = LogLikelihoodRatio + LogPriorRatio; 
    //accept/reject proposal
    if(log( myrand() ) < LogAccProb){//accept
      rho = rhoprop;
      for(int i = 0; i < IC->getSize(); ++i)
	IC->getIndividual(i)->HMMIsBad(true);//rho has changed so current stored loglikelihood is invalid and HMMs need to be updated
    }
    else//reject
      // restore f in Chromosomes
      for( unsigned int j = 0; j < Loci->GetNumberOfChromosomes(); j++ )
	C[j]->SetLociCorr(rho);

    //update sampler object every w updates
    if( !( NumberOfUpdates % w ) ){
      step = TuneRhoSampler.UpdateStepSize( exp(LogAccProb) );
    }

  }//end if global rho model

  else{//non global rho model
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
double Latent::getRhoSamplerAccRate(){
  return TuneRhoSampler.getExpectedAcceptanceRate();
}
double Latent::getRhoSamplerStepsize(){
  return step;
}

void Latent::InitializeOutputFile(std::string *PopulationLabels)
{
  // Header line of paramfile

  if( options->getIndAdmixHierIndicator() ){

    //Pop. Admixture
    for( int i = 0; i < options->getPopulations(); i++ ){
      outputstream << "\""<<PopulationLabels[i] << "\" ";
    }
    //SumIntensities
    if( options->isGlobalRho() )
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
  *avgstream << setprecision(6) << exp(SumLogRho / samples) << " ";
}

void Latent::OutputParams(int iteration){
  //output to logfile for first iteration or every iteration if cout = 0
  if( !options->useCOUT() || iteration == 0 )
    {
      for( int j = 0; j < options->getPopulations(); j++ ){
	Log->width(9);
	Log->write(alpha[0][j], 6);
      }
       if( !options->isGlobalRho() )
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
      if( !options->isGlobalRho() )
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
    if( !options->isGlobalRho() )
      (outputstream) << setprecision(6) << rhobeta << " ";
    else
      (outputstream) << setprecision(6) << rho << " ";

    outputstream << endl;
  }

}



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
double Latent::getSumLogRho(){
  return SumLogRho;
}
const double *Latent::getpoptheta(){
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
float Latent::getEtaSamplerAcceptanceRate(){
  return PopAdmixSampler.getEtaExpectedAcceptanceRate();
}
float Latent::getEtaSamplerStepsize(){
  return PopAdmixSampler.getEtaStepSize();
}
// float Latent::getMuSamplerAcceptanceRate(){
//   return PopAdmixSampler.getMuExpectedAcceptanceRate();
// }
// float Latent::getMuSamplerStepsize(){
//   return PopAdmixSampler.getMuStepSize();
//}
#endif

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
