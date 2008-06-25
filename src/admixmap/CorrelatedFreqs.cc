/** 
 *   ADMIXMAP
 *   CorrelatedFreqs.cc
 *   Class to hold and update allele frequencies in a correlated frequencies model
 *   Copyright (c) 2007 David O'Donnell and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "CorrelatedFreqs.h"
#include "AdmixOptions.h"
#include "InputAdmixData.h"
#include <iomanip>
#include "bclib/misc.h"
#include <string.h>
using bclib::Rand;

CorrelatedFreqs::CorrelatedFreqs(){
  eta = 0.0;
  psi = 0.0;
  tau = 0.0; 
  SumEta = 0.0;
  psi0 = 0.0;

  w = 1; // frequency of tuning sampler
  NumberOfEtaUpdates  = 0;
  etastep = 0;
}

void CorrelatedFreqs::Initialise(Options* const a_options, InputData* const a_data, 
				 Genome *pLoci, bclib::LogWriter &Log, bool MAP){

  AdmixOptions const* options = (AdmixOptions*)a_options;
  InputAdmixData const* data = (InputAdmixData*)a_data;

  AdmixFreqs::Initialise(a_options, a_data, pLoci, Log, MAP);

  // ** settings for sampling of dispersion parameter **
  
  //NOTE: if DispersionSampler is used to sample eta in both historic and correlated allele freq models, we don't need to
  //store psi and tau here, just pass the values to the sampler
  
  // ** set eta priors **
  //rate parameter
  tau = options->getEtaMean() / options->getEtaVar();
  //shape parameter
  psi = tau *options->getEtaMean();
  
  //w = 1;
  
  muSampler = new bclib::MuSampler[NumberOfCompositeLoci];
  
  for(unsigned i = 0; i < NumberOfCompositeLoci; ++i)
    muSampler[i].setDimensions(Populations, Loci->GetNumberOfStates(i), 0.002, 0.00001, 10.0, 0.9);
  
  if(ETASAMPLER==2) {
    
    int* NumStates;
    NumStates = new int[NumberOfCompositeLoci];
    for(unsigned i = 0; i < NumberOfCompositeLoci; ++i) {
      NumStates[i] = Loci->GetNumberOfStates(i);
    }
    bclib::DispersionSampler::setDimensions(NumberOfCompositeLoci, Populations, NumStates);
    delete[] NumStates;
    
    EtaSampler.Initialise( 0.05, 0.01, 1000.0, 0.9);
    EtaSampler.setEtaPrior(psi, tau); 
    
  }
  
  
  //Initialise eta at its prior expectation
  eta = psi/tau;
  SumEta = 0.0;
  //Rescale priorallelefreqs so the columns sum to eta 
  for(unsigned j = 0; j < NumberOfCompositeLoci; j++ ){
    int NumberOfStates = Loci->GetNumberOfStates(j);
    double sum = 0.0;
    for(int s = 0; s < NumberOfStates; ++s)sum+= PriorParams[j][s];
    for(int s = 0; s < NumberOfStates; ++s)PriorParams[j][s]*= eta / sum;
  }
  
  
  if(ETASAMPLER==1) {
    etastep0 = 0.1; // sd of proposal distribution for log eta
    etastep = etastep0;
    NumberOfEtaUpdates = 0;
    
    TuneEtaSampler.SetParameters( etastep0, 0.01, 10, 0.44 );
    
    // ** Open output file for eta **
    if ( options->getIndAdmixHierIndicator()){
      if (strlen( options->getEtaOutputFilename() ) ){
	InitializeEtaOutputFile(options->getEtaOutputFilename(), data->GetPopLabels(), Log); 
      }
      else{
	Log << "Not writing dispersion parameter to file\n";
	//exit(1);
      }
    }
  }
}

void CorrelatedFreqs::PrintPrior(const Vector_s& PopLabels, bclib::LogWriter& Log)const{
  Log << "Gamma prior on allele frequency dispersion parameter with mean and variance:\n"
      << psi/tau << "  " << psi/(tau*tau) << "\n";

}

/// samples allele frequencies and prior parameters.
void CorrelatedFreqs::Update(IndividualCollection*IC , bool afterBurnIn, double coolness){
  // Sample prior parameters
  
  SampleDirichletParams();

  AdmixFreqs::Update(IC, afterBurnIn, coolness);

  // Sample for allele frequency dispersion parameters, eta, conditional on allelefreqs using
  // Metropolis random-walk.

  NumberOfEtaUpdates++;
  if(ETASAMPLER == 1) { // random walk metropolis
    SampleEtaWithRandomWalk(afterBurnIn);
  } // otherwise eta will be sampled in sampleDirichletParams
}

void CorrelatedFreqs::InitializeEtaOutputFile(const char* filename, const Vector_s& PopulationLabels, bclib::LogWriter &Log)
{
  Log.setDisplayMode(bclib::On);
  outputstream.open(filename );
  if( !outputstream )
    {
      Log <<  "ERROR: Couldn't open dispparamfile\n";
      exit( 1 );
    }
  else{
    Log << "Writing dispersion parameter to " << filename << "\n";
    //Dispersion parameters (eta)
    outputstream.delimit(false);

    outputstream << "\"eta\""; 

    outputstream.delimit(true);
    outputstream << bclib::newline;
  }
}

/**
 * Given the number of ancestral populations, sets default values for
 * prior allele freq params.
 */
void CorrelatedFreqs::SetDefaultPriorParams(int i, double defaultpriorparams){
  int NumberOfStates = Loci->GetNumberOfStates(i);

  PriorParams[i] = new double[NumberOfStates];
  // reference prior on allele freqs: all elements of parameter vector set to 0.5
  // this is unrealistic for large haplotypes - should set all elements to sum to 1
  fill(PriorParams[i], PriorParams[i] + NumberOfStates, 0.5);

}

void CorrelatedFreqs::LoadAlleleFreqs(const Matrix_s& New, int i, unsigned row0, bool ){
  // set size of allele freqs array for this locus
  // Freqs array has NumberOfStates - 1 elements for each population
  unsigned NumberOfStates = Loci->GetNumberOfStates(i);

    double sumalpha;
    // allele frequencies are initialised as expectations over the Dirichlet prior distribution, 
    // by dividing each prior parameter by the sum of the parameters.
    for( unsigned j = 0; j < Populations; j++ ){
      //vector<double> NewCol = New.getCol(j);
      //sumalpha = accumulate( NewCol.begin(), NewCol.end(), 0.0, std::plus<double>() );
      sumalpha = 0.0;
      for( unsigned k = 0; k < NumberOfStates ; k++ )sumalpha+= convertValueFromFile(New[row0+k][j+1]);//New.get(row0+k, j+1);
      for( unsigned k = 0; k < NumberOfStates ; k++ )
	Freqs[i][ k + j*NumberOfStates ] = ( convertValueFromFile(New[row0+k][1] )
	  /*New.get( row0+k, +1)*/ ) / sumalpha;
    }

    PriorParams[i] = new double[NumberOfStates];

    for(unsigned row = 0; row < NumberOfStates; ++row)
      PriorParams[i][row] = convertValueFromFile(New[row0+row][1]);//New.get(row0+row, col+1); 

}

/** samples allele/hap freqs at i th composite locus as a conjugate Dirichlet update
 and stores result in array Freqs 
*/
void CorrelatedFreqs::SampleAlleleFreqs(int i, double coolness)
{
  unsigned NumStates = Loci->GetNumberOfStates(i);
  double* temp = new double[NumStates];
  double *freqs = new double[NumStates];
  
  //if there is, the Dirichlet params are common across populations
  for( unsigned j = 0; j < Populations; j++ ){

    // to flatten likelihood when annealing, multiply realized allele counts by coolness
    for(unsigned s = 0; s < NumStates; ++s)
      temp[s] = PriorParams[i][s] + coolness*AlleleCounts[i][s*Populations +j];

    Rand::gendirichlet(NumStates, temp, Freqs[i]+j*NumStates);
    for(unsigned s = 0; s < NumStates; ++s){
	if(Freqs[i][j*NumStates+s]==0.0) Freqs[i][j*NumStates+s] = 0.000001;
	if(Freqs[i][j*NumStates+s]==1.0) Freqs[i][j*NumStates+s] = 0.999999;
    }
  }
  delete[] freqs;  
  delete[] temp;
}

///sampling of Dirichlet parameters for allelefreqs in historic or correlated allele freq model
void CorrelatedFreqs::SampleDirichletParams() {
  for(unsigned i = 0; i < NumberOfCompositeLoci; ++i){
    muSampler[i].Sample(PriorParams[i], eta, AlleleCounts[i]);
  }

  if(ETASAMPLER == 2) { // hamiltonian sampler
    for(unsigned i = 0; i < NumberOfCompositeLoci; ++i) {
      EtaSampler.addAlphas(i, PriorParams[i]);
      EtaSampler.addCounts(i, AlleleCounts[i]);
    }
    eta = EtaSampler.Sample();
    SumEta += eta;
    //rescale PriorParams so they sum to eta
    for( unsigned i = 0; i < NumberOfCompositeLoci; ++i ) {
      int NumberOfStates = Loci->GetNumberOfStates(i);
      double sum = accumulate(PriorParams[i], PriorParams[i]+NumberOfStates, 0.0, 
			      std::plus<double>());
      sum = eta / sum;
      for(int t = 0; t < NumberOfStates; ++t) {
	PriorParams[i][t] *= sum;
      }
    }
  }
}

void CorrelatedFreqs::SampleEtaWithRandomWalk(bool updateSumEta) {
  double etanew, LogPostRatio = 0.0, LogLikelihoodRatio = 0.0, LogPriorRatio = 0.0, AccProb = 0.0;
  double Denom = 0.0;
  double mineta = 0;
  vector< vector<double> > munew;
  // propose etanew from truncated log-normal distribution.
  do{
    etanew = exp( Rand::gennor( log( eta ), etastep ) );
  }while( etanew > 5000.0 );
  // Prior log-odds ratio (proposal ratio cancels with a part of the prior ratio)   
  LogPriorRatio = ( psi - 1 ) * (log(etanew) - log(eta)) - tau * ( etanew - eta );
  // Log-likelihood ratio; numerator of integrating constant
  LogLikelihoodRatio += 2 * NumberOfCompositeLoci * ( gsl_sf_lngamma( etanew ) - gsl_sf_lngamma( eta ) );
  for(unsigned j = 0; j < NumberOfCompositeLoci; j++ ){
    std::vector<double> mu = GetPriorAlleleFreqs(j,0);
    vector<double>::const_iterator it = min_element(mu.begin(), mu.end());
    
    //mineta is a lower bound for proposal etanew
    if( mineta < 0.1 * eta / *it )
      mineta = 0.1 * eta / *it;
    
    munew.push_back( mu );
    //rescale munew so munew sums to etanew
    for( unsigned l = 0; l < mu.size(); l++ )
      munew[j][l] *= etanew / eta;
    
    const double *SumLogFreqs = GetStatsForEta(j);
    for( unsigned l = 0; (int)l < (*Loci)(j)->GetNumberOfStates(); l++ ){
      // Denominator of integrating constant
      Denom += 2*(gsl_sf_lngamma( mu[l] ) - gsl_sf_lngamma( munew[j][l] ));
      // SumLogFreqs = log phi_1 + log phi_2
      Denom += (munew[j][l] - mu[l])*SumLogFreqs[l];
    }
    delete[] SumLogFreqs;
  }
  LogPostRatio = LogPriorRatio + LogLikelihoodRatio + Denom;
  // Log acceptance probability = Log posterior ratio since the
  // proposal ratio (log-normal) cancels with prior.
  if(LogPostRatio > 0.0)LogPostRatio = 0.0;
  AccProb = 0.0; 
  if(mineta < etanew )AccProb = bclib::xexp(LogPostRatio);
  
#ifdef DEBUGETA
  cout<< "eta = "<< eta<<" eta* = "<<etanew<< " stepsize = "<<etastep<<endl;
  cout << "LogPriorRatio= "<< LogPriorRatio <<" LogLikelihoodRatio= " << LogLikelihoodRatio <<endl
       << "Denom = "<< Denom<< " AccProb= exp("<<LogPostRatio<<") = "<< AccProb<<endl<<endl;
#endif
  
  // Acceptance test.
  if( log( Rand::myrand() ) < LogPostRatio && mineta < etanew ){
    eta = etanew;
    UpdatePriorAlleleFreqs( munew );
  }
  
  if( !( NumberOfEtaUpdates % w ) ){
    etastep = TuneEtaSampler.UpdateStepSize( AccProb );
  }
  
  if( updateSumEta )
    SumEta+=eta;
}

const double *CorrelatedFreqs::GetStatsForEta( int locus )const
{
  int NumberOfStates = Loci->GetNumberOfStates(locus);
  double *stats = new double[ NumberOfStates ];
  fill(stats, stats+NumberOfStates, 0.0);

  for(unsigned k = 0; k < Populations; ++k){
    for( int i = 0; i < NumberOfStates; i++ ){
      stats[ i ] += log( Freqs[locus][ i + k*NumberOfStates ] );
    }
  }

  return stats;
}

///sets PriorParams to sum to eta after sampling of eta
void CorrelatedFreqs::UpdatePriorAlleleFreqs(const vector<vector<double> >& mu)
{
  for( unsigned i = 0; i < NumberOfCompositeLoci; i++ ){
    for(int h = 0; h < Loci->GetNumberOfStates(i); ++h)
      PriorParams[i][h] = mu[i][h];
  }
}

void CorrelatedFreqs::OutputErgodicAvg( int samples, std::ofstream *avgstream)const
{
  avgstream->width(9);
  *avgstream << setprecision(6) << SumEta / samples << "\t";
}

void CorrelatedFreqs::OutputParams(){
  if(outputstream.is_open()){
    OutputParams(outputstream);
    outputstream << bclib::newline;
  }
}

void CorrelatedFreqs::OutputParams(bclib::Delimitedostream& os)const{
  cout.width(9);
  os << setprecision(6) << eta ;

}

void CorrelatedFreqs::PrintAcceptanceRates(bclib::LogWriter& Log){
  //TODO maybe write averages over loci here?
  Log<< "Expected acceptance rates in sampler for allele frequency proportion parameters: \n";
  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++){
    if(Loci->GetNumberOfStates(i)>2)
      Log << muSampler[i].getAcceptanceRate() << " ";
  }
  Log<< "Expected acceptance rate in sampler for allele frequency dispersion parameter: \n";
  Log << getEtaSamplerAcceptanceRate()
      << "\nwith final step size of \n";

  for(unsigned int i = 0; i < NumberOfCompositeLoci; i++){
    if(Loci->GetNumberOfStates(i)>2)
      Log << muSampler[i].getStepsize() << " " ;
  }
  Log <<  getEtaSamplerStepsize() << "\n" ;
}

float CorrelatedFreqs::getEtaSamplerAcceptanceRate()const{
  if(ETASAMPLER == 1) {
    return TuneEtaSampler.getExpectedAcceptanceRate();
  } else {
    return EtaSampler.getAcceptanceRate();
  }
}

float CorrelatedFreqs::getEtaSamplerStepsize()const{
  if(ETASAMPLER==1) {
    return etastep;
  } else {
    return EtaSampler.getStepsize();
  }
}

