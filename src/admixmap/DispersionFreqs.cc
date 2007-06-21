/** 
 *   ADMIXMAP
 *   DispersionFreqs.cc
 *   Class to hold and update allele frequencies, their prior parameters, allele counts and sums. Also holds and updates dispersion
 *   parameter eta and its prior parameters, for a dispersion model. Also computes Fst if required.
 *   Copyright (c) 2005 - 2007 David O'Donnell, Clive Hoggart and Paul McKeigue
 *  
 * This program is free software distributed WITHOUT ANY WARRANTY. 
 * You can redistribute it and/or modify it under the terms of the GNU General Public License, 
 * version 2 or later, as published by the Free Software Foundation. 
 * See the file COPYING for details.
 * 
 */
#include "DispersionFreqs.h"
#include "AdmixOptions.h"
#include "AdmixFilenames.h"
#include "InputAdmixData.h"
#include "bcppcl/AdaptiveRejection.h"
#include "bcppcl/misc.h"
#include "bcppcl/linalg.h"
#include "bcppcl/MuSampler.h"
#include <iomanip>
#include <math.h>
#include <numeric>


//#define DEBUGETA 1

static double convertValueFromFile(const string s){
  double d = atof(s.c_str());
  if(d < 0.000001) d = 0.000001;//values must be strictly positive
  return d;
}

DispersionFreqs::DispersionFreqs(){
  eta = 0;
  psi = 0;
  tau = 0; 
  SumEta = 0;
  psi0 = 0.0;
  IsHistoricAlleleFreq = false;
  CorrelatedAlleleFreqs = false;
  HistoricAlleleFreqs = 0;
  HistoricAlleleCounts = 0;

  Fst = 0;
  SumFst = 0;
  calculateFST = false;
  MuProposal = 0;
  TuneEtaSampler = 0;
  w = 1; // frequency of tuning sampler
  NumberOfEtaUpdates  = 0;
  etastep = 0;
}

DispersionFreqs::~DispersionFreqs(){

  if(IsHistoricAlleleFreq){
    for(int i = 0; i < NumberOfCompositeLoci; ++i){
      delete[] HistoricAlleleCounts[i];
      delete[] HistoricAlleleFreqs[i];
    }
    delete[] HistoricAlleleCounts;
    delete[] HistoricAlleleFreqs;
  }

  delete[] Fst;
  delete[] SumFst;
  delete[] psi;
  delete[] tau;
  delete[] SumEta;
  delete[] MuProposal;

#if ETASAMPLER == 1 
if( IsHistoricAlleleFreq || CorrelatedAlleleFreqs ) {
  delete[] etastep;
  delete[] TuneEtaSampler;
}
#elif ETASAMPLER == 2
if( IsHistoricAlleleFreq || CorrelatedAlleleFreqs ) {
  delete[] etastep;
  delete[] EtaSampler;
}
#endif

}

// ************** Initialisation and loading of data  *******************
void DispersionFreqs::Initialise(AdmixOptions* const options, InputAdmixData* const data, 
				 Genome *pLoci, LogWriter &Log, bool MAP){
  //initialise Freqs, PriorAlleleFreqs, HistoricAlleleFreqs etc
  Loci = pLoci;
  Populations = options->getPopulations();
  NumberOfCompositeLoci = Loci->GetNumberOfCompositeLoci();
  //set model indicators
  hapmixmodel = options->getHapMixModelIndicator();
  RandomAlleleFreqs = !options->getFixedAlleleFreqs();
  CorrelatedAlleleFreqs = options->getCorrelatedAlleleFreqs();

  
  LoadAlleleFreqs(options, data, Log);
  Log.setDisplayMode(On);
  //open allelefreqoutputfile
  if(IsRandom() ){
    const char* s = options->getAlleleFreqOutputFilename();
    if(strlen(s))
      OpenOutputFile(s);
  }
  
  // set which sampler will be used for allele freqs
  // current version uses conjugate sampler if annealing without thermo integration
  if( (options->getThermoIndicator() && !options->getTestOneIndivIndicator()) ||
      //using default allele freqs or CAF model
      ( !strlen(options->getAlleleFreqFilename()) &&
	!strlen(options->getHistoricalAlleleFreqFilename()) && 
	!strlen(options->getPriorAlleleFreqFilename()) && 
	!options->getCorrelatedAlleleFreqs() ) ) {
    FREQSAMPLER = FREQ_HAMILTONIAN_SAMPLER;
  } else {
    FREQSAMPLER = FREQ_CONJUGATE_SAMPLER;
  }
  
  if(MAP)
    setAlleleFreqsMAP();    
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
    if(RandomAlleleFreqs){
      if (FREQSAMPLER==FREQ_HAMILTONIAN_SAMPLER){
	//set up samplers for allelefreqs
	FreqSampler.push_back(new AlleleFreqSampler(Loci->GetNumberOfStates(i), Populations, 
						    PriorParams[i], false));
      }
    }
    //set AlleleProbs pointers in CompositeLocus objects to point to Freqs
    //initialise AlleleProbsMAP pointer to 0
    //allocate HapPairProbs and calculate them using AlleleProbs
    (*Loci)(i)->InitialiseHapPairProbs(Freqs[i], false);
    //If using Chib algorithm, allocate HapPairProbsMAP and copy values in HapPairProbs
    if(MAP) {
      //set AlleleProbsMAP in Composite Locus to point to AlleleFreqsMAP
      (*Loci)(i)->setAlleleProbsMAP(AlleleFreqsMAP[i]);
      (*Loci)(i)->InitialiseHapPairProbsMAP();
    }
  }//end comp locus loop
  
  // ** settings for sampling of dispersion parameter **
  if( IsHistoricAlleleFreq || CorrelatedAlleleFreqs) {
    unsigned dim = 1;
    if(IsHistoricAlleleFreq) {
      dim = Populations;
    } 
    // ** dispersion parameter(s) and priors **
    eta = new double[ dim ];//dispersion parameters
    psi = new double[ dim ];//gamma prior shape parameter
    tau = new double[ dim ];//gamma prior rate parameter
    SumEta = new double[ dim ];//running sums
    
    //NOTE: if DispersionSampler is used to sample eta in both historic and correlated allele freq models, we don't need to
    //store psi and tau here, just pass the values to the sampler
    
    // ** set eta priors **
    if( strlen(options->getEtaPriorFilename()) ){
      //specified by user in file
      Log << "Loading gamma prior parameters for allele frequency dispersion from "
	  << options->getEtaPriorFilename() << ".\n";
      const DataMatrix& etaprior = data->getEtaPriorMatrix();
      
      for( unsigned k = 0; k < dim; k++ ){
	psi[k] = etaprior.get( k, 0 );
	tau[k] = etaprior.get( k, 1 );
      }
    }
    else{//specified by user, otherwise default of psi = 3.0, tau = 0.01
      double rate = options->getEtaMean() / options->getEtaVar();
      double shape = rate *options->getEtaMean();
      
      fill(psi, psi + dim, shape);
      fill(tau, tau + dim, rate);
    }
    
    //w = 1;
    
    if(IsHistoricAlleleFreq) {
      // ** settings for sampling of proportion vector mu of prior on allele freqs
      muSampler = new MuSampler[dim*NumberOfCompositeLoci]; // 
      for(int i = 0; i < NumberOfCompositeLoci; ++i)
	for(int k = 0; k < Populations; ++k)
	  muSampler[i*Populations+k].setDimensions(2, Loci->GetNumberOfStates(i), 0.001, 0.000001, 10.0, 0.9);
    } else {//correlated allele freq model
      muSampler = new MuSampler[NumberOfCompositeLoci];
      
      for(int i = 0; i < NumberOfCompositeLoci; ++i)
	muSampler[i].setDimensions(Populations, Loci->GetNumberOfStates(i), 0.002, 0.00001, 10.0, 0.9);
      
      if(ETASAMPLER==2) {
	EtaSampler = new DispersionSampler[dim];
	
	int* NumStates;
	NumStates = new int[NumberOfCompositeLoci];
	for(int i = 0; i < NumberOfCompositeLoci; ++i) {
	  NumStates[i] = Loci->GetNumberOfStates(i);
	}
	DispersionSampler::setDimensions(NumberOfCompositeLoci, Populations, NumStates);
	delete[] NumStates;
	
	EtaSampler[0].Initialise( 0.05, 0.01, 1000.0, 0.9);
	EtaSampler[0].setEtaPrior(psi[0], tau[0]); 
	
      }
    }
    
    for( unsigned k = 0; k < dim; k++ ){
      //Initialise eta at its prior expectation
      eta[k] = psi[k]/tau[k];
      //Rescale priorallelefreqs so the columns sum to eta 
      for(int j = 0; j < NumberOfCompositeLoci; j++ ){
	int NumberOfStates = Loci->GetNumberOfStates(j);
	double sum = 0.0;
	for(int s = 0; s < NumberOfStates; ++s)sum+= PriorParams[j][k*NumberOfStates+s];
	for(int s = 0; s < NumberOfStates; ++s)PriorParams[j][k*NumberOfStates+s]*= eta[k] / sum;
      }
    }
    
    if(ETASAMPLER==1) {
      etastep0 = 0.1; // sd of proposal distribution for log eta
      etastep = new double[ dim ];
      for(unsigned k = 0; k < dim; ++k) etastep[k] = etastep0;
      NumberOfEtaUpdates = 0;
      TuneEtaSampler = new StepSizeTuner[ dim ];
      for( unsigned k = 0; k < dim; k++ )
	TuneEtaSampler[k].SetParameters( etastep0, 0.01, 10, 0.44 );
    }
    
    // ** Open output file for eta **
    if ( options->getIndAdmixHierIndicator()){
      if (strlen( options->getEtaOutputFilename() ) ){
	InitializeEtaOutputFile(options, data->GetPopLabels(), Log); 
      }
      else{
	Log << "No dispparamfile given\n";
	//exit(1);
      }
    }
    // ** open fst output file if specified **
    if( options->getOutputFST() ){
      OpenFSTFile(options->getResultsDir(),Log);
    }
  } //end if dispersion parameter
  
  AllocateAlleleCountArrays(options->getPopulations());

}

void DispersionFreqs::PrintPrior(const Vector_s& PopLabels, LogWriter& Log)const{
  if(IsHistoricAlleleFreq){//historic allele freq (dispersion) model
    Log << "Gamma prior on dispersion parameters with means and variances:\n";
    for( int k = 0; k < Populations; k++ ){
      Log << PopLabels[k] << ": "
	  << psi[k]/tau[k] << "  " << psi[k]/(tau[k]*tau[k]) << "\n";
    }
    Log << "\n";
  }
  else if( CorrelatedAlleleFreqs){//correlated allele freq model
    Log << "Gamma prior on allele frequency dispersion parameter with mean and variance:\n"
	<< psi[0]/tau[0] << "  " << psi[0]/(tau[0]*tau[0]) << "\n";
  }
}

void DispersionFreqs::LoadAlleleFreqs(AdmixOptions* const options, InputAdmixData* const data_, LogWriter &Log)
{
  int newrow;
  int row = 0;
  const Matrix_s* temporary = 0;

  //allocate frequency arrays
  Freqs.array = new double*[NumberOfCompositeLoci];
  AlleleFreqsMAP.array = Freqs.array;

  int offset = 0;
  bool oldformat = false;//flag for old format allelefreqfile
  bool file = false;//flag to indicate if priors have been supplied in a file

  //old format allelefreqfile
  if( strlen( options->getAlleleFreqFilename() ) ){
    temporary = &(data_->getAlleleFreqData());

    offset = 1;
    oldformat = true;
    file = true;
  }
  else if( strlen( options->getHistoricalAlleleFreqFilename() ) || strlen( options->getPriorAlleleFreqFilename() ) ){
    offset = 0;
    oldformat = false;
    file = true;
    //Historic DispersionFreqs
    if( strlen( options->getHistoricalAlleleFreqFilename() ) ){
      HistoricAlleleFreqs = new double*[NumberOfCompositeLoci];
      HistoricAlleleCounts = new double*[NumberOfCompositeLoci];
      if( options->getOutputFST() ){
	calculateFST = true;
	Fst = alloc2D_d(NumberOfCompositeLoci, Populations);
	SumFst = alloc2D_d(NumberOfCompositeLoci, Populations);
      }
      MuProposal = new std::vector<StepSizeTuner>[NumberOfCompositeLoci];
      temporary = &(data_->getHistoricalAlleleFreqData());
      IsHistoricAlleleFreq = true;
      
    } else {
      //Prior on DispersionFreqs
      temporary = &(data_->getPriorAlleleFreqData());
      IsHistoricAlleleFreq = false;
    }
  }

  bool useinitfile = false;
  if(RandomAlleleFreqs){

    if(hapmixmodel){
 
      file = false;
      if(strlen( options->getAlleleFreqFilename() )){
	  useinitfile=true;
	  LoadInitialAlleleFreqs(options->getAlleleFreqFilename(), Log);
      }
    }
    else
      PriorParams = new double*[NumberOfCompositeLoci];//2D array    "      "     otherwise
  }

  //set static members of CompositeLocus
  CompositeLocus::SetRandomAlleleFreqs(RandomAlleleFreqs);
  CompositeLocus::SetNumberOfPopulations(Populations);

if(!useinitfile)
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){

    Freqs.array[i] = new double[Loci->GetNumberOfStates(i)* Populations];

    if(file){//read allele freqs from file
	  newrow = row + (*Loci)(i)->GetNumberOfStates() - offset;
	  LoadAlleleFreqs( *temporary, i, row+1, oldformat);//row+1 is the first row for this locus (+1 for the header)
	  row = newrow;
    }
    else {  //set default Allele Freqs
      SetDefaultAlleleFreqs(i);
      if(!hapmixmodel && RandomAlleleFreqs){
//prior for hapmix model is set later in Initialise
	// reference prior on allele freqs: all elements of parameter vector set to 0.5
	// this is unrealistic for large haplotypes - should set all elements to sum to 1
	double defaultpriorparams = 0.5;
	SetDefaultPriorParams(i, defaultpriorparams);
      }
    }
  }
}

/**
 * Initialises the frequencies of each haplotype in the ith
 * composite locus, given Dirichlet priors (oldformat=false) or frequencies (oldformat=true)in matrix New.  
 * If oldformat=false, Allele freqs are set to their prior expectation. 
 * If fixed, allele freqs will be fixed at their prior expectations   
 *
 * New - a matrix containing either frequencies or
 *   parameters for the Dirichlet prior distribution of the allele frequencies. The second dimension is the allele number, 
 *   being in the range of zero to one less than the number of states.
 *   The sum of the prior parameters over all alleles in a population 
 *   (sumalpha) can be interpreted as 
 *   the "prior sample size". The first dimension is the population. Thus, for a 
 *   composite locus with four states and European and African 
 *   populations, the matrix might be:
 *
 *             Population
 *
 *            | EUR | AFR |
 *         ---|-----|-----|
 *          0 | 9.0 | 3.0 |
 *   State  1 | 3.0 | 4.0 |
 *          2 | 1.0 | 8.0 |
 *          3 | 2.0 | 1.0 |
 *
 * It is also permissable to have one column, in which case the parameters are (initially) the same across all
 * populations. This is required for a correlated allelefreqs model.
 * If Historic is true, sets "historical allele frequencies", where the model has been specified to allow the 
 * allele freqs in the admixed population 
 * to vary from the historical allele frequencies in the unadmixed ancestral populations that have 
 * been sampled. 
 */
void DispersionFreqs::LoadAlleleFreqs(const Matrix_s& New, int i, unsigned row0, bool ){
  // set size of allele freqs array for this locus
  // Freqs array has NumberOfStates - 1 elements for each population
  unsigned NumberOfStates = Loci->GetNumberOfStates(i);

    double sumalpha;
    // allele frequencies are initialised as expectations over the Dirichlet prior distribution, 
    // by dividing each prior parameter by the sum of the parameters.
    for( int j = 0; j < Populations; j++ ){
      //vector<double> NewCol = New.getCol(j);
      //sumalpha = accumulate( NewCol.begin(), NewCol.end(), 0.0, std::plus<double>() );
      sumalpha = 0.0;
      for( unsigned k = 0; k < NumberOfStates ; k++ )sumalpha+= convertValueFromFile(New[row0+k][j+1]);//New.get(row0+k, j+1);
      for( unsigned k = 0; k < NumberOfStates ; k++ )
	Freqs[i][ k + j*NumberOfStates ] = ( convertValueFromFile(New[row0+k][(!CorrelatedAlleleFreqs)*j+1] )
	  /*New.get( row0+k, (!CorrelatedAlleleFreqs)* j +1)*/ ) / sumalpha;
    }

  if(RandomAlleleFreqs){//no need to allocate remaining arrays if fixed allele freqs model
    PriorParams[i] = new double[NumberOfStates* Populations];
    if(IsHistoricAlleleFreq){
      HistoricAlleleFreqs[i] = new double[(NumberOfStates - 1)* Populations];
      fill(HistoricAlleleFreqs[i],HistoricAlleleFreqs[i] + (NumberOfStates - 1)* Populations, 0.0);
      HistoricAlleleCounts[i] = new double[NumberOfStates* Populations];
      
      for(unsigned row = 0; row < NumberOfStates; ++row)
	for(int col = 0; col < Populations; ++col){
	  double d = convertValueFromFile(New[row0+row][col+1]);
	  HistoricAlleleCounts[i][row*Populations +col] = d;//New.get(row0+row, col+1);
	  PriorParams[i][col*NumberOfStates + row] = d/*New.get(row0+row, col+1)*/ + 0.501; // why add 0.501? 
	}
      // set size of vector MuProposal
      if( NumberOfStates > 2 ){
	MuProposal[i].resize( Populations );
	for( int k = 0; k < Populations; k++ ){
	  MuProposal[i][k].SetParameters( 0.01, 0.001, 10.0, 0.23 );
	}
      }
    }
    else{ // priorallelefreqs model, with or without correlated allelefreqs
      RandomAlleleFreqs = true;
      for(unsigned row = 0; row < NumberOfStates; ++row)
	for(int col = 0; col < Populations; ++col){
	  PriorParams[i][col*NumberOfStates + row] = convertValueFromFile(New[row0+row][col+1]);//New.get(row0+row, col+1); 
	}
    }
}

}

/**
 * Given the number of ancestral populations, sets default values for
 * prior allele freq params.
 */
void DispersionFreqs::SetDefaultPriorParams(int i, double defaultpriorparams){
  int NumberOfStates = Loci->GetNumberOfStates(i);
  if(CorrelatedAlleleFreqs){
    PriorParams[i] = new double[NumberOfStates];
    // reference prior on allele freqs: all elements of parameter vector set to 0.5
    // this is unrealistic for large haplotypes - should set all elements to sum to 1
    fill(PriorParams[i], PriorParams[i] + NumberOfStates, 0.5);
  }
  else{
    PriorParams[i] = new double[NumberOfStates* Populations];
    fill(PriorParams[i], PriorParams[i] + NumberOfStates* Populations, defaultpriorparams);
  }
}

// ************************** Sampling and Updating *****************************************
/// samples allele frequencies and prior parameters.
void DispersionFreqs::Update(IndividualCollection*IC , bool afterBurnIn, double coolness){
  // Sample prior parameters
  
    if(IsHistoricAlleleFreq ){ // TODO: use class DirichletParamSampler 
      for( int i = 0; i < NumberOfCompositeLoci; i++ ){
	if( Loci->GetNumberOfStates(i) == 2 )
	  SampleDirichletParams1D( i);
	else
	  SampleDirichletParamsMultiDim( i);
      }
    }
    else if(CorrelatedAlleleFreqs){
      SampleDirichletParams();
    }
  
  // Sample allele frequencies conditional on Dirichlet priors 
  // then use AlleleProbs to set HapPairProbs in CompositeLocus
  // this is the only point at which SetHapPairProbs is called, apart from when 
  // the composite loci are initialized
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
   const unsigned NumberOfStates = Loci->GetNumberOfStates(i);
   //accumulate summary stats for update of priors
//    double sumlogfreqs1 = 0.0,sumlogfreqs2 = 0.0; 
//     for(int k = 0; k < Populations; ++k){
// 	sumlogfreqs1 += mylog(Freqs[i][k*NumberOfStates]);//allele 1
// 	sumlogfreqs2 += mylog(Freqs[i][k*NumberOfStates + 1]);//allele 2
//       }

   //Sample prior parameters
    if (FREQSAMPLER==FREQ_HAMILTONIAN_SAMPLER && (NumberOfStates > 2 ||
	accumulate(hetCounts[i], hetCounts[i]+Populations*Populations, 0, std::plus<int>()) > 0 )) {
      if(NumberOfStates==2) //shortcut for SNPs
	FreqSampler[i]->SampleSNPFreqs(Freqs[i], AlleleCounts[i], hetCounts[i], i, Populations, 
				       coolness);
      else FreqSampler[i]->SampleAlleleFreqs(Freqs[i], IC, i, NumberOfStates, Populations, 
					     coolness);
    }
    else //if (FREQSAMPLER==FREQ_CONJUGATE_SAMPLER)
      SampleAlleleFreqs(i, coolness);

    if(afterBurnIn)
      (*Loci)(i)->AccumulateAlleleProbs();

    //no need to update alleleprobs, they are the same as Freqs
    //set HapPair probs using updated alleleprobs
    (*Loci)(i)->SetHapPairProbs();

  }
  
  // Sample for allele frequency dispersion parameters, eta, conditional on allelefreqs using
  // Metropolis random-walk.
  if(  IsHistoricAlleleFreq){ 
    NumberOfEtaUpdates++;
    for( int k = 0; k < Populations; k++ ){
      SampleEtaWithRandomWalk(k, afterBurnIn);
    }
  }
  else if(CorrelatedAlleleFreqs){
    NumberOfEtaUpdates++;
    if(ETASAMPLER == 1) { // random walk metropolis
    SampleEtaWithRandomWalk(0, afterBurnIn);
    } // otherwise eta will be sampled in sampleDirichletParams
      
  }
  if( calculateFST && afterBurnIn && IsHistoricAlleleFreq ){
    UpdateFst();
  }
}

/** samples allele/hap freqs at i th composite locus as a conjugate Dirichlet update
 and stores result in array Freqs 
*/
void DispersionFreqs::SampleAlleleFreqs(int i, double coolness)
{
  unsigned NumStates = Loci->GetNumberOfStates(i);
  double* temp = new double[NumStates];
  double *freqs = new double[NumStates];
  
  int c = CorrelatedAlleleFreqs? 0 : 1; //indicates whether correlated allelefreqmodel
  //if there is, the Dirichlet params are common across populations
  for( int j = 0; j < Populations; j++ ){

    // to flatten likelihood when annealing, multiply realized allele counts by coolness
    for(unsigned s = 0; s < NumStates; ++s)
      temp[s] = PriorParams[i][c*j*NumStates + s] + coolness*AlleleCounts[i][s*Populations +j];

    Rand::gendirichlet(NumStates, temp, Freqs[i]+j*NumStates);
    for(unsigned s = 0; s < NumStates; ++s){
	if(Freqs[i][j*NumStates+s]==0.0) Freqs[i][j*NumStates+s] = 0.000001;
	if(Freqs[i][j*NumStates+s]==1.0) Freqs[i][j*NumStates+s] = 0.999999;
    }
    // sample HistoricAlleleFreqs as conjugate Dirichlet update with prior specified by PriorParams
    if( IsHistoricAlleleFreq ){
      for(unsigned s = 0; s < NumStates; ++s)
	temp[s] = PriorParams[i][c*j*NumStates + s] + HistoricAlleleCounts[i][s*Populations+j];
      Rand::gendirichlet(NumStates, temp, freqs);
      for(unsigned s = 0; s < NumStates-1; ++s)HistoricAlleleFreqs[i][s*Populations+j] = freqs[s];
    }
  }
  delete[] freqs;  
  delete[] temp;
}


/**
 * Used when model is specified to allow allele freqs in admixed
 * population to vary from the "historical" allele freqs in the
 * unadmixed population.  Dirichlet distribution for allele freqs at
 * locus with k alleles is specified with (k - 1) frequency parameters
 * (mu) and with a single dispersion parameter (eta) this method
 * samples mu and eta, and updates Dirichlet parameters for the allele
 * frequencies
 */
const double *DispersionFreqs::GetStatsForEta( int locus, int population)const
{
  int NumberOfStates = Loci->GetNumberOfStates(locus);
  double *stats = new double[ NumberOfStates ];
  fill(stats, stats+NumberOfStates, 0.0);
  if(IsHistoricAlleleFreq){
    // calculates sufficient stats for update of dispersion parameter
    double sumHistoric = 0.0;
    for( int i = 0; i < NumberOfStates - 1; i++ ){
      stats[ i ] = log( Freqs[locus][ i + population*NumberOfStates ] ) + log( HistoricAlleleFreqs[locus][ i*Populations + population ] );
      sumHistoric +=  HistoricAlleleFreqs[locus][ i*Populations + population ];
    }
    stats[ NumberOfStates - 1 ] = log( Freqs[locus][NumberOfStates-1 + population*NumberOfStates] ) 
      + log( 1 - sumHistoric );
  }
  else{//correlated allelefreqs model, sum over populations
    for(int k = 0; k < Populations; ++k){
      for( int i = 0; i < NumberOfStates; i++ ){
	stats[ i ] += log( Freqs[locus][ i + k*NumberOfStates ] );
      }
    }
  }

  return stats;
}

///sets PriorParams to sum to eta after sampling of eta
void DispersionFreqs::UpdatePriorAlleleFreqs(int j, const vector<vector<double> >& mu)
{
  for( int i = 0; i < NumberOfCompositeLoci; i++ ){
      for(int h = 0; h < Loci->GetNumberOfStates(i); ++h)
	PriorParams[i][j*Loci->GetNumberOfStates(i)+h] = mu[i][h];
  }
}

// problem here is to sample the Dirichlet parameters of a multinomial-Dirichlet distribution
// for a multi-allelic locus, we sample the Dirichlet proportion parameters conditional on the  
// dispersion parameter and the allele counts in admixed and historic populations
// by a Metropolis random walk
void DispersionFreqs::SampleDirichletParamsMultiDim( int locus)
{
  int NumberOfStates = Loci->GetNumberOfStates(locus);
  double LogAccProb;
  double* mu1 = new double[NumberOfStates];
  double* mu2 = new double[NumberOfStates];// mu1 is current vector of proportion parameters, mu2 is proposal

  for( int j = 0; j < Populations; j++ ){
    double Proposal1=0, Proposal2=0, f1=0, f2=0;

    for(int s = 0; s< NumberOfStates; ++s)
      mu1[s] = PriorParams[locus][j*NumberOfStates +s] / (eta[j]*MuProposal[locus][j].getStepSize());
    // propose mu2 from Dirichlet distribution with vector of expectations given by mu1
    // step size parameter controls variance: small step size gives small variance

    Rand::gendirichlet(NumberOfStates, mu1, mu2);
 
        
    for( int i = 0; i < NumberOfStates; i++ ){
      mu1[i] *= MuProposal[locus][j].getStepSize();
      // priors on proportion parameters are apparently Dirichlet(1.1, ,,, 1,1) 
      // f1 += 0.1 * log( mu1(i) ) + 0.1 * log( 1 - mu1(i) ); 
      //f2 += 0.1 * log( mu2(i) ) + 0.1 * log( 1 - mu2(i) );

      //using Di(1,...,1) prior
      Proposal1 += (eta[j] * mu2[i] - 1) * log( mu1[i] ) - gsl_sf_lngamma( eta[j] * mu2[i] );
      Proposal2 += (eta[j] * mu1[i] - 1) * log( mu2[i] ) - gsl_sf_lngamma( eta[j] * mu1[i] );
    }
      
    for( int k = 0; k < NumberOfStates; k++ ){
      f1 -= 2.0 * gsl_sf_lngamma( mu1[ k ]* eta[j] );
      f2 -= 2.0 * gsl_sf_lngamma( mu2[ k ]* eta[j] );
      f1 += gsl_sf_lngamma( mu1[ k ]* eta[j] + AlleleCounts[locus][k*Populations +j ]);
      f2 += gsl_sf_lngamma( mu2[ k ]* eta[j] + AlleleCounts[locus][k*Populations +j ]);
      f1 += gsl_sf_lngamma( mu1[ k ]* eta[j] + HistoricAlleleCounts[locus][k*Populations+j] );
      f2 += gsl_sf_lngamma( mu2[ k ]* eta[j] + HistoricAlleleCounts[locus][k*Populations+j] );
    }
    LogAccProb = f2-f1-Proposal2 + Proposal1;
    if(LogAccProb > 0.0) LogAccProb = 0.0;

    if( log(Rand::myrand()) < LogAccProb ){
      for(int s = 0; s < NumberOfStates; ++s)PriorParams[locus][j*NumberOfStates + s] = mu2[s]*eta[j];
    }
    MuProposal[locus][j].UpdateStepSize(exp(LogAccProb));
  }
  delete[] mu1;
  delete[] mu2;
}

void DispersionFreqs::SampleDirichletParams1D( int locus)
  // with a dispersion model, we sample PriorParams conditional on the observed counts 
  // (with the realized allele freqs integrated out) from a distribution that is 
  // proportional to the product of two binomial-beta likelihoods (for a diallelic locus)
  // we sample the proportion parameter of the beta distribution
  // conditional on the realized allele counts in the admixed population, 
  // the allele counts in the historic population, and the dispersion parameter
  // using an adaptive rejection sampler
  //Note: here NumberOfStates == 2  
{
  double lefttruncation = 0.1;//should be smaller
  MuSamplerArgs MuParameters;
  int* counts0 = new int[2 * Populations];
  double* counts1 = new double[2 * Populations];

  // Construct adaptive rejection sampler for mu.
  for(int i = 0; i < 2; ++i)//loop over the two states/alleles
    for(int j = 0; j < Populations; ++j) {
      counts0[i + j*2] = AlleleCounts[locus][i*Populations +j];
      counts1[i+ j*2] = HistoricAlleleCounts[locus][i*Populations +j];
    }

  AdaptiveRejection SampleMu;
  SampleMu.Initialise(true, true, 1.0, lefttruncation, fMu, dfMu );

  for( int j = 0; j < Populations; j++ ){

    MuParameters.eta = eta[j];
    MuParameters.K = j;
    MuParameters.counts = counts0;
    MuParameters.counts1 = counts1;

    SampleMu.setUpperBound( eta[j] - lefttruncation );//set upper limit for sampler
    //?? wrong  - should have upper limit of 1.0

    PriorParams[locus][ j*2 ] = SampleMu.Sample(&MuParameters, ddfMu); //first state/allele
    //?? wrong should be eta * SampleMu.Sample(...
    // Last (second) prior frequency parameter is determined by sum of mu's = eta.
    PriorParams[locus][ j*2 +1 ] = eta[j] - PriorParams[locus][ j*2 ];
  }
  delete[] counts0;
  delete[] counts1;
}

///sampling of Dirichlet parameters for allelefreqs in historic or correlated allele freq model
void DispersionFreqs::SampleDirichletParams() {
  if(IsHistoricAlleleFreq){
    for(int k = 0; k < Populations; ++k){
      for(int i = 0; i < NumberOfCompositeLoci; ++i){
	int NumberOfStates = Loci->GetNumberOfStates(i);
	int* counts = new int[2*NumberOfStates];
	for(int j = 0; j < NumberOfStates; ++j)
	  {
	    counts[j*NumberOfStates] = AlleleCounts[i][j*Populations +k];
	    counts[j*NumberOfStates+1] = (int)HistoricAlleleCounts[i][j*Populations+k];
	  }
	
	muSampler[i*Populations + k].Sample(PriorParams[i] + (k*NumberOfStates), eta[k], counts);

	//EtaSampler[k].addAlphas(i, PriorParams[i] + (k*NumberOfStates));
	//EtaSampler[k].addCounts(i, counts);
	delete[] counts;
      }

      //sample eta
      //eta[k] = EtaSampler[k].Sample();
      //SumEta[k]+=eta[k];
      //rescale PriorParams so they sum to eta
      //for( int i = 0; i < NumberOfCompositeLoci; i++ ){
      //double sum = accumulate(PriorParams[i][k*NumberOfStates], PriorParams[i][k*NumberOfStates]+NumberOfStates, 
      //			0.0, std::plus<double>());
      //sum = eta[k] / sum;
      //for(int t = 0; t < NumberOfStates; ++t){
      //  PriorParams[i][t] *= sum;
      //}
      //}

    }
  }
  
  else if(CorrelatedAlleleFreqs) {
    for(int i = 0; i < NumberOfCompositeLoci; ++i){
      muSampler[i].Sample(PriorParams[i], eta[0], AlleleCounts[i]);
    }

    if(ETASAMPLER == 2) { // hamiltonian sampler
      for(int i = 0; i < NumberOfCompositeLoci; ++i) {
	EtaSampler[0].addAlphas(i, PriorParams[i]);
	EtaSampler[0].addCounts(i, AlleleCounts[i]);
      }
      eta[0] = EtaSampler[0].Sample();
      SumEta[0] += eta[0];
      //rescale PriorParams so they sum to eta
      for( int i = 0; i < NumberOfCompositeLoci; ++i ) {
	int NumberOfStates = Loci->GetNumberOfStates(i);
	double sum = accumulate(PriorParams[i], PriorParams[i]+NumberOfStates, 0.0, 
				std::plus<double>());
	sum = eta[0] / sum;
	for(int t = 0; t < NumberOfStates; ++t) {
	  PriorParams[i][t] *= sum;
	}
      }
    }

  }
}

void DispersionFreqs::SampleEtaWithRandomWalk(int k, bool updateSumEta) {
  double etanew, LogPostRatio = 0.0, LogLikelihoodRatio = 0.0, LogPriorRatio = 0.0, AccProb = 0.0;
  double Denom = 0.0;
  double mineta = 0;
  vector< vector<double> > munew;
  // propose etanew from truncated log-normal distribution.
  do{
    etanew = exp( Rand::gennor( log( eta[k] ), etastep[k] ) );
  }while( etanew > 5000.0 );
  // Prior log-odds ratio (proposal ratio cancels with a part of the prior ratio)   
  LogPriorRatio = ( psi[k] - 1 ) * (log(etanew) - log(eta[k])) - tau[k] * ( etanew - eta[k] );
  // Log-likelihood ratio; numerator of integrating constant
  LogLikelihoodRatio += 2 * NumberOfCompositeLoci * ( gsl_sf_lngamma( etanew ) - gsl_sf_lngamma( eta[k] ) );
  for(int j = 0; j < NumberOfCompositeLoci; j++ ){
    std::vector<double> mu = GetPriorAlleleFreqs(j,k);
    vector<double>::const_iterator it = min_element(mu.begin(), mu.end());
    
    //mineta is a lower bound for proposal etanew
    if( mineta < 0.1 * eta[k] / *it )
      mineta = 0.1 * eta[k] / *it;
    
    munew.push_back( mu );
    //rescale munew so munew sums to etanew
    for( unsigned l = 0; l < mu.size(); l++ )
      munew[j][l] *= etanew / eta[k];
    
    const double *SumLogFreqs = GetStatsForEta(j,k);
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
  if(mineta < etanew )AccProb = xexp(LogPostRatio);
  
#ifdef DEBUGETA
  cout<< "eta["<<k<<"] = "<< eta[k]<<" eta* = "<<etanew<< " stepsize = "<<etastep[k]<<endl;
  cout << "LogPriorRatio= "<< LogPriorRatio <<" LogLikelihoodRatio= " << LogLikelihoodRatio <<endl
       << "Denom = "<< Denom<< " AccProb= exp("<<LogPostRatio<<") = "<< AccProb<<endl<<endl;
#endif
  
  // Acceptance test.
  if( log( Rand::myrand() ) < LogPostRatio && mineta < etanew ){
    eta[k] = etanew;
    UpdatePriorAlleleFreqs( k, munew );
  }
  
  if( !( NumberOfEtaUpdates % w ) ){
    etastep[k] = TuneEtaSampler[k].UpdateStepSize( AccProb );
  }
  
  if( updateSumEta )
    SumEta[k]+=eta[k];
}

// ******************** Output **************************
void DispersionFreqs::InitializeEtaOutputFile(const AdmixOptions* const options, const Vector_s& PopulationLabels, LogWriter &Log)
{
  Log.setDisplayMode(On);
  outputstream.open( options->getEtaOutputFilename(), ios::out );
  if( !outputstream )
    {
      Log <<  "ERROR: Couldn't open dispparamfile\n";
      exit( 1 );
    }
  else{
    Log << "Writing dispersion parameters to " << options->getEtaOutputFilename() << "\n";
    //Dispersion parameters (eta)
    if(IsHistoricAlleleFreq){
      for( int k = 0; k < Populations; k++ ){
	outputstream << "\"eta." << PopulationLabels[k]<< "\"\t"; 
      }
      outputstream << endl;
    }
    else if(CorrelatedAlleleFreqs){
      outputstream << "\"eta\"" <<endl; 
    }
  }
}

void DispersionFreqs::OutputErgodicAvg( int samples, std::ofstream *avgstream)const
{
  if( IsHistoricAlleleFreq ){
    for( int j = 0; j < Populations; j++ ){
      avgstream->width(9);
      *avgstream << setprecision(6) << SumEta[j] / samples << "\t";
    }
  }
  else if(CorrelatedAlleleFreqs){
    avgstream->width(9);
    *avgstream << setprecision(6) << SumEta[0] / samples << "\t";
  }
}

void DispersionFreqs::OutputEta(int iteration, const AdmixOptions *options, LogWriter &Log){
  if( IsHistoricAlleleFreq ){
    //output to logfile
    if( iteration == -1 )
      {
	Log.setDisplayMode(Off);
	Log.setPrecision(6);
	for( int j = 0; j < Populations; j++ ){
	  //Log.width(9);
	  Log << eta[j];
	}
      }
    //output to screen
    if( options->getDisplayLevel()>2 )
      {
	for( int j = 0; j < Populations; j++ ){
	  cout.width(9);
	  cout << setprecision(6) << eta[j] << " ";
	}
      }
    //Output to paramfile after BurnIn
    if( iteration > options->getBurnIn() ){
      outputstream << eta[0];
      for( int j = 1; j < Populations; j++ ){
	outputstream << "\t" << eta[j];
      }

      outputstream << endl;
    }
  }
  else if(CorrelatedAlleleFreqs){
    
    //output to logfile
    if( iteration == -1 )
      {
	Log.setDisplayMode(Off);
	Log.setPrecision(6);
	Log.width(9);
	Log << eta[0];
      }
    //output to screen
    if( options->getDisplayLevel()>2 )
      {
	cout.width(9);
	cout << setprecision(6) << eta[0] << " ";
      }
    //Output to paramfile after BurnIn
    if( iteration > options->getBurnIn() ){
      outputstream << eta[0]<< endl;
    }
  }
}

// *** FST functions ****************************
void DispersionFreqs::OpenFSTFile(const string& ResultsDir, LogWriter &Log){
  Log.setDisplayMode(On);
  const string FSTfilename = ResultsDir + "/" + FST_OUTPUT_FILE;
  Log << "Writing ergodic averages of FSTs to: " << FSTfilename << "\n";
  fstoutputstream.open( FSTfilename.c_str(), ios::out );
  if( !fstoutputstream ){
    Log << "ERROR: Couldn't open fstoutputfile\n";
    exit( 1 );
  }
}

void DispersionFreqs::UpdateFst()
{
  for(int locus = 0; locus < NumberOfCompositeLoci; ++locus){
    int NumberOfStates = Loci->GetNumberOfStates(locus);
    double q_admix,q_parental,f,H_admix, H_parental, H_combined, pbar;
    for( int k = 0; k < Populations; k++ ){
      H_admix = 0;
      H_parental = 0;
      H_combined = 0;
      q_admix=1.0;
      
      for( int i = 0; i < NumberOfStates - 1; i++ ){
	H_admix += Freqs[locus][ i + k*NumberOfStates ] * Freqs[locus][ i + k*NumberOfStates ];
	H_parental += HistoricAlleleFreqs[locus][ i*Populations+ k ] * HistoricAlleleFreqs[locus][ i*Populations+ k ];
	pbar = 0.5 * ( Freqs[locus][ i + k*NumberOfStates ] + HistoricAlleleFreqs[locus][ i*Populations+ k ] );
	H_combined += pbar * pbar;
	q_admix -= Freqs[locus][i + k*NumberOfStates];
      }
      
      H_admix += q_admix * q_admix;
      double sumHistoric = 0.0;
      for( int i = 0; i < NumberOfStates - 1; i++ ){
	sumHistoric +=  HistoricAlleleFreqs[locus][ i*Populations + k ];
      }
      q_parental = 1 - sumHistoric;
      H_parental += q_parental * q_parental;
      pbar = 0.5 * ( q_admix + q_parental );
      H_combined += pbar * pbar;
      
      H_combined = 1 - H_combined;
      H_admix = 1 - H_admix;
      H_parental = 1 - H_parental;
      f = ( H_combined - 0.5 * ( H_admix + H_parental ) ) / H_combined;
      Fst[locus][k] = 2*f / ( 1 + f );
    }
    for(int i=0;i < Populations;++i)
      SumFst[locus][i] += Fst[locus][i];
  }
}
void DispersionFreqs::OutputFST(){
  for( int j = 0; j < NumberOfCompositeLoci; j++ ){
    fstoutputstream << (*Loci)(j)->GetLabel(0);
    for(int k=0; k<Populations; ++k)fstoutputstream << " " << Fst[j][k];
    fstoutputstream << endl;
  }
}

// ********** log density and derivatives for Adaptive Rejection sampling of PriorParams  ***********

double fMu( double alpha, const void* const args )
{
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int pop = parameters->K;// number of populations
  double eta = parameters->eta;
  const int *counts0 = parameters->counts;
  const double *counts1 = parameters->counts1;

  //counts 0 are counts for admixed pop, counts1 for historic/unadmixed pop
  double f = 0.0, mu = alpha / eta;
  try{
    double logprior = 0.1 * log( mu ) + 0.1 * log( 1 - mu  );//Beta(1,1) prior
    f = logprior - 2 * lngamma( alpha ) - 2 * lngamma( eta - alpha );
    
    f += lngamma( alpha+counts0[pop*2] ) + lngamma( eta-alpha+counts0[1+pop*2] );//state 1 + state 2, admixed pop
    f += lngamma( alpha+counts1[pop*2] ) + lngamma( eta-alpha+counts1[1+pop*2] );//historic pop
  }
  catch (string s){
    throw string("Error in ALleleFreqs::fMu - " + s);
  }
  return f;
}

double dfMu( double alpha, const void* const args )
{
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int pop = parameters->K;// number of populations
  double eta = parameters->eta;
  const int *counts0 = parameters->counts;
  const double *counts1 = parameters->counts1;

  double x;
  double logprior = 0.1 / alpha - 0.1 / ( eta - alpha );//Beta(1,1) prior
  double f = logprior;
  x = eta - alpha;
  try{
    if(alpha < 0)throw string("arg mu to digamma is negative\n"); 
    if(x < 0)throw string("arg x to digamma is negative\n"); 
    f += 2 * ( digamma(x) - digamma(alpha) );
    
    //admixed pop
    //first state/allele
    x = alpha + counts0[pop*2];
    f += digamma(x);
    //second state/allele
    x = eta - alpha + counts0[1+pop*2];
    f -= digamma(x);
    
    //unadmixed pop
    x = alpha + counts1[pop*2];
    f += digamma(x);
    x = eta - alpha + counts1[1+pop*2];
    f -= digamma(x);
  }
  catch(string s){
    throw string("Error in dfMu in allelefreqs.cc - " + s);
  }
  return f;
}

double ddfMu( double alpha, const void* const args ) {
  // 
  const MuSamplerArgs* parameters = (const MuSamplerArgs*) args;
  int pop = parameters->K;// number of populations
  double eta = parameters->eta;
  const int *counts0 = parameters->counts;
  const double *counts1 = parameters->counts1;

  double x;
  double prior = -0.1 / (alpha*alpha) - 0.1 / (( eta - alpha ) * ( eta - alpha ) );
  double f = prior;
  try{
    x = eta - alpha;
    f -= 2 * ( trigamma(x) + trigamma(alpha) );
    
    x = alpha + counts0[pop*2];
    f += trigamma(x);
    x = eta - alpha + counts0[1+pop*2];
    f += trigamma(x);
    
    x = alpha + counts1[pop*2];
    f += trigamma(x);
    x = eta - alpha + counts1[1+pop*2];
    f += trigamma(x);
  }
  catch(string s){
    throw string("Error in ddfMu in allelefreqs.cc - " + s);
  }
  return f;
}

float DispersionFreqs::getEtaSamplerAcceptanceRate(int k)const{
  if(ETASAMPLER == 1) {
    return TuneEtaSampler[k].getExpectedAcceptanceRate();
  } else {
    return EtaSampler[k].getAcceptanceRate();
  }
}

float DispersionFreqs::getEtaSamplerStepsize(int k)const{
  if(ETASAMPLER==1) {
    return etastep[k];
  } else {
    return EtaSampler[k].getStepsize();
  }
}

float DispersionFreqs::getAlphaSamplerAcceptanceRate(int i)const{
  return muSampler[i].getAcceptanceRate();
}

float DispersionFreqs::getAlphaSamplerStepsize(int i)const{
  return muSampler[i].getStepsize();
}


